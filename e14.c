#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <stdlib.h>

mpz_t base_p3t;
mpz_t p3t;

size_t maxe = 0;

// Prints out a as long as it has less than max_digits
void safeprint(mpz_t a, size_t max_digits) {
    size_t alog10 = mpz_sizeinbase(a, 10) - 1;

    if (alog10 < max_digits) {
        gmp_printf("%Zd", a);
    } else {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mod_2exp(tmp, a, 8);
        size_t mod = mpz_get_ui(tmp);
        gmp_printf("# ≡ %ld mod 256", mod);
    }
}

/// Handles the values near the base cases
void single(mpz_t a, mpz_t b) {
    // gmp_printf("b = %Zd\n", b);

    assert(mpz_cmp_ui(b, sizeof(size_t) * 8) <= 0);
    size_t bi = mpz_get_ui(b);

    size_t amod2 = mpz_odd_p(a);

    switch (bi) {
        case 0: {
            if (amod2 == 0) {
                printf("HALT(");
                safeprint(a, 50);
                printf(", %ld)\n", bi);
                exit(0);
            }
            if (amod2 == 1) {
                mpz_div_2exp(b, a, 1);
                mpz_mul_ui(b, b, 6);
                mpz_sub_ui(b, b, 23);
                mpz_set_ui(a, 10);
            }
        } break;
        case 1: {
            if (amod2 == 0) {
                mpz_div_2exp(b, a, 1);
                mpz_mul_ui(b, b, 6);
                mpz_sub_ui(b, b, 23);
                mpz_set_ui(a, 10);
            }
            if (amod2 == 1) {
                printf("HALT(");
                safeprint(a, 50);
                printf(", %ld)\n", bi);
                exit(0);
            }
        } break;
        case 2: {
            printf("HALT(");
            safeprint(a, 50);
            printf(", %ld)\n", bi);
        } break;
        case 3: {
            if (amod2 == 0) {
                mpz_div_2exp(b, a, 1);
                mpz_mul_ui(b, b, 6);
                mpz_sub_ui(b, b, 20);
                mpz_set_ui(a, 10);
            }
            if (amod2 == 1) {
                mpz_div_2exp(b, a, 1);
                mpz_mul_ui(b, b, 6);
                mpz_sub_ui(b, b, 18);
                mpz_set_ui(a, 10);
            }
        } break;
        default: {
            // 2n -> 3n
            // 2n + 1 -> 3n + 1
            // 2n/2 * 3 = 3n
            // 2n+1/2 *3 = 
            mpz_mul_ui(a, a, 3);
            mpz_div_2exp(a, a, 1);
            mpz_sub_ui(b, b, 4);
            return;
        } break;
    }
    return;
}

/// Directly simulates t steps of the map
size_t direct(mpz_t a, size_t t) {

    size_t b = 0;
    mpz_t tmp; mpz_init(tmp);
    for (int i = 0; i < t; i++) {
        b += 4;

        mpz_div_2exp(tmp, a, 1); // a += (a >> 1);
        mpz_add(a, a, tmp);
    }
    mpz_clear(tmp);
    return b;
}
/// Accelerated simulation of 2^e steps
void accel_pow(mpz_t a, mpz_t b, size_t e) {
    // Good
    if (e < 8) {
        mpz_set(p3t, base_p3t);
        size_t b_2 = direct(a, 1 << e);
        mpz_sub_ui(b, b, b_2);
        return;
    }
    assert(e - 1 < sizeof(mp_bitcnt_t) * 8); // If you're running this algorithm on something larger than 2^64 it'll probably be too slow anyways.

    // Get last 2^{e-1} bits of a
    // Unless I'm stupid this should be equivalent
    mpz_t a_1;
    mpz_init(a_1);
    size_t t = 1 << (e - 1);
    mpz_mod_2exp(a_1, a, t); // a&m

    accel_pow(a_1, b, e - 1);
    if (e > maxe) {
        maxe = e;
        printf("Current e = %ld\r", maxe); fflush(stdout);
    }

    mpz_mul(p3t, p3t, p3t);
    mpz_mod_2exp(p3t, p3t, (mpz_get_ui(b) / 4) + 1);

    mpz_div_2exp(a, a, t); // a = a >> t
    mpz_mul(a, a, p3t); // a = (a >> t) * p3t
    mpz_add(a, a, a_1); // a = (a >> t) * p3t + a_1
    mpz_mod_2exp(a, a, (mpz_get_ui(b) / 4) + 1);

    mpz_mod_2exp(a_1, a, t); // a&m
    accel_pow(a_1, b, e - 1);

    mpz_mul(p3t, p3t, p3t);
    mpz_mod_2exp(p3t, p3t, (mpz_get_ui(b) / 4) + 1);

    mpz_div_2exp(a, a, t);
    mpz_mul(a, a, p3t);
    mpz_add(a, a, a_1);
    mpz_mod_2exp(a, a, (mpz_get_ui(b) / 4) + 1);

    mpz_clear(a_1);
}

/// Starts from the config and simulates until it hits the next reset (where b increases, or the machine halts)
/// This modifies the input config.
void next_reset(mpz_t a, mpz_t b) {
    maxe = 0;
    // t is the maximum amount of steps that direct() can be simulated without danger of b going negative.
    mpz_init(p3t);
    mpz_init(base_p3t);
    mpz_mul_2exp(base_p3t, base_p3t, 6);
    mpz_ui_pow_ui(base_p3t, 3, 1 << 6);
    mpz_t t;
    mpz_init_set(t, b);
    mpz_div_2exp(t, t, 2); // t = (b) / 4
    // gmp_printf("t = %Zd\n", t);

    while (mpz_cmp_ui(t, 0) > 0) { // loops until t == 0
        // e = floor(log_2(t))
        // e is the maximum exponent for which accel_pow can be run before we run into danger of hitting a negative b value.
        size_t e = mpz_sizeinbase(t, 2) - 1; // equivalent to floor(log2(t))
        accel_pow(a, b, e);

        mpz_set(t, b);
        mpz_div_2exp(t, t, 2); // t = (b - 1) / 2
        e = mpz_sizeinbase(t, 2) - 1;

        printf("  "); safeprint(a, 50);
        gmp_printf(" %Zd %Zd %ld\n", b, t, e);
    }
    assert(mpz_cmp_ui(b, 0) >= 0);
    while (mpz_cmp_ui(b, 5) <= 0) {
        single(a, b);
    }
}

#include <sys/time.h>

int sim(int init_a, int init_b) {
    mpz_t a, b;
    mpz_init_set_ui(a, init_a);
    mpz_init_set_ui(b, init_b);

    printf("START = A(%d, %d)\n", init_a, init_b);

    while (1) {
        struct timeval start, end;
        
        gettimeofday(&start, NULL);
        next_reset(a, b);
        gettimeofday(&end, NULL);

        gmp_printf("A(%Zd, ", a);
        safeprint(b, 50);
        printf(") - Computed in %lf seconds\n", ((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1000000.0);
    }
}

int main(int argc, char* argv[]) {
    mpz_t test_1;
    sim(10,1);

}