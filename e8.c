#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <stdlib.h>

mpz_t base_p3t;
mpz_t p3t;

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
    // Added this just because it made it a bit easier to visualize.
    // a always gets overwritten so it's fine to alias it like this.
    #define k a

    assert(mpz_cmp_ui(b, sizeof(size_t) * 8) <= 0);
    size_t bi = mpz_get_ui(b);

    if (mpz_odd_p(a)) { // a = 1 mod 2
        mpz_div_2exp(k, a, 1);

        if (bi >= 1) {
            // -> A(3k + 4, b - 1)
            mpz_mul_ui(a, k, 3);
            mpz_add_ui(a, a, 4);
            mpz_sub_ui(b, b, 1);
            return;
        } else {
            // -> A(4, 6k + 10)
            mpz_mul_ui(b, k, 6);
            mpz_add_ui(b, b, 10);
            mpz_set_ui(a, 4);
            return;
        }
    } else { // a = 0 mod 2
        mpz_div_2exp(k, a, 1);

        switch (bi) {
            case 0: {
                printf("HALT\n");
                exit(0);
            } break;
            case 1: {
                // -> A(4, 6k + 8)
                mpz_mul_ui(b, k, 6);
                mpz_add_ui(b, b, 8);
                mpz_set_ui(a, 4);
                return;
            } break;
            default: {
                // -> A(3k + 3, 2)
                mpz_mul_ui(a, k, 3);
                mpz_add_ui(a, a, 3);
                mpz_sub_ui(b, b, 2);
                return;
            }
        }
    }
    #undef k
}

/// Directly simulates t steps of the map
size_t direct(mpz_t a, size_t t) {

    size_t b = 0;
    mpz_t tmp; mpz_init(tmp);
    for (int i = 0; i < t; i++) {
        size_t r = mpz_odd_p(a);
        b += 2 - r;

        mpz_div_2exp(tmp, a, 1); // a += (a >> 1) + 3;
        mpz_add(a, a, tmp);
        mpz_add_ui(a, a, 3);
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

    mpz_mul(p3t, p3t, p3t);
    mpz_mod_2exp(p3t, p3t, mpz_get_ui(b));

    mpz_div_2exp(a, a, t); // a = a >> t
    mpz_mul(a, a, p3t); // a = (a >> t) * p3t
    mpz_add(a, a, a_1); // a = (a >> t) * p3t + a_1
    mpz_mod_2exp(a, a, mpz_get_ui(b));

    mpz_mod_2exp(a_1, a, t); // a&m
    accel_pow(a_1, b, e - 1);

    mpz_mul(p3t, p3t, p3t);
    mpz_mod_2exp(p3t, p3t, mpz_get_ui(b));

    mpz_div_2exp(a, a, t);
    mpz_mul(a, a, p3t);
    mpz_add(a, a, a_1);
    mpz_mod_2exp(a, a, mpz_get_ui(b));

    mpz_clear(a_1);
}

/// Starts from the config and simulates until it hits the next reset (where b increases, or the machine halts)
/// This modifies the input config.
void next_reset(mpz_t a, mpz_t b) {
    // t is the maximum amount of steps that direct() can be simulated without danger of b going negative.
    mpz_init(p3t);
    mpz_init(base_p3t);
    mpz_mul_2exp(base_p3t, base_p3t, 6);
    mpz_ui_pow_ui(base_p3t, 3, 1 << 6);
    mpz_t t;
    mpz_init_set(t, b);
    mpz_sub_ui(t, t, 1);
    mpz_div_2exp(t, t, 1); // t = (b - 1) / 2

    while (mpz_cmp_ui(t, 0)) { // loops until t == 0
        // e = floor(log_2(t))
        // e is the maximum exponent for which accel_pow can be run before we run into danger of hitting a negative b value.
        size_t e = mpz_sizeinbase(t, 2) - 1; // equivalent to floor(log2(t))
        accel_pow(a, b, e);

        mpz_set(t, b);
        mpz_sub_ui(t, t, 1);
        mpz_div_2exp(t, t, 1); // t = (b - 1) / 2

        printf("  "); safeprint(a, 50);
        gmp_printf(" %Zd %Zd %ld\n", b, t, e);
    }
    assert(mpz_cmp_ui(b, 0) >= 0);
    while (mpz_cmp_ui(a, 4) != 0) {
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
    sim(4,4);
}