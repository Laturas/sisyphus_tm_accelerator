#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <stdlib.h> // Only for the exit(0) function

mpz_t base_p3t;
mpz_t p3t;

size_t maxe = 0;
size_t goal_e = 0;
size_t exponent_cap = 0;

/// Prints out a as long as it has less than max_digits
/// Otherwise prints the value mod 256
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
        mpz_clear(tmp);
    }
}

void mod2e128_a(mpz_t a) {
    mpz_mod_2exp(a, a, 8);
}

/// Directly simulates t steps of the map
size_t direct(mpz_t a, unsigned long long t) {
    mpz_t tmp; mpz_init(tmp);
    for (int i = 0; i < t; i++) {
        mpz_div_2exp(tmp, a, 1); // a += (a >> 1);
        mpz_add(a, a, tmp);
    }
    mpz_clear(tmp);
    return t;
}
/// Accelerated simulation of 2^e steps
size_t accel_pow(mpz_t a, size_t e) {
    // Good
    if (e < 8) {
        mpz_set(p3t, base_p3t);
        return direct(a, (size_t)1 << e);
    }
    assert(e - 1 < sizeof(mp_bitcnt_t) * 8); // If you're running this algorithm on something larger than 2^64 it'll probably be too slow anyways.

    // Get last 2^{e-1} bits of a
    // Unless I'm stupid this should be equivalent
    mpz_t a_1;
    mpz_init(a_1);
    size_t t = (size_t)1 << (e - 1);
    mpz_mod_2exp(a_1, a, t); // a&m
    size_t processed_iters = accel_pow(a_1, e - 1);
    if (e > maxe) {
        maxe = e;
        printf("Current e = %ld/%ld\r", maxe, goal_e); fflush(stdout);
    }

    mpz_mul(p3t, p3t, p3t);
    // mpz_mod_2exp(p3t, p3t, exponent_cap);

    mpz_div_2exp(a, a, t); // a = a >> t
    mpz_mul(a, a, p3t); // a = (a >> t) * p3t
    mpz_add(a, a, a_1); // a = (a >> t) * p3t + a_1
    // mpz_mod_2exp(a, a, exponent_cap);

    mpz_mod_2exp(a_1, a, t); // a&m
    processed_iters += accel_pow(a_1, e - 1);

    mpz_mul(p3t, p3t, p3t);
    // mpz_mod_2exp(p3t, p3t, exponent_cap);

    mpz_div_2exp(a, a, t);
    mpz_mul(a, a, p3t);
    mpz_add(a, a, a_1);
    // mpz_mod_2exp(a, a, exponent_cap);

    mpz_clear(a_1);
    return processed_iters;
}

const int tab64[64] = {
    63,  0, 58,  1, 59, 47, 53,  2,
    60, 39, 48, 27, 54, 33, 42,  3,
    61, 51, 37, 40, 49, 18, 28, 20,
    55, 30, 34, 11, 43, 14, 22,  4,
    62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21,
    56, 45, 25, 31, 35, 16,  9, 12,
    44, 24, 15,  8, 23,  7,  6,  5};
int log2_64 (size_t value) {
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return tab64[((size_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}

/// Starts from the config and simulates until it hits the next reset (where b increases, or the machine halts)
/// This modifies the input config.
void next_reset(mpz_t a, unsigned long long iterations) {
    maxe = 0;
    mpz_init(p3t);
    mpz_init(base_p3t);
    mpz_mul_2exp(base_p3t, base_p3t, 6);
    mpz_ui_pow_ui(base_p3t, 3, 1 << 6);

    while (iterations > 4) { // loops until t == 0
        maxe = 0;
        // e = floor(log_2(t))
        // e is the maximum exponent for which accel_pow can be run before we run into danger of hitting a negative b value.
        size_t e = log2_64(iterations);
        goal_e = e;
        iterations -= accel_pow(a, e);

        e = log2_64(iterations);
    }
    direct(a, iterations);
    mpz_mod_2exp(a, a, 8);
}

// Linux specific behavior. You can get rid of this as long as you remove the timing information.
// There's something equivalent on Windows
#include <sys/time.h>
#include <string.h>

double process_memory() {
    FILE *fp;
    char line[256];
    long vm_rss = 0;

    fp = fopen("/proc/self/status", "r");
    if (fp == NULL) {
        perror("Error opening /proc/self/status");
        return -1;
    }

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            sscanf(line, "VmRSS: %ld kB", &vm_rss);
            break;
        }
    }

    fclose(fp);
    return vm_rss / 1024.0;
}

typedef struct {
    char initialized;
    /// Upper half
    mpz_t x;
    /// Lower half
    mpz_t y;
} hydra_chunk;

#define MAX_CHUNKS 255
#define SPLIT_2EXP_POW 19

hydra_chunk chunk_list[MAX_CHUNKS] = {0};
size_t chunk_elements = 0;
hydra_chunk p3t_chunk = {0};

/// Splits a given number into two halves and stores it into a hydra chunk
void push_chunk(mpz_t to_chunk) {
    assert(chunk_elements < MAX_CHUNKS);

    mpz_init(chunk_list[chunk_elements].x);
    mpz_div_2exp(chunk_list[chunk_elements].x, to_chunk, SPLIT_2EXP_POW);

    mpz_init(chunk_list[chunk_elements].y);
    mpz_mod_2exp(chunk_list[chunk_elements].y, to_chunk, SPLIT_2EXP_POW);

    chunk_elements++;
    return;
}

/// Gets the lower 2^m bits of the given chunk.
void get_lower_bits(mpz_t OUT_lower_bits, size_t chunk) {
    // Every time an element is added to the chunk, the elements need to be multiplied by p3t one more time.
    // chunk_elements - chunk_position = # of times it needs to be multiplied
    // 1 chunk element - pos 0 = 1 p3t multiplication
    size_t p3tmuls = chunk_elements - chunk;

    
}

void LSHYDRA(int init_a, unsigned long long iterations) {
    // Top level function. Should function like accel_pow but works with chunks instead of the entire values.

    /// feed represents the bits that will be passed to the lower accel_pow calls.
    mpz_t feed; mpz_init(feed);

    while (iterations > 0) {
        mpz_set_ui(feed, 0);

        for (int i = 0; i < chunk_elements; i++) {
            mpz_t lower_bits; mpz_init(lower_bits);

            get_lower_bits(lower_bits, i);
            mpz_add(feed, feed, lower_bits);

            mpz_clear(lower_bits);
        }

        iterations -= accel_pow(feed, (size_t)1 << SPLIT_2EXP_POW);
        push_chunk(feed);
    }
}

int HYDRA(int init_a, unsigned long long iterations) {
    mpz_t a;
    mpz_init_set_ui(a, init_a);

    printf("Calculating H^{%d}(%llu)\n", init_a, iterations);

    struct timeval start, end;
    
    gettimeofday(&start, NULL);
    next_reset(a, iterations);
    gettimeofday(&end, NULL);

    gmp_printf("  H^{%llu} (%d) = %Zd mod 2^128\n", iterations, init_a, a);
    printf("  Computed in %lf seconds using %.3lf MB\n", ((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1000000.0, process_memory());
}

int main(int argc, char* argv[]) {
    HYDRA(3, 1 << 19);

    mpz_t a;
    mpz_init_set_ui(a, 3);

    direct(a, 1 << 19);
    mpz_mod_2exp(a, a, 8);
    gmp_printf("%Zd mod 256\n", a);
}
