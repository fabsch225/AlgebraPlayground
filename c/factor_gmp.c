/* factor_gmp.c
 *
 * Factorization utility using GMP.
 * - Trial division by small primes
 * - Miller-Rabin primality test via mpz_probab_prime_p
 * - Pollard-Brent (Pollard-Rho variant) with random polynomial f(x)=x^2+c
 *
 * Compile:
 *   gcc -O2 -std=c11 factor_gmp.c -o factor_gmp -lgmp
 *
 * Usage:
 *   ./factor_gmp <decimal-or-hex-number>
 * Example:
 *   ./factor_gmp "2^255-19"
 *   ./factor_gmp 0x1234abcd...
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>   // for getpid()

/* -------- Parameters -------- */
#define SMALL_PRIME_LIMIT 1000000   /* generate small primes up to this for trial division */
#define MR_REPS 25                  /* Miller-Rabin reps for mpz_probab_prime_p */

/* -------- small primes sieve -------- */
static int *sieve_list = NULL;
static size_t sieve_count = 0;

static void make_sieve(uint32_t limit) {
    uint32_t n = limit;
    char *is_composite = calloc(n + 1, 1);
    size_t cap = 0;
    for (uint32_t i = 2; i <= n; ++i) {
        if (!is_composite[i]) {
            ++cap;
            if (i <= n / i) {
                for (uint32_t j = i * i; j <= n; j += i) is_composite[j] = 1;
            }
        }
    }
    sieve_list = malloc(sizeof(int) * cap);
    size_t idx = 0;
    for (uint32_t i = 2; i <= n; ++i) if (!is_composite[i]) sieve_list[idx++] = (int)i;
    sieve_count = idx;
    free(is_composite);
}

/* --------- utilities --------- */
static void rand_state_init(gmp_randstate_t rstate) {
    unsigned long seed;
    seed = (unsigned long)time(NULL) ^ (unsigned long)getpid();
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed);
}

/* ---------- Miller-Rabin primality check ---------- */
static int is_probable_prime(const mpz_t n) {
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 1;
    if (mpz_even_p(n)) return 0;
    int r = mpz_probab_prime_p(n, MR_REPS);
    return r != 0;
}

/* ---------- Pollard-Brent implementation ---------- */
static void pollard_brent(mpz_t factor, const mpz_t n, gmp_randstate_t rstate) {
    if (mpz_cmp_ui(n, 1) == 0) { mpz_set_ui(factor, 0); return; }
    if (mpz_even_p(n)) { mpz_set_ui(factor, 2); return; }

    mpz_t y, c, m, g, r, q, x, ys, tmp, absdiff;
    mpz_inits(y, c, m, g, r, q, x, ys, tmp, absdiff, NULL);
    mpz_set_ui(m, 128);
    mpz_set_ui(g, 1);
    mpz_set_ui(r, 1);

    mpz_urandomm(y, rstate, n); mpz_add_ui(y, y, 1);
    mpz_urandomm(c, rstate, n); mpz_add_ui(c, c, 1);

    while (mpz_cmp_ui(g, 1) == 0) {
        mpz_set(x, y);
        for (mpz_set_ui(tmp, 0); mpz_cmp(tmp, r) < 0; mpz_add_ui(tmp, tmp, 1)) {
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, n);
        }
        mpz_set_ui(q, 1);
        unsigned long k = 0;
        unsigned long rr = mpz_get_ui(r);
        while (k < rr && mpz_cmp_ui(g, 1) == 0) {
            unsigned long limit = ( (unsigned long) mpz_get_ui(m) < (rr - k) ) ? mpz_get_ui(m) : (rr - k);
            mpz_set(ys, y);
            for (unsigned long i = 0; i < limit; ++i) {
                mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, n);
                mpz_sub(absdiff, x, y); mpz_abs(absdiff, absdiff);
                mpz_mul(q, q, absdiff); mpz_mod(q, q, n);
            }
            mpz_gcd(g, q, n);
            k += limit;
        }
        if (mpz_cmp_ui(g, 1) > 0) break;
        mpz_mul_ui(r, r, 2);
    }

    if (mpz_cmp_ui(g, 0) == 0 || mpz_cmp_ui(g, 1) == 0 || mpz_cmp(g, n) == 0) {
        mpz_set_ui(factor, 0);
    } else {
        mpz_set(factor, g);
    }

    mpz_clears(y, c, m, g, r, q, x, ys, tmp, absdiff, NULL);
}

/* ---------- factor vector helpers ---------- */
typedef struct {
    mpz_t *list;
    size_t len;
    size_t cap;
} factor_vec;

static void factor_vec_init(factor_vec *v) {
    v->len = 0; v->cap = 8;
    v->list = malloc(sizeof(mpz_t) * v->cap);
    for (size_t i = 0; i < v->cap; ++i) mpz_init(v->list[i]);
}

static void factor_vec_push(factor_vec *v, const mpz_t x) {
    if (v->len == v->cap) {
        size_t newcap = v->cap * 2;
        v->list = realloc(v->list, sizeof(mpz_t) * newcap);
        for (size_t i = v->cap; i < newcap; ++i) mpz_init(v->list[i]);
        v->cap = newcap;
    }
    mpz_set(v->list[v->len++], x);
}

static void factor_vec_clear(factor_vec *v) {
    for (size_t i = 0; i < v->cap; ++i) mpz_clear(v->list[i]);
    free(v->list);
    v->list = NULL;
    v->len = v->cap = 0;
}

/* ---------- recursive factorization ---------- */
static void factor_recursive(factor_vec *out, mpz_t n, gmp_randstate_t rstate) {
    if (mpz_cmp_ui(n, 1) == 0) return;
    if (is_probable_prime(n)) {
        factor_vec_push(out, n);
        return;
    }

    for (size_t i = 0; i < sieve_count; ++i) {
        int p = sieve_list[i];
        if ((mpz_cmp_ui(n, p) < 0)) break;
        if (mpz_divisible_ui_p(n, p)) {
            while (mpz_divisible_ui_p(n, p)) {
                mpz_t pz; mpz_init_set_ui(pz, p);
                factor_vec_push(out, pz);
                mpz_clear(pz);
                mpz_tdiv_q_ui(n, n, p);
            }
            factor_recursive(out, n, rstate);
            return;
        }
    }

    if (mpz_cmp_ui(n, 1) == 0) return;
    if (is_probable_prime(n)) {
        factor_vec_push(out, n);
        return;
    }

    mpz_t factor;
    mpz_init(factor);
    int attempts = 0;
    while (attempts < 50) {
        pollard_brent(factor, n, rstate);
        if (mpz_cmp_ui(factor, 0) != 0 && mpz_cmp(factor, n) != 0) break;
        attempts++;
    }
    if (mpz_cmp_ui(factor, 0) == 0 || mpz_cmp(factor, n) == 0) {
        factor_vec_push(out, n);
        mpz_clear(factor);
        return;
    }

    mpz_t cofactor;
    mpz_init(cofactor);
    mpz_tdiv_q(cofactor, n, factor);

    factor_recursive(out, factor, rstate);
    factor_recursive(out, cofactor, rstate);

    mpz_clear(factor);
    mpz_clear(cofactor);
}

/* ---------- main ---------- */
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <integer (decimal or 0xhex)>\n", argv[0]);
        return 1;
    }

    make_sieve(SMALL_PRIME_LIMIT);

    mpz_t N;
    mpz_init(N);
    if (strncasecmp(argv[1], "0x", 2) == 0) {
        if (mpz_set_str(N, argv[1]+2, 16) != 0) { fprintf(stderr, "Bad hex input\n"); return 1; }
    } else if (strchr(argv[1], '^') != NULL) {
        unsigned long k = 0;
        if (sscanf(argv[1], "2^%lu-19", &k) == 1) {
            mpz_set_ui(N, 1);
            mpz_mul_2exp(N, N, k);
            mpz_sub_ui(N, N, 19);
        } else {
            if (mpz_set_str(N, argv[1], 10) != 0) { fprintf(stderr, "Bad input\n"); return 1; }
        }
    } else {
        if (mpz_set_str(N, argv[1], 10) != 0) { fprintf(stderr, "Bad input\n"); return 1; }
    }

    gmp_printf("Input N = %Zd\n", N);

    gmp_randstate_t rstate;
    rand_state_init(rstate);

    if (is_probable_prime(N)) {
        gmp_printf("N is probably prime (Miller-Rabin)\n");
        //mpz_clear(N);
        //gmp_randclear(rstate);
        //return 0;
    }

    factor_vec out;
    factor_vec_init(&out);

    mpz_t ncopy;
    mpz_init_set(ncopy, N);

    factor_recursive(&out, ncopy, rstate);

    printf("Factors found:\n");
    for (size_t i = 0; i < out.len; ++i) {
        gmp_printf("  %Zd\n", out.list[i]);
    }

    factor_vec_clear(&out);
    mpz_clear(ncopy);
    mpz_clear(N);
    gmp_randclear(rstate);

    free(sieve_list);
    return 0;
}
