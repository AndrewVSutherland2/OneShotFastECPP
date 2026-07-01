#ifndef _CORNACCHIA_H_
#define _CORNACCHIA_H_

#include <gmp.h>

/*
    Cornacchia for the ECPP CM-discriminant search: with p fixed, solve the norm
    equation  4p = t^2 + d*v^2  for many small d = |D| (D = -d < 0 a discriminant,
    so d == 0 or 3 mod 4).  The roles are inverted relative to classpoly (which
    fixes D and searches over primes p), so everything that depends only on p is
    precomputed once in a context.

    This is the BASELINE: a per-call Tonelli-Shanks modular square root plus the
    Euclidean descent, using GMP mpz.  It is the correctness reference and the
    ground truth for benchmarking the multiplicative / mpn-level optimizations.

    A context is single-thread; use one per thread for parallel search.
*/

typedef struct {
    mpz_t p, p4, L;          // p, 4p, floor(sqrt(4p)) = floor(2 sqrt p)
    int pmod8;               // p mod 8
    // Tonelli-Shanks data (depends only on p):
    mpz_t Q;                 // odd part of p-1, p-1 = Q*2^S
    int   S;
    mpz_t z;                 // z = n^Q mod p for least non-residue n (p == 1 mod 4)
    mpz_t e;                 // (Q+1)/2
    mpz_t pp1q;              // (p+1)/4   (p == 3 mod 4 fast path)
    // scratch (kept in the context so solve() allocates nothing):
    mpz_t x0, a, b, r, c, t0, t1, R, cc, bb;
    mpz_t two_p;                                  // 2p
    mp_limb_t *da, *db, *dr, *dq;                 // mpn Euclidean-descent scratch
    mp_size_t dcap;
    // counters (for benchmarking / introspection)
    unsigned long n_sqrt, n_descent;
} cornacchia_ctx;

void cornacchia_init  (cornacchia_ctx *ctx, const mpz_t p);
void cornacchia_clear (cornacchia_ctx *ctx);

// out = sqrt(a) mod p for a in [0,p); returns 1 if a is a QR (out set), else 0.
int  cornacchia_sqrtmodp (cornacchia_ctx *ctx, mpz_t out, const mpz_t a);

// Solve 4p = t^2 + d*v^2 (d = |D| > 0, d == 0 or 3 mod 4).  Returns 1 and sets
// t>=0, v>=0 if solvable (with the CM parity t == d mod 2), else 0.
int  cornacchia_solve (cornacchia_ctx *ctx, unsigned long d, mpz_t t, mpz_t v);

// As cornacchia_solve, but with x0 = a square root of -d mod p supplied (computed
// multiplicatively from a factor base).  Caller must guarantee (-d/p) = 1.
int  cornacchia_solve_x0 (cornacchia_ctx *ctx, unsigned long d, const mpz_t x0, mpz_t t, mpz_t v);

#endif
