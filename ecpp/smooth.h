#ifndef _SMOOTH_H_
#define _SMOOTH_H_

#include <stdint.h>
#include <stddef.h>
#include <gmp.h>

/*
    Batched n^4-smoothness testing for one-shot ECPP (component b).

    For a fixed target prime p we must decide, for many candidate curve orders
    N = p+1 -/+ t (each ~ the size of p), whether N has a y-smooth divisor m with
    L < m < L*r (r = least prime of m), where y = n^4, n = bitlength(p), and
    L = isqrt(p)+1+isqrt(4*isqrt(p)) = (p^{1/4}+1)^2.  A valid m exists iff the
    y-smooth PART of N exceeds L, so the core primitive is: given a batch of N_i,
    compute their y-smooth parts.

    Method (Bernstein, "How to find smooth parts of integers"): with the single
    big integer P = prod_{prime q <= y} q, a remainder tree yields r_i = P mod N_i
    for the whole batch by reducing P once at the tree root and descending.  Then
    the y-smooth part of N_i is gcd(N_i, r_i^(2^s) mod N_i) with 2^s >= bitlen(N_i)
    (deterministic and exact: for a prime q | N_i, q | r_i iff q | P iff q <= y).

    Building/storing P is the only step that scales with y (|P| ~ 1.44*y bits), so
    P is built once per y and cached to disk; the per-batch cost is one big
    reduction P mod (prod N_i) plus cheap descent, amortized over the whole batch.
*/

// ---- parallel segmented prime sieve ----
// All primes in (plo, y] into *out (malloc'd, ascending; includes 2 iff plo < 2).
// Returns the count.  nthreads<=0 => omp default.
uint64_t sieve_primes_range (uint64_t plo, uint64_t y, uint64_t **out, int nthreads);

// ---- prime product  P = prod_{prime q in (lo, y]} q  (lo = 0: all q <= y) ----
typedef struct {
    mpz_t P;              // the product (squarefree)
    uint64_t y;           // upper smoothness bound
    uint64_t lo;          // lower bound (exclusive); 0 for a full product
    uint64_t nprimes;     // number of primes in the product (for introspection)
} smooth_base;

// Build P = prod_{q<=y} q by a segmented sieve + parallel product tree.
// nthreads<=0 => use omp default.
void smooth_base_build (smooth_base *sb, uint64_t y, int nthreads);

// Build the SEGMENT product P = prod_{lo < q <= y} q (a progressive-smoothness
// ladder rung; segments for different bit-lengths are shared since they depend
// only on (lo, y]).
void smooth_base_build_range (smooth_base *sb, uint64_t lo, uint64_t y, int nthreads);

// Cache P to / from disk (raw GMP export, with a small header).  Return 1 on ok.
int  smooth_base_save (const smooth_base *sb, const char *path);
int  smooth_base_load (smooth_base *sb, const char *path);   // sets sb from file

void smooth_base_clear (smooth_base *sb);

// Cheap corruption check: a valid P (y>=47) is divisible by the primorial of the
// first 15 primes.  Returns 1 if P passes, 0 if it looks corrupt.
int  smooth_base_selfcheck (const smooth_base *sb);

// ---- batched smooth-part extraction ----
// Given N[0..n) (each > 1), set S[i] = the y-smooth part of N[i] (S[i] | N[i],
// S[i] holds exactly the prime-power factors of N[i] whose prime is <= y).
// S[] must be an array of n initialized mpz_t.  nthreads<=0 => omp default.
void smooth_parts (const smooth_base *sb, const mpz_t *N, size_t n, mpz_t *S,
                   int nthreads);

// ---- certificate helpers ----
// Certificate size bounds for prime p:
//   *L    = isqrt(p) + 1 + isqrt(4 isqrt p)          (lower bound, must have m>L)
//   *Hass = p + 1 + isqrt(4 p)                        (Hasse upper bound, m<=Hass)
//   *n    = bitlength(p);  *n2 = n^2;  *n4 = n^4      (smoothness bound = n4)
void cert_bounds (const mpz_t p, mpz_t L, mpz_t Hass,
                  unsigned long *n, uint64_t *n2, uint64_t *n4);

// Factor a y-smooth integer S (all prime factors <= y ~ 2^34) into distinct
// primes pr[] with exponents ex[]; returns the number of distinct primes (<=64).
// Uses trial division of tiny primes then Pollard-rho; S must really be smooth.
int  factor_smooth (const mpz_t S, uint64_t *pr, int *ex, int cap);

// Given the y-smooth part S of a candidate N with S>L, choose m | S with
// L < m < L*r (r = least prime of m), and list the distinct primes of m lying in
// (n2, n4] into qs[] (ascending).  Sets *m and *nq (count of qs).  Returns 1 on
// success, 0 if no valid m could be formed (rare edge cases).  Primes of S above
// n4 (possible when S came from a rounded-up cached prime product) are skipped.
int  build_m (mpz_t m, uint64_t *qs, int *nq,
              const mpz_t S, const mpz_t L, uint64_t n2, uint64_t n4);
int  build_m2 (mpz_t m, uint64_t *qs, int *nq,                 // odd_only=1: skip the
              const mpz_t S, const mpz_t L, uint64_t n2, uint64_t n4, int odd_only);  // prime 2 (see smooth.c)

// Near-miss top-up: S is the y'-smooth part of N for some y' < n4.  Search the
// cofactor N/S (whose prime factors all exceed y') for primes q <= n4 with
// bounded Pollard rho -- a factor q <= n4 ~ 2^34 costs ~sqrt(q) ~ 2^17 mulmods,
// so this is milliseconds per candidate -- and fold q^v_q(N) into S.  Stops
// when the remaining cofactor (whp) has no factor <= n4 or the iteration budget
// is exhausted.  Never overshoots: only primes <= n4 are multiplied in.
void smooth_topup (mpz_t S, const mpz_t N, uint64_t n4);

#endif
