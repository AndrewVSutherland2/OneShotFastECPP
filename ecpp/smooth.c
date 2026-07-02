#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <gmp.h>
#include "smooth.h"

/* ===================================================================== *
 *  Segmented sieve: produce all primes <= y into a uint64 array.        *
 * ===================================================================== */

// Simple sieve of odd primes in [3, lim]; returns count, fills *out (malloc'd).
static uint32_t *base_odd_primes (uint64_t lim, size_t *cnt)
{
    if ( lim < 3 ) { *cnt = 0; return malloc (sizeof(uint32_t)); }
    size_t half = (size_t)((lim - 1) / 2) + 1;             // index i <-> 2i+1
    uint8_t *comp = calloc (half, 1);                      // comp[i]=1 => 2i+1 composite
    for ( uint64_t i = 1 ; (2*i+1)*(2*i+1) <= lim ; i++ ) {
        if ( comp[i] ) continue;
        uint64_t q = 2*i + 1;
        for ( uint64_t m = q*q ; m <= lim ; m += 2*q ) comp[(m-1)/2] = 1;
    }
    size_t n = 0;
    for ( size_t i = 1 ; i < half ; i++ ) if ( ! comp[i] ) n++;
    uint32_t *pr = malloc (n * sizeof(uint32_t));
    size_t k = 0;
    for ( size_t i = 1 ; i < half ; i++ ) if ( ! comp[i] ) pr[k++] = (uint32_t)(2*i+1);
    free (comp);
    *cnt = n;
    return pr;
}

// Sieve odd numbers in [lo,hi) (lo,hi odd bounds; lo>=3), collecting primes into
// a caller-provided growing buffer.  Chunked for cache friendliness.
#define CHUNK_ODDS (1u<<18)
static void sieve_collect (uint64_t lo, uint64_t hi, const uint32_t *base, size_t nbase,
                           uint64_t **vec, size_t *vn, size_t *vcap)
{
    if ( lo < 3 ) lo = 3;
    lo |= 1;                                               // make odd
    uint8_t *bits = malloc (CHUNK_ODDS);                   // one byte per odd (simple, L2-resident)
    for ( uint64_t c0 = lo ; c0 < hi ; c0 += 2*CHUNK_ODDS ) {
        uint64_t c1 = c0 + 2*CHUNK_ODDS;  if ( c1 > hi ) c1 = hi;
        size_t nodd = (size_t)((c1 - c0 + 1) / 2);         // odds in [c0,c1)
        memset (bits, 0, nodd);
        for ( size_t bi = 0 ; bi < nbase ; bi++ ) {
            uint64_t q = base[bi];
            if ( q*q >= c1 ) break;
            // smallest odd multiple k*q >= max(c0, q*q), k odd
            uint64_t k = (c0 + q - 1) / q;  if ( k < q ) k = q;
            if ( ! (k & 1) ) k++;
            uint64_t m = k * q;
            for ( ; m < c1 ; m += 2*q ) bits[(m - c0)/2] = 1;
        }
        if ( *vn + nodd > *vcap ) {
            while ( *vn + nodd > *vcap ) *vcap = *vcap ? 2*(*vcap) : (1u<<20);
            *vec = realloc (*vec, *vcap * sizeof(uint64_t));
        }
        for ( size_t i = 0 ; i < nodd ; i++ )
            if ( ! bits[i] ) (*vec)[(*vn)++] = c0 + 2*i;
    }
    free (bits);
}

// All primes in (plo, y] into *out (malloc'd, ascending; incl. 2 iff plo < 2).
static uint64_t sieve_primes_range (uint64_t plo, uint64_t y, uint64_t **out, int nth)
{
    size_t nbase;
    uint64_t blim = (uint64_t) sqrtl ((long double) y) + 2;
    uint32_t *base = base_odd_primes (blim, &nbase);
    uint64_t start = plo < 3 ? 3 : plo + 1;                // first candidate > plo, >= 3

    if ( nth <= 0 ) nth = omp_get_max_threads ();
    // Split [start,y] into nth contiguous odd bands; each thread collects its band.
    uint64_t **bv = calloc (nth, sizeof(*bv));
    size_t *bn = calloc (nth, sizeof(*bn));
    #pragma omp parallel num_threads(nth)
    {
        int tid = omp_get_thread_num ();
        int T   = omp_get_num_threads ();
        uint64_t span = (y >= start) ? (y - start + 1) : 0;
        uint64_t lo = start + (uint64_t)((__uint128_t) span * tid / T);
        uint64_t hi = start + (uint64_t)((__uint128_t) span * (tid+1) / T);
        if ( tid == T-1 ) hi = y + 1;
        size_t cap = 0, n = 0;  uint64_t *v = NULL;
        if ( lo < hi && span ) sieve_collect (lo, hi, base, nbase, &v, &n, &cap);
        bv[tid] = v;  bn[tid] = n;
    }
    // concatenate: 2 (if in range), then bands in order
    uint64_t total = 1;
    for ( int t = 0 ; t < nth ; t++ ) total += bn[t];
    uint64_t *all = malloc (total * sizeof(uint64_t));
    size_t k = 0;
    if ( plo < 2 && y >= 2 ) all[k++] = 2;
    for ( int t = 0 ; t < nth ; t++ ) {
        if ( bn[t] ) { memcpy (all + k, bv[t], bn[t]*sizeof(uint64_t)); k += bn[t]; }
        free (bv[t]);
    }
    free (bv); free (bn); free (base);
    *out = all;
    return k;
}

/* ===================================================================== *
 *  Product tree: P = prod of an ascending prime array.                  *
 * ===================================================================== */

// serial product of pr[lo..hi)
static void prod_serial (mpz_t out, const uint64_t *pr, size_t lo, size_t hi)
{
    if ( hi <= lo ) { mpz_set_ui (out, 1); return; }
    if ( hi - lo <= 16 ) {
        mpz_set_ui (out, pr[lo]);
        for ( size_t i = lo+1 ; i < hi ; i++ ) mpz_mul_ui (out, out, pr[i]);
        return;
    }
    size_t mid = lo + (hi - lo)/2;
    mpz_t a, b;  mpz_init (a);  mpz_init (b);
    prod_serial (a, pr, lo, mid);
    prod_serial (b, pr, mid, hi);
    mpz_mul (out, a, b);
    mpz_clear (a);  mpz_clear (b);
}

// parallel product via OpenMP tasks (spawn only for ranges > cutoff)
static void prod_par (mpz_t out, const uint64_t *pr, size_t lo, size_t hi, size_t cutoff)
{
    if ( hi - lo <= cutoff ) { prod_serial (out, pr, lo, hi); return; }
    size_t mid = lo + (hi - lo)/2;
    mpz_t a, b;  mpz_init (a);  mpz_init (b);
    #pragma omp task shared(a)
    prod_par (a, pr, lo, mid, cutoff);
    prod_par (b, pr, mid, hi, cutoff);
    #pragma omp taskwait
    mpz_mul (out, a, b);
    mpz_clear (a);  mpz_clear (b);
}

void smooth_base_build_range (smooth_base *sb, uint64_t lo, uint64_t y, int nthreads)
{
    sb->y = y;  sb->lo = lo;
    uint64_t *pr;
    uint64_t np = sieve_primes_range (lo, y, &pr, nthreads);
    sb->nprimes = np;
    mpz_init (sb->P);
    if ( nthreads <= 0 ) nthreads = omp_get_max_threads ();
    size_t cutoff = np / (size_t)(nthreads * 8) + 1;       // ~8 tasks/thread
    if ( cutoff < 4096 ) cutoff = 4096;
    #pragma omp parallel num_threads(nthreads)
    {
        #pragma omp single
        prod_par (sb->P, pr, 0, np, cutoff);
    }
    free (pr);
}

void smooth_base_build (smooth_base *sb, uint64_t y, int nthreads)
{ smooth_base_build_range (sb, 0, y, nthreads); }

void smooth_base_clear (smooth_base *sb)
{
    mpz_clear (sb->P);
    sb->y = 0;  sb->nprimes = 0;
}

int smooth_base_selfcheck (const smooth_base *sb)
{
    if ( sb->lo < 2 ) {
        // full product: divisible by 2*3*...*47 (fits in 64 bits) when y>=47
        unsigned long primorial15 = 614889782588491410UL;
        if ( sb->y < 47 ) return 1;
        return mpz_divisible_ui_p (sb->P, primorial15) != 0;
    }
    // segment (lo, y]: the first prime above lo must divide P
    mpz_t q;  mpz_init_set_ui (q, sb->lo);
    mpz_nextprime (q, q);
    int ok = mpz_cmp_ui (q, sb->y) > 0 || mpz_divisible_p (sb->P, q);
    mpz_clear (q);
    return ok;
}

/* ===================================================================== *
 *  Disk cache of P.                                                     *
 * ===================================================================== */
#define SMOOTH_MAGIC   0x534d5031554c3121ULL  /* "SMP1..." v1: full products, no lo field */
#define SMOOTH_MAGIC2  0x534d5032554c3121ULL  /* "SMP2..." v2: header carries lo (segments) */

// Cache format: 4x uint64 header {magic, y, nprimes, nlimbs} then nlimbs native
// mp_limb_t.  NOT mpz_out_raw -- its 4-byte size field overflows for P > 2^31
// bytes (i.e. p >= ~384-bit).  Native limbs + a 64-bit count have no such cap.
// Reads/writes in 1 GiB chunks so a single fwrite/fread never exceeds SSIZE_MAX.
static int rw_chunks (void *buf, size_t nbytes, FILE *f, int write)
{
    char *p = buf;  size_t chunk = 1ULL << 30;
    while ( nbytes ) {
        size_t k = nbytes < chunk ? nbytes : chunk;
        size_t got = write ? fwrite (p, 1, k, f) : fread (p, 1, k, f);
        if ( got != k ) return 0;
        p += k;  nbytes -= k;
    }
    return 1;
}

int smooth_base_save (const smooth_base *sb, const char *path)
{
    FILE *f = fopen (path, "wb");
    if ( ! f ) return 0;
    uint64_t nl = mpz_size (sb->P);
    uint64_t hdr[5] = { SMOOTH_MAGIC2, sb->y, sb->nprimes, nl, sb->lo };
    if ( fwrite (hdr, sizeof(hdr), 1, f) != 1 ) { fclose (f); return 0; }
    int ok = rw_chunks ((void*) mpz_limbs_read (sb->P), nl * sizeof(mp_limb_t), f, 1);
    fclose (f);
    return ok;
}

int smooth_base_load (smooth_base *sb, const char *path)
{
    FILE *f = fopen (path, "rb");
    if ( ! f ) return 0;
    uint64_t hdr[4];
    if ( fread (hdr, sizeof(hdr), 1, f) != 1 || (hdr[0] != SMOOTH_MAGIC && hdr[0] != SMOOTH_MAGIC2) ) { fclose (f); return 0; }
    sb->y = hdr[1];  sb->nprimes = hdr[2];  sb->lo = 0;
    uint64_t nl = hdr[3];
    if ( hdr[0] == SMOOTH_MAGIC2 && fread (&sb->lo, sizeof(uint64_t), 1, f) != 1 ) { fclose (f); return 0; }
    mpz_init (sb->P);
    mp_limb_t *ld = mpz_limbs_write (sb->P, nl);
    int ok = rw_chunks (ld, nl * sizeof(mp_limb_t), f, 0);
    mpz_limbs_finish (sb->P, nl);
    fclose (f);
    return ok;
}

/* ===================================================================== *
 *  Remainder tree + smooth-part extraction.                             *
 * ===================================================================== */

// y-smooth part of N given r = P mod N: gcd(N, r^(2^s) mod N), 2^s >= bitlen(N).
static void smooth_part_from_rem (mpz_t S, const mpz_t N, const mpz_t r)
{
    if ( mpz_sgn (r) == 0 ) { mpz_set (S, N); return; }    // N | P  => fully squarefree-smooth
    unsigned long bl = mpz_sizeinbase (N, 2);
    int s = 0;  while ( (1UL << s) < bl ) s++;
    mpz_t t;  mpz_init_set (t, r);
    for ( int k = 0 ; k < s ; k++ ) { mpz_mul (t, t, t);  mpz_mod (t, t, N); }
    mpz_gcd (S, t, N);
    mpz_clear (t);
}

// rem = P mod X, parallelized for |P| >> |X|.  Splits P into G limb-aligned
// chunks P = sum_g P_g * 2^(g*chunk*bits); reduces each chunk mod X in parallel
// (zero-copy limb aliasing) and recombines with 2^(g*chunk*bits) mod X.  Same
// total work as one big division, spread over G cores.
static void reduce_big_mod (mpz_t rem, const mpz_t P, const mpz_t X, int nth)
{
    mp_size_t Ln = mpz_size (P);
    mp_size_t Xn = mpz_size (X);
    if ( nth <= 1 || Ln < Xn * 4 || Ln < 8 ) { mpz_mod (rem, P, X); return; }
    int G = nth;  if ( (mp_size_t) G > Ln ) G = (int) Ln;
    mp_srcptr Pd = mpz_limbs_read (P);
    mp_size_t chunk = (Ln + G - 1) / G;                     // limbs per chunk

    mpz_t step;  mpz_init (step);                           // 2^(chunk*bits) mod X
    mpz_set_ui (step, 1);  mpz_mul_2exp (step, step, (mp_bitcnt_t) chunk * GMP_NUMB_BITS);
    mpz_mod (step, step, X);
    mpz_t *pw = malloc (G * sizeof(mpz_t));                 // pw[g] = 2^(g*chunk*bits) mod X
    mpz_init_set_ui (pw[0], 1);
    for ( int g = 1 ; g < G ; g++ ) { mpz_init (pw[g]); mpz_mul (pw[g], pw[g-1], step); mpz_mod (pw[g], pw[g], X); }

    mpz_t *part = malloc (G * sizeof(mpz_t));
    #pragma omp parallel for num_threads(nth) schedule(static)
    for ( int g = 0 ; g < G ; g++ ) {
        mpz_init (part[g]);
        mp_size_t lo = (mp_size_t) g * chunk;
        if ( lo >= Ln ) continue;                           // part[g] = 0
        mp_size_t len = chunk;  if ( lo + len > Ln ) len = Ln - lo;
        while ( len > 0 && Pd[lo + len - 1] == 0 ) len--;   // normalize high limb
        if ( len == 0 ) continue;
        mpz_t cg;  mpz_roinit_n (cg, Pd + lo, len);          // read-only alias, no copy
        mpz_mod (part[g], cg, X);
        mpz_mul (part[g], part[g], pw[g]);
        mpz_mod (part[g], part[g], X);
    }
    mpz_set_ui (rem, 0);
    for ( int g = 0 ; g < G ; g++ ) mpz_add (rem, rem, part[g]);
    mpz_mod (rem, rem, X);
    for ( int g = 0 ; g < G ; g++ ) { mpz_clear (part[g]); mpz_clear (pw[g]); }
    free (part);  free (pw);  mpz_clear (step);
}

void smooth_parts (const smooth_base *sb, const mpz_t *N, size_t n, mpz_t *S, int nthreads)
{
    if ( n == 0 ) return;
    if ( nthreads <= 0 ) nthreads = omp_get_max_threads ();

    size_t sz = 1;  while ( sz < n ) sz <<= 1;             // pad to power of two
    int nlev = 0;  { size_t s2 = sz;  while ( s2 ) { nlev++; s2 >>= 1; } }  // = log2(sz)+1

    // Build product-tree levels.  lev[0]=leaves (copies of N, pad 1), lev[top]=root.
    mpz_t **lev = malloc (nlev * sizeof(*lev));
    size_t *cnt = malloc (nlev * sizeof(*cnt));
    cnt[0] = sz;
    lev[0] = malloc (sz * sizeof(mpz_t));
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for ( size_t i = 0 ; i < sz ; i++ ) {
        mpz_init (lev[0][i]);
        if ( i < n ) mpz_set (lev[0][i], N[i]);  else mpz_set_ui (lev[0][i], 1);
    }
    for ( int d = 1 ; d < nlev ; d++ ) {
        cnt[d] = cnt[d-1] / 2;
        lev[d] = malloc (cnt[d] * sizeof(mpz_t));
        #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 64)
        for ( size_t j = 0 ; j < cnt[d] ; j++ ) {
            mpz_init (lev[d][j]);
            mpz_mul (lev[d][j], lev[d-1][2*j], lev[d-1][2*j+1]);
        }
    }

    // Reduce P once at the root, then descend computing P mod (each subproduct).
    mpz_t *rem = malloc (cnt[nlev-1] * sizeof(mpz_t));     // top level: 1 node
    mpz_init (rem[0]);
    reduce_big_mod (rem[0], sb->P, lev[nlev-1][0], nthreads);   // <-- the one big reduction (parallel)
    for ( int d = nlev - 2 ; d >= 0 ; d-- ) {
        mpz_t *rc = malloc (cnt[d] * sizeof(mpz_t));
        #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 64)
        for ( size_t j = 0 ; j < cnt[d] ; j++ ) {
            mpz_init (rc[j]);
            mpz_mod (rc[j], rem[j >> 1], lev[d][j]);
        }
        for ( size_t j = 0 ; j < cnt[d+1] ; j++ ) mpz_clear (rem[j]);
        free (rem);
        rem = rc;
    }
    // rem[i] = P mod N[i] for i<n.  Extract smooth parts (parallel).
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 32)
    for ( size_t i = 0 ; i < n ; i++ )
        smooth_part_from_rem (S[i], N[i], rem[i]);

    for ( size_t i = 0 ; i < sz ; i++ ) mpz_clear (rem[i]);
    free (rem);
    for ( int d = 0 ; d < nlev ; d++ ) {
        for ( size_t j = 0 ; j < cnt[d] ; j++ ) mpz_clear (lev[d][j]);
        free (lev[d]);
    }
    free (lev);  free (cnt);
}

/* ===================================================================== *
 *  Certificate helpers.                                                 *
 * ===================================================================== */

void cert_bounds (const mpz_t p, mpz_t L, mpz_t Hass,
                  unsigned long *n, uint64_t *n2, uint64_t *n4)
{
    unsigned long nn = mpz_sizeinbase (p, 2);
    *n = nn;  *n2 = (uint64_t) nn * nn;  *n4 = (*n2) * (*n2);
    mpz_t sp, t;  mpz_init (sp);  mpz_init (t);
    mpz_sqrt (sp, p);                                      // isqrt(p)
    mpz_mul_2exp (t, sp, 2);  mpz_sqrt (t, t);             // isqrt(4 isqrt p)
    mpz_add (L, sp, t);  mpz_add_ui (L, L, 1);             // L = isqrt(p)+1+isqrt(4 isqrt p)
    mpz_mul_2exp (t, p, 2);  mpz_sqrt (t, t);              // isqrt(4p)
    mpz_add (Hass, p, t);  mpz_add_ui (Hass, Hass, 1);     // p+1+isqrt(4p)
    mpz_clear (sp);  mpz_clear (t);
}

// Brent's rho: put a nontrivial factor of n (composite) into d.
static void rho_factor (mpz_t d, const mpz_t n)
{
    if ( mpz_even_p (n) ) { mpz_set_ui (d, 2); return; }
    mpz_t x, y, c, t, q, g;
    mpz_init (x); mpz_init (y); mpz_init (c); mpz_init (t); mpz_init (q); mpz_init (g);
    unsigned long seed = 2;
    for (;;) {
        mpz_set_ui (c, seed);  mpz_set_ui (y, 2);  mpz_set_ui (q, 1);
        mpz_set_ui (g, 1);
        long r = 1, m = 128;
        mpz_t xs, ys;  mpz_init (xs);  mpz_init (ys);
        while ( mpz_cmp_ui (g, 1) == 0 ) {
            mpz_set (x, y);
            for ( long i = 0 ; i < r ; i++ ) { mpz_mul (t, y, y); mpz_add (t, t, c); mpz_mod (y, t, n); }
            long k = 0;
            while ( k < r && mpz_cmp_ui (g, 1) == 0 ) {
                mpz_set (ys, y);
                long lim = (r - k < m) ? (r - k) : m;
                for ( long i = 0 ; i < lim ; i++ ) {
                    mpz_mul (t, y, y); mpz_add (t, t, c); mpz_mod (y, t, n);
                    mpz_sub (t, x, y);  mpz_abs (t, t);
                    mpz_mul (q, q, t);  mpz_mod (q, q, n);
                }
                mpz_gcd (g, q, n);
                k += m;
            }
            r <<= 1;
        }
        if ( mpz_cmp (g, n) == 0 ) {                       // backtrack
            mpz_set (y, ys);
            do { mpz_mul (t, y, y); mpz_add (t, t, c); mpz_mod (y, t, n);
                 mpz_sub (t, x, y); mpz_abs (t, t); mpz_gcd (g, t, n); }
            while ( mpz_cmp_ui (g, 1) == 0 );
        }
        mpz_clear (xs);  mpz_clear (ys);
        if ( mpz_cmp (g, n) != 0 ) { mpz_set (d, g); break; }
        seed++;                                            // retry with a new polynomial
    }
    mpz_clear (x); mpz_clear (y); mpz_clear (c); mpz_clear (t); mpz_clear (q); mpz_clear (g);
}

// recursive factor of a smooth m into a list of prime factors (with multiplicity)
static void fac_rec (const mpz_t m, uint64_t *fac, int *nf, int cap)
{
    if ( mpz_cmp_ui (m, 1) == 0 ) return;
    if ( mpz_probab_prime_p (m, 25) ) {
        if ( *nf < cap ) fac[(*nf)++] = mpz_get_ui (m);
        return;
    }
    mpz_t d, q;  mpz_init (d);  mpz_init (q);
    rho_factor (d, m);
    mpz_divexact (q, m, d);
    fac_rec (d, fac, nf, cap);
    fac_rec (q, fac, nf, cap);
    mpz_clear (d);  mpz_clear (q);
}

static int cmp_u64 (const void *a, const void *b)
{ uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;  return (x>y) - (x<y); }

int factor_smooth (const mpz_t S, uint64_t *pr, int *ex, int cap)
{
    uint64_t fac[256];  int nf = 0;
    mpz_t m;  mpz_init_set (m, S);
    // strip small primes by trial division (fast path; keeps rho inputs prime-ish)
    for ( uint64_t q = 2 ; q < 65536 && mpz_cmp_ui (m, 1) > 0 ; q += (q==2?1:2) ) {
        if ( mpz_divisible_ui_p (m, q) ) {
            while ( mpz_divisible_ui_p (m, q) ) { mpz_divexact_ui (m, m, q); if ( nf < 256 ) fac[nf++] = q; }
        }
        if ( mpz_cmp_ui (m, (unsigned long)q*q) < 0 && mpz_cmp_ui (m,1) > 0 ) {   // remainder is prime
            if ( nf < 256 ) fac[nf++] = mpz_get_ui (m);
            mpz_set_ui (m, 1);
            break;
        }
    }
    if ( mpz_cmp_ui (m, 1) > 0 ) fac_rec (m, fac, &nf, 256);
    mpz_clear (m);
    qsort (fac, nf, sizeof(uint64_t), cmp_u64);
    int k = 0;
    for ( int i = 0 ; i < nf ; ) {
        int j = i;  while ( j < nf && fac[j] == fac[i] ) j++;
        if ( k < cap ) { pr[k] = fac[i];  ex[k] = j - i;  k++; }
        i = j;
    }
    return k;
}

// Brent rho with an iteration budget; returns 1 with a nontrivial factor in d,
// 0 if the budget ran out (whp no factor <= ~budget^2 exists).
static int rho_bounded (mpz_t d, const mpz_t n, long budget)
{
    if ( mpz_even_p (n) ) { mpz_set_ui (d, 2); return 1; }
    mpz_t x, y, c, t, q, g, ys;
    mpz_inits (x, y, c, t, q, g, ys, NULL);
    long spent = 0;  int found = 0;
    for ( unsigned long seed = 2 ; spent < budget && ! found ; seed++ ) {
        mpz_set_ui (c, seed);  mpz_set_ui (y, 2);  mpz_set_ui (q, 1);  mpz_set_ui (g, 1);
        long r = 1, mstep = 128;
        while ( mpz_cmp_ui (g, 1) == 0 && spent < budget ) {
            mpz_set (x, y);
            for ( long i = 0 ; i < r ; i++ ) { mpz_mul (t, y, y); mpz_add (t, t, c); mpz_mod (y, t, n); }
            spent += r;
            long k = 0;
            while ( k < r && mpz_cmp_ui (g, 1) == 0 ) {
                mpz_set (ys, y);
                long lim = (r - k < mstep) ? (r - k) : mstep;
                for ( long i = 0 ; i < lim ; i++ ) {
                    mpz_mul (t, y, y); mpz_add (t, t, c); mpz_mod (y, t, n);
                    mpz_sub (t, x, y);  mpz_abs (t, t);
                    mpz_mul (q, q, t);  mpz_mod (q, q, n);
                }
                spent += lim;
                mpz_gcd (g, q, n);
                k += mstep;
            }
            r <<= 1;
        }
        if ( mpz_cmp_ui (g, 1) > 0 && mpz_cmp (g, n) < 0 ) { mpz_set (d, g);  found = 1; }
        else if ( mpz_cmp (g, n) == 0 ) {                       // backtrack for the exact factor
            do { mpz_mul (t, ys, ys); mpz_add (t, t, c); mpz_mod (ys, t, n);
                 mpz_sub (t, x, ys); mpz_abs (t, t); mpz_gcd (g, t, n);  spent++; }
            while ( mpz_cmp_ui (g, 1) == 0 && spent < budget );
            if ( mpz_cmp_ui (g, 1) > 0 && mpz_cmp (g, n) < 0 ) { mpz_set (d, g);  found = 1; }
        }
    }
    mpz_clears (x, y, c, t, q, g, ys, NULL);
    return found;
}

void smooth_topup (mpz_t S, const mpz_t N, uint64_t n4)
{
    mpz_t c, d, q;  mpz_inits (c, d, q, NULL);
    mpz_divexact (c, N, S);                                   // cofactor: all prime factors > y'
    long budget = 3L << 18;                                   // total rho iterations; a <=2^34 factor
    long cap = 1L << 18;                                      // costs ~2^17 expected, so cap each try
    while ( mpz_cmp_ui (c, 1) > 0 && budget > 0 ) {
        if ( mpz_probab_prime_p (c, 15) ) {                   // prime cofactor: take it iff <= n4
            if ( mpz_cmp_ui (c, n4) <= 0 ) mpz_mul (S, S, c);
            break;
        }
        if ( ! rho_bounded (d, c, budget < cap ? budget : cap) ) break;   // whp no factor <= n4 remains
        budget -= 1L << 17;                                   // charge one factor's expected cost
        // d may be composite (product of small primes): peel its primes
        uint64_t pr[64];  int ex[64];
        if ( mpz_sizeinbase (d, 2) <= 80 ) {                  // small enough to factor fully
            int k = factor_smooth (d, pr, ex, 64);
            for ( int i = 0 ; i < k ; i++ )
                for ( int e = 0 ; e < ex[i] ; e++ ) {
                    if ( pr[i] <= n4 ) mpz_mul_ui (S, S, pr[i]);
                    mpz_divexact_ui (c, c, pr[i]);
                }
            // fold in any extra powers of the found primes remaining in c
            for ( int i = 0 ; i < k ; i++ )
                while ( mpz_divisible_ui_p (c, pr[i]) ) {
                    mpz_divexact_ui (c, c, pr[i]);
                    if ( pr[i] <= n4 ) mpz_mul_ui (S, S, pr[i]);
                }
        } else {                                              // big composite chunk: just remove it
            mpz_divexact (c, c, d);
        }
    }
    mpz_clears (c, d, q, NULL);
}

int build_m (mpz_t m, uint64_t *qs, int *nq,
             const mpz_t S, const mpz_t L, uint64_t n2, uint64_t n4)
{
    uint64_t pr[64];  int ex[64];
    int k = factor_smooth (S, pr, ex, 64);
    if ( k == 0 ) return 0;

    // Accumulate prime factors largest-first until the product first exceeds L.
    // Primes above n4 are skipped: S may have been extracted with a smoothness
    // bound above n4 (a rounded-up cached prime product); such primes cannot
    // appear in a valid certificate, so m must clear L without them.
    mpz_set_ui (m, 1);
    uint64_t rleast = 0;  int done = 0;
    for ( int j = k - 1 ; j >= 0 && ! done ; j-- ) {
        if ( pr[j] > n4 ) continue;
        for ( int e = 0 ; e < ex[j] && ! done ; e++ ) {
            mpz_mul_ui (m, m, pr[j]);  rleast = pr[j];
            if ( mpz_cmp (m, L) > 0 ) done = 1;
        }
    }
    if ( ! done ) return 0;                                // n4-smooth part <= L

    mpz_t Lr;  mpz_init (Lr);  mpz_mul_ui (Lr, L, rleast);
    int ok = mpz_cmp (m, Lr) < 0;                          // need m < L*r
    mpz_clear (Lr);
    if ( ! ok ) return 0;                                  // rare boundary case

    *nq = 0;
    for ( int j = 0 ; j < k ; j++ )
        if ( pr[j] > n2 && pr[j] < n4 && mpz_divisible_ui_p (m, pr[j]) )
            qs[(*nq)++] = pr[j];                           // ascending (pr[] ascending)
    return 1;
}
