#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <gmp.h>
#include "zp_poly.h"
#include "fproot.h"

// Above this degree, the gcd and exact-division steps of the equal-degree
// splitting are delegated to the zp_poly library, whose half-gcd is
// sub-quadratic (ours is schoolbook and dominates for class polynomials with
// h(D) in the thousands).  The exponentiation stays on the Montgomery/
// Kronecker/Barrett path here, which is faster than zp_poly's at these sizes.
#define FPROOT_ZPGCD_MIN 1024

#define MAXLIMB 64          // supports s up to 31 (p up to ~1900 bits)

/* ===================================================================== *
 *  Montgomery F_p arithmetic (multi-limb, s padded so REDC covers        *
 *  Kronecker product coefficients < 2^34 * p^2 < p*R).                   *
 * ===================================================================== */

void fp_init (fp_ctx *C, const mpz_t p)
{
    C->n = (int) mpz_sizeinbase (p, 2);
    C->s = (C->n + 34 + 63) / 64;             // 64 s >= n + 34
    assert (2 * C->s + 2 <= MAXLIMB);
    mpz_init_set (C->pz, p);

    size_t cnt;
    C->p  = calloc (C->s, sizeof(mp_limb_t));
    mpz_export (C->p, &cnt, -1, sizeof(mp_limb_t), 0, 0, p);

    mp_limb_t p0 = C->p[0], inv = p0;         // p is odd
    for ( int k = 0 ; k < 6 ; k++ ) inv *= 2 - p0 * inv;   // -> p0^{-1} mod 2^64
    C->pinv = -inv;

    mpz_t t;  mpz_init (t);
    mpz_setbit (t, 64 * C->s);  mpz_mod (t, t, p);         // R mod p
    C->R1 = calloc (C->s, sizeof(mp_limb_t));
    mpz_export (C->R1, &cnt, -1, sizeof(mp_limb_t), 0, 0, t);
    mpz_mul (t, t, t);  mpz_mod (t, t, p);                 // R^2 mod p
    C->R2 = calloc (C->s, sizeof(mp_limb_t));
    mpz_export (C->R2, &cnt, -1, sizeof(mp_limb_t), 0, 0, t);
    mpz_clear (t);
}

void fp_clear (fp_ctx *C)
{ free (C->p);  free (C->R1);  free (C->R2);  mpz_clear (C->pz); }

// REDC: T (2s limbs, = a*b) -> r = T * R^{-1} mod p  (s limbs, in [0,p)).
static void redc (const fp_ctx *C, mp_limb_t *r, mp_limb_t *T)
{
    int s = C->s;  const mp_limb_t *p = C->p;
    mp_limb_t topcarry = 0;
    for ( int i = 0 ; i < s ; i++ ) {
        mp_limb_t m = T[i] * C->pinv;
        mp_limb_t c = mpn_addmul_1 (T + i, p, s, m);       // T[i..i+s-1] += m*p
        mp_size_t k = i + s;  mp_limb_t cc = c;
        while ( cc && k < 2*s ) { mp_limb_t old = T[k];  T[k] = old + cc;  cc = (T[k] < old);  k++; }
        topcarry += cc;                                    // carry out of limb 2s-1 (total < 2pR => <=1)
    }
    if ( topcarry || mpn_cmp (T + s, p, s) >= 0 ) mpn_sub_n (r, T + s, p, s);
    else mpn_copyi (r, T + s, s);
}

void fp_mul (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b)
{
    mp_limb_t T[MAXLIMB];
    mpn_mul_n (T, a, b, C->s);
    redc (C, r, T);
}

void fp_add (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b)
{
    mpn_add_n (r, a, b, C->s);                             // no carry-out: 2p < R (slack)
    if ( mpn_cmp (r, C->p, C->s) >= 0 ) mpn_sub_n (r, r, C->p, C->s);
}

void fp_sub (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b)
{
    if ( mpn_sub_n (r, a, b, C->s) ) mpn_add_n (r, r, C->p, C->s);   // borrow => add p
}

void fp_set_mpz (const fp_ctx *C, mp_limb_t *r, const mpz_t a)
{
    mp_limb_t t[MAXLIMB];  size_t cnt;
    memset (t, 0, C->s * sizeof(mp_limb_t));
    mpz_export (t, &cnt, -1, sizeof(mp_limb_t), 0, 0, a);  // a in [0,p) -> t (s limbs)
    fp_mul (C, r, t, C->R2);                               // t * R^2 * R^{-1} = t*R  (Montgomery)
}

void fp_get_mpz (const fp_ctx *C, mpz_t r, const mp_limb_t *a)
{
    mp_limb_t T[MAXLIMB];
    mpn_copyi (T, a, C->s);  memset (T + C->s, 0, C->s * sizeof(mp_limb_t));
    mp_limb_t out[MAXLIMB];
    redc (C, out, T);                                      // aR * R^{-1} = a
    mpz_import (r, C->s, -1, sizeof(mp_limb_t), 0, 0, out);
}

void fp_set_ui (const fp_ctx *C, mp_limb_t *r, unsigned long a)
{ mpz_t t; mpz_init_set_ui (t, a); mpz_mod (t, t, C->pz); fp_set_mpz (C, r, t); mpz_clear (t); }

int fp_is_zero (const fp_ctx *C, const mp_limb_t *a)
{ for ( int i = 0 ; i < C->s ; i++ ) if ( a[i] ) return 0;  return 1; }

int fp_equal (const fp_ctx *C, const mp_limb_t *a, const mp_limb_t *b)
{ return mpn_cmp (a, b, C->s) == 0; }

void fp_inv (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a)
{
    // r = a^(p-2) mod p, via square-and-multiply in Montgomery form
    mpz_t e;  mpz_init (e);  mpz_sub_ui (e, C->pz, 2);
    mp_limb_t acc[MAXLIMB], base[MAXLIMB];
    mpn_copyi (acc, C->R1, C->s);                          // 1 (Montgomery)
    mpn_copyi (base, a, C->s);
    for ( int i = 0 ; i < (int) mpz_sizeinbase (e, 2) ; i++ ) {
        if ( mpz_tstbit (e, i) ) fp_mul (C, acc, acc, base);
        fp_mul (C, base, base, base);
    }
    mpn_copyi (r, acc, C->s);
    mpz_clear (e);
}

/* ===================================================================== *
 *  Polynomials over F_p (flat arrays of s-limb Montgomery coefficients).*
 * ===================================================================== */

#define CF(f,i) ((f)->c + (size_t)(i) * s)          // pointer to coeff i (needs `s` in scope)

void fpoly_init (const fp_ctx *C, fp_poly *f, int cap)
{ f->c = malloc ((size_t) cap * C->s * sizeof(mp_limb_t));  f->cap = cap;  f->deg = -1; }

void fpoly_clear (fp_poly *f) { free (f->c);  f->c = NULL;  f->cap = 0;  f->deg = -1; }

static void fpoly_copy (const fp_ctx *C, fp_poly *dst, const fp_poly *src)
{
    int s = C->s;
    dst->deg = src->deg;
    if ( src->deg >= 0 ) mpn_copyi (dst->c, src->c, (size_t)(src->deg + 1) * s);
}

static void fpoly_trim (const fp_ctx *C, fp_poly *f)
{ int s = C->s;  while ( f->deg >= 0 && fp_is_zero (C, CF(f, f->deg)) ) f->deg--; }

static int ceil_log2 (int x) { int b = 0;  while ( (1 << b) < x ) b++;  return b; }

#define BARRETT_MIN 48          // use fast reduction for h of degree >= this

// g = a*b (full product), Montgomery coeffs, via Kronecker substitution.
static void fpoly_mul_kron (const fp_ctx *C, fp_poly *g, const fp_poly *a, const fp_poly *b)
{
    int s = C->s, da = a->deg, db = b->deg;
    if ( da < 0 || db < 0 ) { g->deg = -1; return; }
    int nterms = (da < db ? da : db) + 1;
    int w = (2 * C->n + ceil_log2 (nterms) + 1 + 63) / 64;  assert (w <= 2 * s);
    size_t alen = (size_t)(da + 1) * w, blen = (size_t)(db + 1) * w;
    mp_limb_t *pa = calloc (alen, sizeof(mp_limb_t)), *pb = calloc (blen, sizeof(mp_limb_t));
    for ( int i = 0 ; i <= da ; i++ ) mpn_copyi (pa + (size_t) i * w, CF(a, i), s);
    for ( int i = 0 ; i <= db ; i++ ) mpn_copyi (pb + (size_t) i * w, CF(b, i), s);
    mp_limb_t *rr = malloc ((alen + blen) * sizeof(mp_limb_t));
    if ( alen >= blen ) mpn_mul (rr, pa, alen, pb, blen);  else mpn_mul (rr, pb, blen, pa, alen);
    int gd = da + db;
    for ( int i = 0 ; i <= gd ; i++ ) {
        mp_limb_t T[MAXLIMB];  memset (T, 0, (size_t) 2 * s * sizeof(mp_limb_t));
        mpn_copyi (T, rr + (size_t) i * w, w);
        redc (C, CF(g, i), T);
    }
    g->deg = gd;  free (pa);  free (pb);  free (rr);  fpoly_trim (C, g);
}

// ---- task-parallel products (OpenMP): Karatsuba over polynomial halves ----
// The ~n squarings of the EDS exponentiation are sequentially dependent, so the
// parallelism has to live inside each product: split f = f0 + x^m f1 and compute
// the three Karatsuba half-products as OpenMP tasks (recursively: 3^depth leaves,
// each a serial Kronecker product big enough to keep GMP's FFT efficient).
#define FPROOT_PAR_MIN 2048     // parallelize when both halves have degree >= this

static int par_depth (int d)
{ return d >= 32768 ? 3 : d >= 8192 ? 2 : d >= 2*FPROOT_PAR_MIN ? 1 : 0; }

// read-only view of coefficients [lo, hi] of f (trimmed)
static fp_poly fpoly_view (const fp_ctx *C, const fp_poly *f, int lo, int hi)
{
    int s = C->s;
    fp_poly v;
    if ( hi > f->deg ) hi = f->deg;
    v.c = f->c + (size_t) lo * s;  v.deg = hi - lo;  v.cap = 0;
    while ( v.deg >= 0 && fp_is_zero (C, v.c + (size_t) v.deg * s) ) v.deg--;
    return v;
}

// g = a*b (full product; g distinct from a,b), parallel Karatsuba
static void fpoly_mul_par_rec (const fp_ctx *C, fp_poly *g, const fp_poly *a, const fp_poly *b, int depth)
{
    int s = C->s, da = a->deg, db = b->deg;
    if ( da < 0 || db < 0 ) { g->deg = -1; return; }
    int mind = da < db ? da : db;
    if ( depth <= 0 || mind < 2 * FPROOT_PAR_MIN ) { fpoly_mul_kron (C, g, a, b); return; }
    int m = ((da > db ? da : db) + 1) / 2;                   // split in coefficient index
    fp_poly a0 = fpoly_view (C, a, 0, m - 1), a1 = fpoly_view (C, a, m, da);
    fp_poly b0 = fpoly_view (C, b, 0, m - 1), b1 = fpoly_view (C, b, m, db);
    // asum = a0+a1, bsum = b0+b1
    int dasum = a0.deg > a1.deg ? a0.deg : a1.deg, dbsum = b0.deg > b1.deg ? b0.deg : b1.deg;
    fp_poly as, bs, z0, z1, z2;
    fpoly_init (C, &as, dasum + 2);  fpoly_init (C, &bs, dbsum + 2);
    for ( int i = 0 ; i <= dasum ; i++ ) {
        if ( i <= a0.deg && i <= a1.deg ) fp_add (C, as.c + (size_t)i*s, a0.c + (size_t)i*s, a1.c + (size_t)i*s);
        else if ( i <= a0.deg ) mpn_copyi (as.c + (size_t)i*s, a0.c + (size_t)i*s, s);
        else mpn_copyi (as.c + (size_t)i*s, a1.c + (size_t)i*s, s);
    }
    as.deg = dasum;  fpoly_trim (C, &as);
    for ( int i = 0 ; i <= dbsum ; i++ ) {
        if ( i <= b0.deg && i <= b1.deg ) fp_add (C, bs.c + (size_t)i*s, b0.c + (size_t)i*s, b1.c + (size_t)i*s);
        else if ( i <= b0.deg ) mpn_copyi (bs.c + (size_t)i*s, b0.c + (size_t)i*s, s);
        else mpn_copyi (bs.c + (size_t)i*s, b1.c + (size_t)i*s, s);
    }
    bs.deg = dbsum;  fpoly_trim (C, &bs);

    fpoly_init (C, &z0, a0.deg + b0.deg + 2);
    fpoly_init (C, &z1, as.deg + bs.deg + 2);
    fpoly_init (C, &z2, a1.deg + b1.deg + 2);
    #pragma omp task shared(z0, a0, b0)
    fpoly_mul_par_rec (C, &z0, &a0, &b0, depth - 1);
    #pragma omp task shared(z2, a1, b1)
    fpoly_mul_par_rec (C, &z2, &a1, &b1, depth - 1);
    fpoly_mul_par_rec (C, &z1, &as, &bs, depth - 1);
    #pragma omp taskwait

    // g = z0 + (z1 - z0 - z2) x^m + z2 x^{2m}
    int dg = da + db;
    for ( int i = 0 ; i <= dg ; i++ ) memset (g->c + (size_t)i*s, 0, (size_t)s*sizeof(mp_limb_t));
    for ( int i = 0 ; i <= z0.deg ; i++ ) mpn_copyi (g->c + (size_t)i*s, z0.c + (size_t)i*s, s);
    for ( int i = 0 ; i <= z2.deg ; i++ ) mpn_copyi (g->c + (size_t)(2*m + i)*s, z2.c + (size_t)i*s, s);
    for ( int i = 0 ; i <= z1.deg ; i++ ) {
        mp_limb_t t[MAXLIMB];
        mpn_copyi (t, z1.c + (size_t)i*s, s);
        if ( i <= z0.deg ) fp_sub (C, t, t, z0.c + (size_t)i*s);
        if ( i <= z2.deg ) fp_sub (C, t, t, z2.c + (size_t)i*s);
        fp_add (C, g->c + (size_t)(m + i)*s, g->c + (size_t)(m + i)*s, t);
    }
    g->deg = dg;  fpoly_trim (C, g);
    fpoly_clear (&as);  fpoly_clear (&bs);
    fpoly_clear (&z0);  fpoly_clear (&z1);  fpoly_clear (&z2);
}

// entry point: spins up the OpenMP team when the product is large enough
static void fpoly_mul_par (const fp_ctx *C, fp_poly *g, const fp_poly *a, const fp_poly *b)
{
    int mind = a->deg < b->deg ? a->deg : b->deg;
    int depth = par_depth (mind);
    if ( ! depth ) { fpoly_mul_kron (C, g, a, b); return; }
    #pragma omp parallel
    #pragma omp single nowait
    fpoly_mul_par_rec (C, g, a, b, depth);
}

// g = x^m f(1/x): reversal of f viewed as a degree-m polynomial.
static void fpoly_rev (const fp_ctx *C, fp_poly *g, const fp_poly *f, int m)
{
    int s = C->s;
    for ( int i = 0 ; i <= m ; i++ ) {
        int src = m - i;
        if ( src >= 0 && src <= f->deg ) mpn_copyi (CF(g, i), CF(f, src), s);
        else memset (CF(g, i), 0, (size_t) s * sizeof(mp_limb_t));
    }
    g->deg = m;  fpoly_trim (C, g);
}

static void fpoly_truncate (const fp_ctx *C, fp_poly *f, int k)
{ if ( f->deg >= k ) f->deg = k - 1;  fpoly_trim (C, f); }

// hi = h^{-1} mod x^k  (h->c[0] a unit), by Newton iteration.
static void fpoly_inv_series (const fp_ctx *C, fp_poly *hi, const fp_poly *h, int k)
{
    int s = C->s;
    mp_limb_t two[MAXLIMB], z[MAXLIMB];  fp_set_ui (C, two, 2);  memset (z, 0, (size_t) s * sizeof(mp_limb_t));
    fp_inv (C, CF(hi, 0), CF(h, 0));  hi->deg = 0;             // 1/h0
    fp_poly ht, t, u;  int cap = 2 * k + 4;
    fpoly_init (C, &ht, cap);  fpoly_init (C, &t, cap);  fpoly_init (C, &u, cap);
    for ( int prec = 1 ; prec < k ; ) {
        int np = 2 * prec;  if ( np > k ) np = k;
        fpoly_copy (C, &ht, h);  fpoly_truncate (C, &ht, np);
        fpoly_mul_kron (C, &t, &ht, hi);  fpoly_truncate (C, &t, np);      // t = h*hi mod x^np
        u.deg = (t.deg < 0) ? 0 : t.deg;                                   // u = 2 - t
        for ( int i = 0 ; i <= u.deg ; i++ ) {
            if ( i == 0 ) { if ( t.deg >= 0 ) fp_sub (C, CF(&u,0), two, CF(&t,0)); else mpn_copyi (CF(&u,0), two, s); }
            else fp_sub (C, CF(&u, i), z, CF(&t, i));
        }
        fpoly_trim (C, &u);
        fpoly_mul_kron (C, &t, hi, &u);  fpoly_truncate (C, &t, np);       // hi = hi*(2-t) mod x^np
        fpoly_copy (C, hi, &t);
        prec = np;
    }
    fpoly_clear (&ht);  fpoly_clear (&t);  fpoly_clear (&u);
}

// r = u mod h, h monic degree d, using hi = rev_d(h)^{-1} mod x^{>=deg u-d+1}.
static void fpoly_redbarrett (const fp_ctx *C, fp_poly *r, const fp_poly *u,
                              const fp_poly *h, const fp_poly *hi)
{
    int s = C->s, d = h->deg, m = u->deg;  (void) s;
    if ( m < d ) { fpoly_copy (C, r, u); return; }
    int kq = m - d;                                           // quotient degree
    // q* = (rev_m(u) * hi) mod x^{kq+1}: only the low kq+1 terms of each factor matter,
    // so truncate first -> a short product (degree <= 2 kq) instead of degree m+deg(hi).
    fp_poly ru, hit, qs, q, qh;
    fpoly_init (C, &ru, m + 2);  fpoly_init (C, &hit, hi->deg + 2);   // hit holds a full copy of hi
    fpoly_init (C, &qs, 2 * kq + 2);  fpoly_init (C, &q, kq + 2);  fpoly_init (C, &qh, m + 2);
    fpoly_rev (C, &ru, u, m);  fpoly_truncate (C, &ru, kq + 1);
    fpoly_copy (C, &hit, hi);  fpoly_truncate (C, &hit, kq + 1);
    fpoly_mul_par (C, &qs, &ru, &hit);  fpoly_truncate (C, &qs, kq + 1);   // q* (degree <= kq)
    fpoly_rev (C, &q, &qs, kq);                                // q = rev_{kq}(q*)
    fpoly_mul_par (C, &qh, &q, h);                             // q*h  (degree m)
    fpoly_copy (C, r, u);
    for ( int i = 0 ; i <= qh.deg ; i++ ) fp_sub (C, CF(r, i), CF(r, i), CF(&qh, i));   // r = u - q*h
    fpoly_trim (C, r);
    fpoly_clear (&ru);  fpoly_clear (&hit);  fpoly_clear (&qs);  fpoly_clear (&q);  fpoly_clear (&qh);
}

// g = f^2 mod h   (g != f;  h monic, degree d >= 1;  Montgomery coeffs).
// hi != NULL => Barrett reduction (fast); NULL => schoolbook.
static void fpoly_sqrmod (const fp_ctx *C, fp_poly *g, const fp_poly *f, const fp_poly *h, const fp_poly *hi)
{
    int s = C->s, d = h->deg, df = f->deg;
    if ( df < 0 ) { g->deg = -1; return; }
    if ( par_depth (df) ) {                                 // large: task-parallel Karatsuba square
        fpoly_mul_par (C, g, f, f);
        goto reduce;
    }
    // Kronecker slot width: 64 w >= 2n + ceil(log2(#terms)) + 1,  #terms <= df+1
    int bbits = 2 * C->n + ceil_log2 (df + 1) + 1;
    int w = (bbits + 63) / 64;
    assert (w <= 2 * s);
    // pack f -> integer with coeff i at limb offset i*w
    int flen = (df + 1) * w;
    mp_limb_t *pk = calloc ((size_t) flen, sizeof(mp_limb_t));
    for ( int i = 0 ; i <= df ; i++ ) mpn_copyi (pk + (size_t) i * w, CF(f, i), s);
    mp_limb_t *sq = malloc ((size_t) 2 * flen * sizeof(mp_limb_t));
    mpn_sqr (sq, pk, flen);                                  // F^2, 2*flen limbs
    // extract coeff i in [0,2df], REDC to Montgomery
    int gd = 2 * df;
    for ( int i = 0 ; i <= gd ; i++ ) {
        mp_limb_t T[MAXLIMB];
        memset (T, 0, (size_t) 2 * s * sizeof(mp_limb_t));
        mpn_copyi (T, sq + (size_t) i * w, w);              // c_i < 2^{64w} occupies w limbs
        redc (C, CF(g, i), T);                              // (f^2)_i in Montgomery
    }
    g->deg = gd;
    free (pk);  free (sq);
reduce:
    if ( g->deg < d ) { fpoly_trim (C, g);  return; }       // no reduction needed
    if ( hi ) {                                             // fast (Barrett) reduction
        fp_poly r;  fpoly_init (C, &r, g->deg + 2);
        fpoly_redbarrett (C, &r, g, h, hi);
        fpoly_copy (C, g, &r);  fpoly_clear (&r);
    } else {                                                // schoolbook: x^d ≡ -h'(x)
        for ( int i = g->deg ; i >= d ; i-- ) {
            mp_limb_t *fi = CF(g, i);
            if ( fp_is_zero (C, fi) ) continue;
            for ( int j = 0 ; j < d ; j++ ) {
                mp_limb_t t[MAXLIMB];
                fp_mul (C, t, fi, CF(h, j));
                fp_sub (C, CF(g, i - d + j), CF(g, i - d + j), t);
            }
        }
        g->deg = d - 1;
    }
    fpoly_trim (C, g);
}

// f <- (x+a)*f mod h, in place.  h monic degree d.  a in Montgomery form.
// Shift + scalar-mult + subtract lead*h -- no polynomial multiplication.
static void fpoly_mul_linear (const fp_ctx *C, fp_poly *f, const mp_limb_t *a, const fp_poly *h)
{
    int s = C->s, d = h->deg, df = f->deg;
    if ( df < 0 ) return;
    int has_lead = (df == d - 1);
    mp_limb_t lead[MAXLIMB];
    if ( has_lead ) mpn_copyi (lead, CF(f, d - 1), s);
    int newdeg = has_lead ? d - 1 : df + 1;
    for ( int i = newdeg ; i >= 0 ; i-- ) {                 // high->low: shift source f_{i-1} still original
        mp_limb_t acc[MAXLIMB], t[MAXLIMB];
        if ( i <= df ) fp_mul (C, acc, a, CF(f, i));  else memset (acc, 0, (size_t) s * sizeof(mp_limb_t));
        if ( i - 1 >= 0 && i - 1 <= df ) fp_add (C, acc, acc, CF(f, i - 1));
        if ( has_lead && i < d ) { fp_mul (C, t, lead, CF(h, i));  fp_sub (C, acc, acc, t); }
        mpn_copyi (CF(f, i), acc, s);
    }
    f->deg = newdeg;
    fpoly_trim (C, f);
}

// b = (x+a)^e mod h  via left-to-right binary exponentiation.  Precomputes the
// reciprocal of h once so every squaring uses the fast (Barrett) reduction.
static void fpoly_powmod_linear (const fp_ctx *C, fp_poly *b, const mp_limb_t *a,
                                 const mpz_t e, const fp_poly *h)
{
    int s = C->s, d = h->deg;
    int use_barrett = (d >= BARRETT_MIN);
    fp_poly hrev, hi;
    if ( use_barrett ) {
        fpoly_init (C, &hrev, d + 2);  fpoly_init (C, &hi, d + 1);
        fpoly_rev (C, &hrev, h, d);                         // rev_d(h), constant term 1
        fpoly_inv_series (C, &hi, &hrev, d - 1);            // rev_d(h)^{-1} mod x^{d-1}
    }
    mpn_copyi (CF(b, 0), C->R1, s);  b->deg = 0;            // b = 1
    fp_poly t;  fpoly_init (C, &t, 2 * d + 2);
    for ( int i = (int) mpz_sizeinbase (e, 2) - 1 ; i >= 0 ; i-- ) {
        fpoly_sqrmod (C, &t, b, h, use_barrett ? &hi : NULL);
        fpoly_copy (C, b, &t);
        if ( mpz_tstbit (e, i) ) fpoly_mul_linear (C, b, a, h);
    }
    fpoly_clear (&t);
    if ( use_barrett ) { fpoly_clear (&hrev);  fpoly_clear (&hi); }
}

// r = f mod g (over F_p); g nonzero.
static void fpoly_rem (const fp_ctx *C, fp_poly *r, const fp_poly *f, const fp_poly *g)
{
    int s = C->s, dg = g->deg;
    fpoly_copy (C, r, f);
    if ( dg == 0 ) { r->deg = -1; return; }
    mp_limb_t ginv[MAXLIMB];  fp_inv (C, ginv, CF(g, dg));
    for ( int i = r->deg ; i >= dg ; i-- ) {
        if ( fp_is_zero (C, CF(r, i)) ) continue;
        mp_limb_t factor[MAXLIMB];  fp_mul (C, factor, CF(r, i), ginv);
        for ( int j = 0 ; j <= dg ; j++ ) {
            mp_limb_t t[MAXLIMB];  fp_mul (C, t, factor, CF(g, j));
            fp_sub (C, CF(r, i - dg + j), CF(r, i - dg + j), t);
        }
    }
    r->deg = dg - 1;
    fpoly_trim (C, r);
}

// q = f / g  (exact; g | f).  h monic not required.
static void fpoly_divexact (const fp_ctx *C, fp_poly *q, const fp_poly *f, const fp_poly *g)
{
    int s = C->s, dg = g->deg, df = f->deg;
    fp_poly rem;  fpoly_init (C, &rem, df + 2);  fpoly_copy (C, &rem, f);
    q->deg = df - dg;
    mp_limb_t ginv[MAXLIMB];  fp_inv (C, ginv, CF(g, dg));
    for ( int i = df ; i >= dg ; i-- ) {
        mp_limb_t factor[MAXLIMB];  fp_mul (C, factor, CF(&rem, i), ginv);
        mpn_copyi (CF(q, i - dg), factor, s);
        for ( int j = 0 ; j <= dg ; j++ ) {
            mp_limb_t t[MAXLIMB];  fp_mul (C, t, factor, CF(g, j));
            fp_sub (C, CF(&rem, i - dg + j), CF(&rem, i - dg + j), t);
        }
    }
    fpoly_clear (&rem);
    fpoly_trim (C, q);
}

static void fpoly_make_monic (const fp_ctx *C, fp_poly *f)
{
    int s = C->s;
    if ( f->deg < 0 ) return;
    mp_limb_t linv[MAXLIMB];  fp_inv (C, linv, CF(f, f->deg));
    for ( int i = 0 ; i <= f->deg ; i++ ) fp_mul (C, CF(f, i), CF(f, i), linv);
}

// out = gcd(A,B), made monic.
static void fpoly_gcd (const fp_ctx *C, fp_poly *out, const fp_poly *A, const fp_poly *B)
{
    int cap = (A->deg > B->deg ? A->deg : B->deg) + 2;
    fp_poly a, b, r;  fpoly_init (C, &a, cap);  fpoly_init (C, &b, cap);  fpoly_init (C, &r, cap);
    fpoly_copy (C, &a, A);  fpoly_copy (C, &b, B);
    while ( b.deg >= 0 ) {
        fpoly_rem (C, &r, &a, &b);
        fpoly_copy (C, &a, &b);
        fpoly_copy (C, &b, &r);
    }
    if ( a.deg >= 0 ) fpoly_make_monic (C, &a);
    fpoly_copy (C, out, &a);
    fpoly_clear (&a);  fpoly_clear (&b);  fpoly_clear (&r);
}

// ---- zp_poly bridge: sub-quadratic gcd / division for large degrees ----
static mpz_t *fpoly_to_zp (const fp_ctx *C, const fp_poly *f, mpz_t pp)
{
    int s = C->s;
    mpz_t *z = zp_poly_alloc (f->deg, pp);
    for ( int i = 0 ; i <= f->deg ; i++ ) fp_get_mpz (C, z[i], CF(f, i));
    return z;
}

static void fpoly_from_zp (const fp_ctx *C, fp_poly *f, mpz_t *z, int d)
{
    int s = C->s;
    for ( int i = 0 ; i <= d ; i++ ) fp_set_mpz (C, CF(f, i), z[i]);
    f->deg = d;
    fpoly_trim (C, f);
}

// out = gcd(A,B) (monic); zp_poly half-gcd when both degrees are large
static void fpoly_gcd_fast (const fp_ctx *C, fp_poly *out, const fp_poly *A, const fp_poly *B)
{
    int mind = A->deg < B->deg ? A->deg : B->deg;
    if ( mind < FPROOT_ZPGCD_MIN || A->deg < 0 || B->deg < 0 ) { fpoly_gcd (C, out, A, B); return; }
    mpz_t pp;  mpz_init_set (pp, C->pz);
    mpz_t *az = fpoly_to_zp (C, A, pp), *bz = fpoly_to_zp (C, B, pp);
    mpz_t *rz = zp_poly_alloc (mind, pp);
    int dr;
    zp_poly_gcd (rz, &dr, az, A->deg, bz, B->deg, pp);      // returns a monic gcd
    fpoly_from_zp (C, out, rz, dr);
    zp_poly_free (az, A->deg);  zp_poly_free (bz, B->deg);  zp_poly_free (rz, mind);
    mpz_clear (pp);
}

// q = f / g (g | f); zp_poly fast division when the degrees are large
static void fpoly_divexact_fast (const fp_ctx *C, fp_poly *q, const fp_poly *f, const fp_poly *g)
{
    if ( g->deg < FPROOT_ZPGCD_MIN || f->deg - g->deg < FPROOT_ZPGCD_MIN )
        { fpoly_divexact (C, q, f, g); return; }
    mpz_t pp;  mpz_init_set (pp, C->pz);
    mpz_t *fz = fpoly_to_zp (C, f, pp), *gz = fpoly_to_zp (C, g, pp);
    mpz_t *qz = zp_poly_alloc (f->deg - g->deg, pp);
    zp_poly_div (qz, fz, f->deg, gz, g->deg, pp);
    fpoly_from_zp (C, q, qz, f->deg - g->deg);
    zp_poly_free (fz, f->deg);  zp_poly_free (gz, g->deg);  zp_poly_free (qz, f->deg - g->deg);
    mpz_clear (pp);
}

void fp_find_root (const fp_ctx *C, fp_poly *h, mp_limb_t *root, uint64_t seed)
{
    int s = C->s;
    mpz_t e;  mpz_init (e);  mpz_sub_ui (e, C->pz, 1);  mpz_fdiv_q_2exp (e, e, 1);   // (p-1)/2
    int cap = 2 * h->deg + 2;
    fp_poly b, g, q;  fpoly_init (C, &b, cap);  fpoly_init (C, &g, cap);  fpoly_init (C, &q, cap);
    uint64_t st = seed ? seed : 0x9e3779b97f4a7c15ULL;
    mpz_t az;  mpz_init (az);
    mp_limb_t zero[MAXLIMB];  memset (zero, 0, (size_t) s * sizeof(mp_limb_t));

    while ( h->deg > 1 ) {
        mp_limb_t a[MAXLIMB], am[MAXLIMB];
        for ( int i = 0 ; i < s ; i++ ) { st = st * 6364136223846793005ULL + 1442695040888963407ULL; a[i] = st; }
        mpz_import (az, s, -1, sizeof(mp_limb_t), 0, 0, a);  mpz_mod (az, az, C->pz);
        fp_set_mpz (C, am, az);
        fpoly_powmod_linear (C, &b, am, e, h);              // (x+a)^((p-1)/2) mod h
        fp_sub (C, CF(&b, 0), CF(&b, 0), C->R1);            // b -= 1
        fpoly_trim (C, &b);
        fpoly_gcd_fast (C, &g, h, &b);
        int dg = g.deg;
        if ( dg <= 0 || dg >= h->deg ) continue;           // trivial split; new a
        if ( 2 * dg <= h->deg ) fpoly_copy (C, h, &g);      // keep smaller side
        else { fpoly_divexact_fast (C, &q, h, &g);  fpoly_copy (C, h, &q); }
        fpoly_make_monic (C, h);
    }
    fp_sub (C, root, zero, CF(h, 0));                       // h = x - r (monic) => root = -h_0
    mpz_clear (e);  mpz_clear (az);
    fpoly_clear (&b);  fpoly_clear (&g);  fpoly_clear (&q);
}

void fp_pow (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mpz_t e)
{
    mp_limb_t acc[MAXLIMB], base[MAXLIMB];
    mpn_copyi (acc, C->R1, C->s);  mpn_copyi (base, a, C->s);
    for ( int i = 0 ; i < (int) mpz_sizeinbase (e, 2) ; i++ ) {
        if ( mpz_tstbit (e, i) ) fp_mul (C, acc, acc, base);
        fp_mul (C, base, base, base);
    }
    mpn_copyi (r, acc, C->s);
}

// g <- g/(x-r) by synthetic division (r a root of g); degree drops by one.
static void fpoly_div_root (const fp_ctx *C, fp_poly *g, const mp_limb_t *r)
{
    int s = C->s, d = g->deg;
    mp_limb_t carry[MAXLIMB];  mpn_copyi (carry, CF(g, d), s);   // q_{d-1} = g_d
    for ( int i = d - 1 ; i >= 0 ; i-- ) {
        mp_limb_t qi[MAXLIMB], t[MAXLIMB];
        mpn_copyi (qi, carry, s);
        fp_mul (C, t, r, qi);  fp_add (C, carry, CF(g, i), t);   // next carry = g_i + r*q_i (reads old g_i)
        mpn_copyi (CF(g, i), qi, s);                            // store q_i
    }
    g->deg = d - 1;
}

int fp_find_all_roots (const fp_ctx *C, const fp_poly *hin, mp_limb_t *roots, int maxr, uint64_t seed)
{
    int s = C->s;
    if ( hin->deg < 1 ) return 0;
    fp_poly h, xp, g;
    fpoly_init (C, &h, hin->deg + 2);  fpoly_copy (C, &h, hin);
    fpoly_make_monic (C, &h);                                  // powmod/gcd/EDS need h monic
    fpoly_init (C, &xp, hin->deg + 2);  fpoly_init (C, &g, hin->deg + 2);
    mp_limb_t zero[MAXLIMB];  memset (zero, 0, (size_t) s * sizeof(mp_limb_t));
    fpoly_powmod_linear (C, &xp, zero, C->pz, &h);              // x^p mod h
    if ( xp.deg < 0 ) { memset (CF(&xp, 0), 0, (size_t) s * sizeof(mp_limb_t));  xp.deg = 0; }
    if ( xp.deg < 1 ) { memset (CF(&xp, 1), 0, (size_t) s * sizeof(mp_limb_t));  xp.deg = 1; }
    fp_sub (C, CF(&xp, 1), CF(&xp, 1), C->R1);                  // xp -= x
    fpoly_trim (C, &xp);
    fpoly_gcd_fast (C, &g, &h, &xp);                            // completely-split part
    int nr = 0;
    while ( g.deg >= 1 && nr < maxr ) {
        mp_limb_t r[MAXLIMB];
        if ( g.deg == 1 ) {                                    // g1 x + g0 -> root -g0/g1
            mp_limb_t g1inv[MAXLIMB];  fp_inv (C, g1inv, CF(&g, 1));
            fp_mul (C, r, CF(&g, 0), g1inv);  fp_sub (C, r, zero, r);
            mpn_copyi (roots + (size_t) nr * s, r, s);  nr++;  break;
        }
        fp_poly gc;  fpoly_init (C, &gc, g.deg + 2);  fpoly_copy (C, &gc, &g);
        fp_find_root (C, &gc, r, seed + nr);
        mpn_copyi (roots + (size_t) nr * s, r, s);  nr++;
        fpoly_div_root (C, &g, r);
        fpoly_clear (&gc);
    }
    fpoly_clear (&h);  fpoly_clear (&xp);  fpoly_clear (&g);
    return nr;
}

#ifdef FP_TEST_MAIN
/* self-test: Montgomery ops vs GMP mpz over random inputs */
int main (int argc, char **argv)
{
    unsigned long bits = argc > 1 ? strtoul (argv[1], 0, 10) : 256;
    gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, 12345);
    mpz_t p;  mpz_init (p);  mpz_urandomb (p, rs, bits);  mpz_setbit (p, bits-1);  mpz_nextprime (p, p);

    fp_ctx C;  fp_init (&C, p);
    printf ("p %lu bits, s=%d limbs\n", (unsigned long) mpz_sizeinbase(p,2), C.s);

    mpz_t a, b, e1, e2, r;  mpz_inits (a, b, e1, e2, r, NULL);
    mp_limb_t am[MAXLIMB], bm[MAXLIMB], rm[MAXLIMB];
    int bad = 0;
    for ( int it = 0 ; it < 200000 ; it++ ) {
        mpz_urandomm (a, rs, p);  mpz_urandomm (b, rs, p);
        fp_set_mpz (&C, am, a);  fp_set_mpz (&C, bm, b);
        // mul
        fp_mul (&C, rm, am, bm);  fp_get_mpz (&C, r, rm);
        mpz_mul (e1, a, b);  mpz_mod (e1, e1, p);
        if ( mpz_cmp (r, e1) ) { gmp_printf ("MUL bad\n"); bad++; }
        // add
        fp_add (&C, rm, am, bm);  fp_get_mpz (&C, r, rm);
        mpz_add (e1, a, b);  mpz_mod (e1, e1, p);
        if ( mpz_cmp (r, e1) ) { gmp_printf ("ADD bad\n"); bad++; }
        // sub
        fp_sub (&C, rm, am, bm);  fp_get_mpz (&C, r, rm);
        mpz_sub (e1, a, b);  mpz_mod (e1, e1, p);
        if ( mpz_cmp (r, e1) ) { gmp_printf ("SUB bad\n"); bad++; }
        // inv (skip a=0)
        if ( mpz_sgn (a) ) {
            fp_inv (&C, rm, am);  fp_get_mpz (&C, r, rm);
            mpz_invert (e1, a, p);
            if ( mpz_cmp (r, e1) ) { gmp_printf ("INV bad\n"); bad++; }
        }
        if ( bad > 5 ) break;
    }
    printf ("Montgomery self-test: %s\n", bad ? "FAIL" : "OK (200k random mul/add/sub/inv vs GMP)");
    mpz_clears (a, b, e1, e2, r, NULL);  fp_clear (&C);  mpz_clear (p);  gmp_randclear (rs);
    return bad ? 1 : 0;
}
#endif
