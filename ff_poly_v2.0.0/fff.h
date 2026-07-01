#ifndef _FFF_INCLUDE_
#define _FFF_INCLUDE_

/*
    Copyright 2016-2017 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>
#include "ff_poly/ff.h"
#include "ff_poly/ffpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    This is the recommended interface to the ff_poly library.
*/
/*
TODO for version 5.0, switch to mm.h/m64.h, but this interface should be stable (ish)

typedef struct {
    mm_t p;         // odd prime less than or equal to FF_MAX_P
    mm_t pinv;      // -1/p mod B, where B=2^MM_BITS
    mm_t R;         // B mod p (1 in Montgomery rep)
    mm_t R2;        // B^2 mod p
    mm_t R3;        // B^3 mod p
    mm_t g2;        // generator for 2-sylow subgroup of Fp* (also a non-residue)
    mm_t g2inv;     // 1/g2 mod p
    mm_t half;      // B/2 mod p (1/2 in Montgomery rep)
    mm_t third;     // B/3 mod p (1/3 in Montgomery rep)
    mm_t fourth;    // B/4 mod p (1/4 in Montgomery rep)
    mm_t fifth;     // B/5 mod p (1/5 in Montgomery rep)
    mm_t sixth;     // B/6 mod p (1/6 in Montgomery rep)
    mm_t seventh;   // B/7 mod p (1/7 in Montgomery rep)
    mm_t z3;        // primitive cube root of unity or 1 if p=2mod3
    int v2;         // 2-adic valuation of p-1
    int v3;         // 3-adic valuation of p-1
    mm_t m2;        // (p-1)/2^v2   // integer
    mm_t m3;        // (p-1)/3^v3   // integer
} fffctx_t[1];
*/

static inline void fff_init(long p) { ff_setup_ui(p); }
static inline ff_t fff_char(void) { return _ff_p; }
static inline ff_t fff_p(void) { return _ff_p; }
static inline ff_t fff_zero(void) { return _ff_add_identity; }
static inline ff_t fff_one(void) { return _ff_mult_identity; }
static inline ff_t fff_one_half(void) { return _ff_half; }
static inline ff_t fff_one_third(void) { return _ff_third; }
static inline ff_t fff_one_fourth(void) { return _ff_fourth; }
static inline ff_t fff_negone(void) { return _ff_negone; }
static inline ff_t fff_random(void) { ff_t x; ff_random(&x); return x; }
static inline ff_t *fff_zero_vector(ff_t *v, int n) { memset(v,0,n*sizeof(*v)); return v; }

static inline unsigned long fff_to_ui (ff_t a) { return _ff_get_ui(a); }
static inline long fff_to_si (ff_t a) { return _ff_get_ui(a); }
static inline long fff_to_int (ff_t a) { return fff_to_si(a); }
static inline ff_t fff_from_ui (unsigned long n) { ff_t x; _ff_set_ui(x,n); return x; }
static inline ff_t fff_from_si (long n) { ff_t x; _ff_set_i(x,n); return x; }
static inline ff_t fff_from_int (long n) { return fff_from_si(n); }
static inline ff_t fff_from_mpz (mpz_t n) { ff_t x; _ff_set_mpz(x,n); return x; }
static inline ff_t fff_from_str (char *s) { mpz_t X; mpz_init_set_str(X,s,0); ff_t x = fff_from_mpz(X); mpz_clear(X); return x; }
static inline void fff_to_mpz (mpz_t A, ff_t a) { mpz_set_ui (A, fff_to_ui(a)); }

static inline ff_t fff_is_zero(ff_t a) { return _ff_zero(a); }
static inline ff_t fff_is_nonzero(ff_t a) { return ! _ff_zero(a); }
static inline ff_t fff_is_one(ff_t a) { return _ff_one(a); }

static inline ff_t fff_neg(ff_t a) { ff_t b; _ff_neg(b,a); return b; }
static inline ff_t fff_inc(ff_t a) { _ff_inc(a); return a; }
static inline ff_t fff_dec(ff_t a) { _ff_dec(a); return a; }
static inline ff_t fff_add(ff_t a, ff_t b) { ff_t c; _ff_add(c,a,b); return c; }
static inline ff_t fff_sub(ff_t a, ff_t b) { ff_t c; _ff_sub(c,a,b); return c; }
static inline ff_t fff_dbl(ff_t a) { return fff_add(a,a); }
static inline ff_t fff_mul3(ff_t a) { a += a+a; while ( a >= _ff_p ) a -= _ff_p; return a; }
static inline ff_t fff_div2(ff_t a) { ff_t c; _ff_halve(c,a); return c; }
static inline ff_t fff_half(ff_t a) { return fff_div2(a); }
static inline ff_t fff_mul(ff_t a, ff_t b) { ff_t c; _ff_mult(c,a,b); return c; }
static inline ff_t fff_addmul(ff_t a, ff_t b,ff_t c) { return fff_add(a,fff_mul(b,c)); } // OPTIMIZE THIS
static inline ff_t fff_sqr(ff_t a) { ff_t c; _ff_square(c,a); return c; }
static inline ff_t fff_cube(ff_t a) { return fff_mul(a,fff_sqr(a)); }
static inline ff_t fff_inv(ff_t a) { ff_t c; _ff_invert(c,a); return c; }
static inline ff_t fff_exp(ff_t a, uint64_t e) { ff_t x[1]; x[0] = a; ff_exp_ui(x,x,e); return x[0]; }
static inline ff_t fff_fac(long n) { if ( n >= fff_p() ) { return fff_zero(); } ff_t x = fff_one(), a = fff_one(); for ( int i = 2 ; i <= n ; i++ ) { a=fff_inc(a); x = fff_mul(x,a); } return x; }
static inline int fff_inv_array(ff_t b[], ff_t a[], int n) { return ff_parallel_invert_check(b,a,n); }
static inline int fff_is_square(ff_t a) { return ff_residue(a); }
static inline int fff_legendre(ff_t a) { return a ? (ff_residue(a)?1:-1) : 0; }
static inline int fff_sqrt(ff_t *c,ff_t a) { return ff_sqrt(c,&a); }
static inline int fff_invsqrt(ff_t *c,ff_t a) { return ff_invsqrt(c,&a,0); }
static inline int fff_cbrt(ff_t *c,ff_t a) { return ff_cbrt(c,&a); }
static inline int fff_invcbrt(ff_t *c,ff_t a) { ff_t x; return ff_cbrt_invcbrt(&x,c,&a); }

static inline ff_t fff_sum_2_muls(ff_t a1,ff_t a2, ff_t b2, ff_t b1) // returns a1*b1 + a2*b2
    { ff_t c; _ff_sum_2_mults (c,a1,a2,b2,b1); return c; }
static inline ff_t fff_sum_3_muls(ff_t a1,ff_t a2, ff_t a3, ff_t b3, ff_t b2, ff_t b1)
    { ff_t c; _ff_sum_3_mults (c,a1,a2,a3,b3,b2,b1); return c; }
static inline ff_t fff_sum_4_muls(ff_t a1,ff_t a2, ff_t a3, ff_t a4, ff_t b4, ff_t b3, ff_t b2, ff_t b1)
    { ff_t c; _ff_sum_4_mults (c,a1,a2,a3,a4,b4,b3,b2,b1); return c; }
static inline ff_t fff_sum_5_muls(ff_t a1,ff_t a2, ff_t a3, ff_t a4, ff_t a5, ff_t b5, ff_t b4, ff_t b3, ff_t b2, ff_t b1)
    { ff_t c; _ff_sum_5_mults (c,a1,a2,a3,a4,a5,b5,b4,b3,b2,b1); return c; }
static inline ff_t fff_sum_6_muls(ff_t a1,ff_t a2, ff_t a3, ff_t a4, ff_t a5, ff_t a6, ff_t b6, ff_t b5, ff_t b4, ff_t b3, ff_t b2, ff_t b1)
    { ff_t c; _ff_sum_6_mults (c,a1,a2,a3,a4,a5,a6,b6,b5,b4,b3,b2,b1); return c; }

static inline ff_t fff_dot2 (ff_t a1,ff_t b1, ff_t a2, ff_t b2) // returns a1*b1 + a2*b2
    { ff_t c; _ff_sum_2_mults (c,a1,a2,b2,b1); return c; }
static inline ff_t fff_dot_product(ff_t V[], ff_t W[], int r)
    { ff_t c; ff_dot_product(&c,V,W,r); return c; }

static inline int fff_poly_degree (ff_t f[], int d)
    { while ( d >= 0 && fff_is_zero(f[d]) ) d--; return d; }
static inline ff_t *fff_poly_copy (ff_t g[], ff_t f[], int d)
    { memcpy(g,f,(d+1)*sizeof(ff_t)); return g; }
static inline int fff_poly_from_ui (ff_t f[], unsigned long F[], int d) { for ( int i = 0 ; i <= d ; i++ ) f[i] = fff_from_ui (F[i]); return fff_poly_degree(f,d); }
static inline int fff_poly_from_si (ff_t f[], long F[], int d) { for ( int i = 0 ; i <= d ; i++ ) f[i] = fff_from_si (F[i]); return fff_poly_degree(f,d); }
static inline int fff_poly_from_int (ff_t f[], long F[], int d) { return fff_poly_from_si(f,F,d); }
static inline int fff_poly_from_mpz (ff_t f[], mpz_t F[], int d) { for ( int i = 0 ; i <= d ; i++ ) f[i] = fff_from_mpz (F[i]); return fff_poly_degree(f,d); }
static inline unsigned long *fff_poly_to_ui (unsigned long F[], ff_t f[], int d) { for ( int i = 0 ; i <= d ; i++ ) F[i] = fff_to_ui (f[i]); return F; }
static inline long *fff_poly_to_si (long F[], ff_t f[], int d) { for ( int i = 0 ; i <= d ; i++ ) F[i] = fff_to_si (f[i]); return F; }
static inline long *fff_poly_to_int (long F[], ff_t f[], int d) { return fff_poly_to_si (F,f,d); }
static inline mpz_t *fff_poly_to_mpz (mpz_t F[], ff_t f[], int d) { for ( int i = 0 ; i <= d ; i++ ) fff_to_mpz (F[i], f[i]);  return F; }
static inline void fff_poly_print (ff_t f[], int d)
    { ff_poly_print(f, d); }

static inline ff_t *fff_poly_mul (ff_t h[], ff_t f[], int df, ff_t g[], int dg) { ff_poly_mult (h, 0, f, df, g, dg); return h; }
static inline ff_t *fff_poly_sqr (ff_t g[], ff_t f[], int d) { ff_poly_square (g, f, d); return g; }
static inline ff_t *fff_poly_exp (ff_t g[], ff_t f[], int d, uint64_t e)
{
    ff_t *x, _x[d+1];
    if ( !e ) { g[0] = fff_one(); return g; }
    if ( d < 0 ) return g;
    if ( g == f ) { memcpy (_x,f,(d+1)*sizeof(x[0])); x = _x; } else { x = f; }
    memcpy (g,x,(d+1)*sizeof(g[0])); uint64_t dg = d;
    for ( int i = ui_len(e)-2 ; i >= 0 ; i-- ) {
        fff_poly_sqr (g, g, dg); dg *= 2;
        if ( (e&(1UL<<i)) ) { fff_poly_mul (g, g, dg, x, d); dg += d; }
    }
    return g;
}

static inline ff_t fff_linear_root (ff_t f[2])
    { ff_t r = fff_neg(f[0]); return f[1] == fff_one() ? r : fff_mul(r,fff_inv(f[1])); }

static inline ff_t fff_quadratic_disc (ff_t f[3])
    { return fff_sum_2_muls(fff_dbl(f[0]),f[1],f[1],fff_neg(fff_dbl(f[2]))); }

// computes the discriminant -27*f0^2*f3^2 + 18*f0*f1*f2*f3 - 4*f0*f2^3 - 4*f1^3*f3 + f1^2*f2^2 of a general cubic
// = 9*f0*f3(-3*f0*f3+2*f1*f2) + (-4*f0*f2)*f2*2 + f1^2*(f2^2-4*f1*f3) using 9M+15A and 7redc
static inline ff_t fff_cubic_disc (ff_t f[4])
{
    ff_t a1 = fff_mul3(fff_mul(f[0],f[3]))          ;                           // 3*f0*f3
    ff_t b1 = fff_add(fff_neg(a1),fff_dbl(fff_mul(f[1],f[2])));                 // -3*f0*f3+2*f1*f2
    a1 = fff_mul3(a1);                                                          // 9*f0*f3
    ff_t a2 = fff_neg(fff_dbl(fff_dbl(fff_mul(f[0],f[2]))));                    // -4*f0*f2
    ff_t b2 = fff_sqr(f[2]);                                                    // f2^2
    ff_t a3 = fff_sqr(f[1]);                                                    // f1^2
    ff_t b3 = fff_sub(b2,fff_dbl(fff_dbl(fff_mul(f[1],f[3]))));                 // f2^2-4*f1*f3
    return fff_sum_3_muls(a1,a2,a3,b3,b2,b1);
}

static inline ff_t *fff_poly_depress_cubic (ff_t f[3])
    { ff_t t; ff_poly_depress_cubic(&t, f); return f; }

static inline int fff_poly_mod (ff_t g[], ff_t f[], int df, ff_t m[], int dm) { int dg; ff_poly_mod (g, &dg, f, df, m, dm); return dg; }
static inline int fff_poly_exp_mod (ff_t g[], ff_t f[], int df, uint64_t e, ff_t m[], int dm) { int dg; ff_poly_pow_mod (g, &dg, f, df, e, m, dm); return dg; }

static inline ff_t fff_poly_eval (ff_t f[], int d, ff_t x)
    { if ( d < 0 ) return fff_zero(); ff_t y = f[d]; for ( int i = d-1 ; i >= 0 ; i-- ) y = fff_add(fff_mul(y,x),f[i]); return y; }
static inline ff_t *fff_poly_scale (ff_t g[], ff_t f[], int d, ff_t a)
    { ff_poly_scale (g,f,d,a); return g; }
static inline ff_t *fff_poly_translate_small (ff_t g[], ff_t f[], int d, ff_t c)
{
    if ( g != f ) memcpy(g,f,(d+1)*sizeof(ff_t));
    ff_t w = fff_mul(g[d],c);
    for ( int i = d-1 ; i >= 0 ; i-- ) { for ( int j = i ; j < d-1 ; j++ ) { g[j] = fff_addmul(g[j],g[j+1],c); } g[d-1] = fff_add(g[d-1],w); }
    return g;
}
static inline ff_t *fff_poly_translate (ff_t g[], ff_t f[], int d, ff_t a)
    { if ( d < 32 ) fff_poly_translate_small (g,f,d,a); else ff_poly_translate (g,f,d,a); return g; }
static inline ff_t *fff_poly_derivative (ff_t g[], ff_t f[], int d)
    { return ff_poly_derivative (g,0,f,d); }
static inline ff_t *fff_poly_monic (ff_t g[], ff_t f[], int d)
    { return ff_poly_monic (g,0,f,d); }

static inline int fff_poly_count_distinct_roots (ff_t f[], int d)
    { return ff_poly_count_distinct_roots (f,d); }
static inline int fff_poly_distinct_roots (ff_t r[], ff_t f[], int d) 
    { return ff_poly_distinct_roots (r,f,d); }
static inline int fff_poly_find_root (ff_t *r, ff_t f[], int d)
    { return ff_poly_find_root(r,f,d); }
static inline int fff_poly_squarefree (ff_t f[], int d)
    { return ff_poly_squarefree (f,d); }
static inline int fff_poly_gcd (ff_t h[], ff_t f[], int df, ff_t g[], int dg)
    { int dh; return ff_poly_gcd (h, &dh, f, df, g, dg); }
static inline ff_t *fff_poly_remove_root (ff_t g[], ff_t f[], int d, ff_t r)
    { ff_poly_remove_root (g,f,d,&r); return g; }

// computes y[i]=f(x[i]) for 0 <= i < n, where f in Fp[x] has degree d
static inline ff_t *fff_poly_multipoint_eval (ff_t y[], ff_t f[], int d, ff_t x[], int n)
    { ff_poly_multipoint_eval (y,f,d,x,n); return y; }
// computes unique f in Fp[x] of degree d < n for which f(x[i])=f(y[i]) for 0 <= i < n (sets d to exact degree)
static inline ff_t *fff_poly_interpolate (ff_t f[], int *d, ff_t x[], ff_t y[], int n)
    { ff_poly_interpolate (f,d,x,y,n); return f; }
// Given a polynomial f computes f(0)f(1)...f(n-1) using Bostan-Gaudry-Schost O(M(sqrt(n))) algorithm
static inline ff_t fff_poly_multipoint_product (ff_t f[], int d, long n) // returns f(0)*...*f(n-1)
    { ff_t y; ff_poly_multipoint_product (&y,f,d,n); return y; }

// Computes B=A(x) for an r x r matrix A whose entries are polynomials of degree <= d at x.
static inline ff_t *fff_poly_matrix_eval (ff_t B[], ff_t A[], int r, int d, ff_t x)
    { ff_poly_matrix_eval (B, A, r, d, &x); return B; }
// Computes B(x)=A(a*x) for an r x r matrix A whose entries are polynomials of degree <= d.
static inline ff_t *fff_poly_matrix_scale (ff_t B[], ff_t A[], int r, int d, ff_t a)
    { ff_poly_matrix_scale (B, A, r, d, a); return A; }
// Computes B(x)=A(x+a) for an r x r matrix A whose entries are polynomials of degree <= d.
static inline ff_t *fff_poly_matrix_translate (ff_t B[], ff_t A[], int r, int d, ff_t a)
    { ff_poly_matrix_translate (B, A, r, d, a); return A; }
// Computes C(x)=A(x)*B(x) where A and B are r x r matrices whose entries are polynomials of degree <= d (C will have entries of degree <= 2d).
static inline ff_t *fff_poly_matrix_mul (ff_t C[], ff_t A[], ff_t B[], int r, int d)
    { ff_poly_matrix_mult (C, A, B, r, d); return C; }
// Given an r x r matrix A(x) whose entries are polynomials of degree <= d computes A(0)A(1)...A(n-1) using BGS
static inline ff_t *fff_poly_matrix_multipoint_product (ff_t A[], ff_t B[], int r, int d, long n)
    { ff_poly_matrix_multipoint_product (A,B,r,d,n); return A; }
// Given an r x r matrix A(x) whose entries are polynomials of degree <= d computes A(0)A(1)...A(n-1) using BGS but scales degree up by 2^k first
static inline ff_t *fff_poly_matrix_multipoint_product_scale (ff_t A[], ff_t B[], int r, int d, long n, int k)
    { ff_poly_matrix_multipoint_product_scale (A,B,r,d,n,k); return A; }

// Replaces A with its transpose
static inline ff_t *fff_poly_matrix_transpose_inplace (ff_t A[], int r, int d)
{
    for ( int i = 0, n = d+1 ; i < r-1 ; i++ ) for ( int j = i+1 ; j < r ; j++ )
        { ff_t x[n]; memcpy(x,A+n*(i*r+j),sizeof(x)); memcpy (A+n*(i*r+j),A+n*(j*r+i),sizeof(x)); memcpy (A+n*(j*r+i),x,sizeof(x)); } 
    return A;
}
// Sets B to the transpose of A, aliasing okay
static inline ff_t *fff_poly_matrix_transpose (ff_t B[], ff_t A[], int r, int d)
    { if ( B != A ) { memmove(B,A,(d+1)*r*r*sizeof(A[0])); } return fff_poly_matrix_transpose_inplace (B, r, d); }

// Computes C=A*B where A and B are r x r matrices
static inline ff_t *fff_matrix_mul (ff_t C[], ff_t A[], ff_t B[], int r)
    { ff_matrix_mult (C, A, B, r); return C; }
// Computes B=A*B where A and B are r x r matrices
static inline ff_t *fff_matrix_lmul_inplace (ff_t B[], ff_t A[], int r)
    { ff_t C[r*r]; ff_matrix_mult (C, A, B, r); memcpy(B,C,r*r*sizeof(ff_t)); return B; }
// Computes B=B*A where A and B are r x r matrices
static inline ff_t *fff_matrix_rmul_inplace (ff_t B[], ff_t A[], int r)
    { ff_t C[r*r]; ff_matrix_mult (C, B, A, r); memcpy(B,C,r*r*sizeof(ff_t)); return B; }
// Computes W=A*V where A is an r x s matrix and V is a vector of length s (no overlap)
static inline ff_t *fff_matrix_vector_mul (ff_t W[], ff_t A[], ff_t V[], int r, int s)
    { for ( int i = 0 ; i < r ; i++ ) { W[i] = fff_dot_product(A+i*s,V,s); } return W; }
// Computes V=A*V where A is an r x r matrix and V is a vector of length r
static inline ff_t *fff_matrix_vector_mul_inplace (ff_t V[], ff_t A[], int r)
    { ff_t W[r]; for ( int i = 0 ; i < r ; i++ ) { W[i] = fff_dot_product(A+i*r,V,r); } memcpy(V,W,r*sizeof(ff_t)); return V; }
// Computes W=V*A where A is an r x s matrix and V is a vector of length r (no overlap)
static inline ff_t *fff_vector_matrix_mul (ff_t W[], ff_t V[], ff_t A[], int r, int s)
    { for ( int j = 0 ; j < s ; j++ ) { ff_t v[r]; for ( int i = 0 ; i < r ; i++ ) { v[i] = A[i*s+j]; } W[j] = fff_dot_product(v,V,r); } return W; }
// Computes V=V*A where A is an r x r matrix and V is a vector of length r
static inline ff_t *fff_vector_matrix_mul_inplace (ff_t V[], ff_t A[], int r)
    { ff_t W[r]; fff_vector_matrix_mul (W, V, A, r, r); memcpy (V,W,sizeof(W)); return V; }
// Replaces A with its transpose
static inline ff_t *fff_matrix_transpose_inplace (ff_t A[], int r)
    { for ( int i = 0 ; i < r-1 ; i++ ) for ( int j = i+1 ; j < r ; j++ ) { ff_t x = A[i*r+j]; A[i*r+j] = A[j*r+i]; A[j*r+i] = x; }  return A; }
// Sets B to the transpose of A, aliasing okay
static inline ff_t *fff_matrix_transpose (ff_t B[], ff_t A[], int r)
    { if ( B != A ) { memmove(B,A,r*r*sizeof(A[0])); } return fff_matrix_transpose_inplace (B, r); }

static inline int fff_norm_equation (long *t, long *v, long D) // computes solution to 4p=t^2-v^2D (assuming one exists)
    { return ff_norm_equation (t,v,D); }

#ifdef __cplusplus
}
#endif

#endif
