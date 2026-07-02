/*

zp_poly_exp.c: exponentiation functions for zp_poly

Copyright (C) 2010-2012, David Harvey and Andrew V. Sutherland

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <assert.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"
#include "zp_roots.h"

static int zp_poly_algo = ZP_POLY_MOD_ALGO_AUTO;


// computes h =f^n mod g (both f and h are zero padded to degree deg(g) -1)
// h and f may coincide
// TODO: change to use more efficient exponentiation, e.g. sliding window
void zp_poly_mod_pow (zp_t h[], zp_t f[], mpz_t n, zp_poly_mod_t mod)
{
    zp_t *g;
    register int t;
    
    if ( h == f ) g = zp_poly_mod_alloc (mod); else g = h;
    zp_poly_mod_set_one (g, mod);
    for ( t = mpz_sizeinbase (n, 2) - 1 ; t >= 0 ; t-- ) {
        zp_poly_mod_mul (g,g,g,mod);
        if ( mpz_tstbit(n,t) ) zp_poly_mod_mul (g, g, f, mod);
    }
    if ( h == f ) { zp_poly_mod_copy (h, g, mod); zp_poly_mod_free (g, mod); }
}

// computes h = x^n mod g (zero padded to degree deg(g) -1)
void zp_poly_mod_pow_xn (zp_t h[], zp_t n, zp_poly_mod_t mod)
{
#if ZP_POLY_EXP_TIME_TRACE
    clock_t start, end;
    int e;
#endif
    register int t;
    
#if ZP_POLY_EXP_TIME_TRACE
    printf ("zp_poly_mod_pow_xn deg=%d ebits=%ld\n", mod->d, mpz_sizeinbase (n,2)-ui_bits(mod->d));     // don't count the bits up to d, these are easy
    e = -1;
#endif
    zp_poly_mod_set_one (h, mod);
    for ( t = mpz_sizeinbase (n, 2) - 1 ; t >= 0 ; t-- )  {
        zp_poly_mod_mulx (h,h,h,mpz_tstbit(n,t),mod);
#if ZP_POLY_EXP_TIME_TRACE
        if ( e >= 0 ) { e++; end = clock(); printf ("    %6d ebits in %10.3f secs  %10.6f secs/bit\r", e, (double)(end-start)/CLOCKS_PER_SEC,  (double)(end-start)/(CLOCKS_PER_SEC*e)); fflush(stdout); }
        if ( e<0 && mpz_sizeinbase (n, 2)-t > ui_bits (mod->d) ) { e = 0; start = clock(); }
#endif
    }
#if ZP_POLY_EXP_TIME_TRACE
    puts ("");
#endif
}

// h = (x+a) * f mod g, h may alias f
static inline void _zp_poly_mod_mul_by_xa (mpz_t h[], mpz_t f[], mpz_t a, zp_poly_mod_t mod, mpz_t w1, mpz_t w2)
{
    register int i;

    mpz_set (w1, f[mod->d-1]);
    for ( i = mod->d-1 ; i >= 0 ; i-- ) {
        // h[i] = a*h[i] - h[d-1]*g[i] + h[i-1]  (for i > 0)
        mpz_mul (w2, a, f[i]);  mpz_submul (w2, w1, mod->g[i]);
        if ( i ) mpz_add (w2, w2, f[i-1]);
        mpz_mod (h[i], w2, mod->p);
    }
}

// computes h = (x+a)^n mod g (zero padded to degree deg(g) -1)
void zp_poly_mod_pow_xan (zp_t h[], zp_t a, mpz_t n, zp_poly_mod_t mod)
{
#if ZP_POLY_EXP_TIME_TRACE
    clock_t start, end;
#endif
    register int t;
    mpz_t w1, w2;
    
#if ZP_POLY_EXP_TIME_TRACE
    printf ("zp_poly_mod_pow_xan deg=%d ebits=%ld\n", mod->d, mpz_sizeinbase (n,2));
    start = clock();
#endif
    mpz_init_modp_w (w1, mod->p);  mpz_init_modp_w (w2, mod->p);
    zp_poly_mod_set_one (h, mod);
    for ( t = mpz_sizeinbase (n, 2) - 1 ; t >= 0 ; t-- ) {
        zp_poly_mod_mul (h,h,h,mod);
        if ( mpz_tstbit(n,t) ) _zp_poly_mod_mul_by_xa (h, h, a, mod, w1, w2);
#if ZP_POLY_EXP_TIME_TRACE
        end = clock();
        printf ("    %6ld ebits in %10.3f secs  %10.6f secs/bit\r", mpz_sizeinbase (n, 2)-t, (double)(end-start)/CLOCKS_PER_SEC,  (double)(end-start)/(CLOCKS_PER_SEC*(mpz_sizeinbase(n,2)-t))); fflush(stdout);
#endif
    }
#if ZP_POLY_EXP_TIME_TRACE
    puts ("");
#endif
    mpz_clear (w1);  mpz_clear (w2);
}

// h = (x^3+ax+b) * f mod g, where g has degree > 3.  h may alias f
static inline void _zp_poly_mod_mul_by_x3axb (mpz_t h[], mpz_t f[], mpz_t a, mpz_t b, zp_poly_mod_t mod, mpz_t w, mpz_t q[3])
{
    mpz_t *g;
    register int i, d;

    g = mod->g;  d = mod->d;
    assert (d > 3);
    // compute q = [(x^3+a*x+b)*f] / g with degree 2
    mpz_set (q[2], f[d-1]);
    mpz_mul_mod (w, q[2], g[d-1], mod->p);  mpz_sub_mod (q[1], f[d-2], w, mod->p);
    mpz_mul (w, a, f[d-1]); mpz_add(w, w, f[d-3]);  mpz_submul (w, q[2], g[d-2]);  mpz_submul (w, q[1], g[d-1]); mpz_mod (q[0], w, mod->p);
    // compute (x^3+a*x+b) * f - q * g of degree d-1, working from the top down
    for ( i = d-1 ; i >= 3 ; i-- ) {
        // deg i coefficient is w = f[i-3] + a*f[i-1] + b*f[i] - q[2]*g[i-2] - q[1]*g[i-1] - q[0]*g[i]
        mpz_mul (w, a, f[i-1]);  mpz_add (w, w, f[i-3]);  mpz_addmul (w, b, f[i]); mpz_submul (w, q[2], g[i-2]); mpz_submul (w, q[1], g[i-1]); mpz_submul (w, q[0], g[i]);
        mpz_mod (h[i], w, mod->p);
    }
    // deg 2 coefficient
    mpz_mul (w, a, f[1]);  mpz_addmul (w, b, f[2]); mpz_submul (w, q[2], g[0]); mpz_submul (w, q[1], g[1]); mpz_submul (w, q[0], g[2]);
    mpz_mod (h[2], w, mod->p);
    // deg 1 coefficient
    mpz_mul (w, a, f[0]);  mpz_addmul (w, b, f[1]); mpz_submul (w, q[1], g[0]); mpz_submul (w, q[0], g[1]);
    mpz_mod (h[1], w, mod->p);
    // deg 0 coefficient
    mpz_mul (w, b, f[0]); mpz_submul (w, q[0], g[0]);
    mpz_mod (h[0], w, mod->p);
}

void zp_poly_mod_mul_by_x3axb (zp_t h[], zp_t f[], zp_t a, zp_t b, zp_poly_mod_t mod)
{
    mpz_t w, q[3];
    
    mpz_init_modp_w (w, mod->p);  mpz_init_modp (q[0], mod->p); mpz_init_modp (q[1], mod->p); mpz_init_modp (q[2], mod->p);
    if ( mod->d > 3 ) {
        _zp_poly_mod_mul_by_x3axb (h, f, a, b, mod,  w, q);
    } else {
        mpz_set_ui (w, 3);
        zp_poly_mod_pow_xn (q, w, mod);
        if ( mod->d > 1 ) { mpz_add_mod (q[1], q[1], a, mod->p); } else { mpz_mul_mod (w, a, mod->g[0], mod->p);  mpz_sub_mod (q[0], q[0], w, mod->p); }
        mpz_add_mod (q[0], q[0], b, mod->p);
        zp_poly_mod_mul (h, f, q, mod);
    }
    mpz_clear (w); mpz_clear (q[0]); mpz_clear (q[1]); mpz_clear (q[2]);
}

// computes h = (x^3+ax+b)^n mod g (zero padded to degree deg(g) -1)
void zp_poly_mod_pow_x3axbn (zp_t h[], zp_t a, zp_t b, mpz_t n, zp_poly_mod_t mod)
{
    register int t;
    mpz_t w, q[3];
    
    mpz_init_modp_w (w, mod->p);  mpz_init_modp (q[0], mod->p); mpz_init_modp (q[1], mod->p); mpz_init_modp (q[2], mod->p);
    zp_poly_mod_set_one (h, mod);
    if ( mod->d > 3 ) {
        for ( t = mpz_sizeinbase (n, 2) - 1 ; t >= 0 ; t-- ) {
            zp_poly_mod_mul (h,h,h,mod);
            if ( mpz_tstbit(n,t) ) _zp_poly_mod_mul_by_x3axb (h, h, a, b, mod,  w, q);
        }
    } else {
        // compute q = x^3 + ax + b mod g
        mpz_set_ui (w, 3);
        zp_poly_mod_pow_xn (q, w, mod);
        if ( mod->d > 1 ) { mpz_add_mod (q[1], q[1], a, mod->p); } else { mpz_mul_mod (w, a, mod->g[0], mod->p);  mpz_sub_mod (q[0], q[0], w, mod->p); }
        mpz_add_mod (q[0], q[0], b, mod->p);
        // we could easily make this more efficient by optimizing for low degree q
        for ( t = mpz_sizeinbase (n, 2) - 1 ; t >= 0 ; t-- ) {
            zp_poly_mod_mul (h,h,h,mod);
            if ( mpz_tstbit(n,t) ) zp_poly_mod_mul (h, h, q, mod);
        }       
    }
    mpz_clear (w);  mpz_clear (q[0]);  mpz_clear (q[1]);  mpz_clear (q[2]);
}

// computes h = gcd(x^p-x,f) = monic product of the linear factors of a monic f in Z/pZ[x].
void zp_poly_linear_factor_product (zp_t h[], int *dh, zp_t f[], int d, zp_mod_t p)
{
#if ZP_POLY_EXP_TIME_TRACE
    clock_t start, end;
#endif
    zp_t *g;
    zp_poly_mod_t mod;
    int dg;
    
    // handle d < 2
    assert ( zp_poly_is_monic (f,d) && d >= 0);         // factoring 0 is undefined
    if ( !d ) { mpz_set_one(h[*dh=0]); return; }
    if ( d== 1 ) { zp_poly_copy (h,f,*dh=d); return; }
    
    // compute g=x^p mod f, then g -= x, then h = gcd(f,g)
    g = zp_poly_alloc (d-1, p);
#if ZP_POLY_EXP_TIME_TRACE
    printf ("begin zp_poly_linear_factor_product: deg=%d, pbits=%ld, algo=%d\n", d, zp_bits(p), zp_poly_algo);
    start = clock();
#endif
    zp_poly_mod_init (mod, f, d, p, zp_poly_algo);
#if ZP_POLY_EXP_TIME_TRACE
    end = clock();
    printf("zp_poly_mod_init deg %d took %.3f secs\n", d, (double) (end-start)/CLOCKS_PER_SEC);
#endif
    zp_poly_mod_pow_xn (g, p->p, mod);
    zp_poly_mod_clear (mod);
    zp_dec (g[1], p);
    dg = zp_poly_degree (g, d-1);
#if ZP_POLY_EXP_TIME_TRACE
    start = clock();
#endif
    zp_poly_gcd (h, dh, f, d, g, dg, p);
#if ZP_POLY_EXP_TIME_TRACE
    end = clock();
    printf ("zp_poly_gcd (deg %d, deg %d) = deg %d took %.3f secs\n", d, dg, *dh, (double) (end-start)/CLOCKS_PER_SEC);
    printf ("end zp_poly_linear_factor_product: deg=%d -> %d, pbits=%ld, algo=%d\n", d, *dh, zp_bits(p), zp_poly_algo);
#endif
    zp_poly_free (g, d-1);
}

// finds a root of a nonconstant monic poly that splits completely in Z/pZ[x] into *distinct* linear factors
// uses a randomized algorithm if d > 2, trashes h in the process
void zp_poly_find_split_root (zp_t r, zp_t h[], int d, zp_mod_t p)
{
#if ZP_POLY_EXP_TIME_TRACE
    clock_t start, end;
#endif
    zp_poly_mod_t mod;
    zp_t a;
    mpz_t n;
    int i, dh;
    
    assert (d > 0 && zp_poly_is_monic (h, d) );
    switch (d) {
    case 1: zp_neg (r, h[0], p); return;
    case 2: if ( ! zp_poly_quadratic_roots (h, h, p) ) { fprintf (stderr, "irreducible quadratic in zp_poly_split_roots!\n"); abort(); }    zp_set (r, h[0]); return;
    case 3: if ( ! zp_poly_cubic_roots (h, h, p) ) { fprintf (stderr, "irreducible cubic in zp_poly_split_roots!\n"); abort(); }  zp_set (r, h[0]); return;
    }
    
#if ZP_POLY_EXP_TIME_TRACE
    printf ("begin zp_poly_find_split_root deg=%d, pbits=%ld, algo=%d\n", d, zp_bits(p), zp_poly_algo);
    start = clock();
#endif
    zp_poly_mod_init (mod, h, d, p, zp_poly_algo);
#if ZP_POLY_EXP_TIME_TRACE
    end = clock();
    printf("zp_poly_mod_init deg %d took %.3f secs\n", d, (double) (end-start)/CLOCKS_PER_SEC);
#endif
    zp_init (a, p); mpz_init_modp (n, p->p);
    mpz_sub_ui(n, p->p, 1); mpz_div_2exp (n, n, 1);     // n = (p-1)/2
    i = 0;
    do {
        zp_random (a, p);
        zp_poly_mod_pow_xan (h, a, n, mod);  zp_dec (h[0], p); dh = zp_poly_degree (h, d-1);
        if ( dh > 0 ) {
#if ZP_POLY_EXP_TIME_TRACE
            start = clock();
            printf ("gcd(deg %d,deg %d) took ", mod->d, dh);
#endif
            zp_poly_gcd (h, &dh, mod->g, mod->d, h, dh, p);
#if ZP_POLY_EXP_TIME_TRACE
            end = clock();
            printf ("%.3f secs\n", (double) (end-start)/CLOCKS_PER_SEC);
#endif
        } else {
            i += d;
            if ( i > 100 ) { printf ("Unable to split degree %d poly in zp_poly_find_split_root, are you sure the input poly really splits?\n", d); zp_poly_print (h, d); abort(); }
        }
    } while ( dh <= 0 );

    // be sure to use the smaller factor when we recurse
    if ( dh > d-dh ) { zp_poly_div (h, mod->g, mod->d , h, dh, p); dh = d-dh; }
    zp_poly_mod_clear (mod);
    mpz_clear (a);  mpz_clear (n);
    zp_poly_find_split_root (r, h, dh, p);
#if ZP_POLY_EXP_TIME_TRACE
    printf ("end zp_poly_find_split_root: deg=%d, pbits=%ld, algo=%d\n", d, zp_bits(p), zp_poly_algo);
#endif
}

// finds min(d,k) roots of a nonconstant monic poly that splits completely in Z/pZ[x] into *distinct* linear factors
// uses a randomized algorithm if d > 2, trashes h in the process, h cannot overlap r
void zp_poly_find_split_roots (zp_t r[], int k, zp_t h[], int d, zp_mod_t p)
{
#if ZP_POLY_EXP_TIME_TRACE
    clock_t start, end;
#endif
    zp_poly_mod_t mod;
    zp_t *g;
    zp_t a;
    mpz_t n;
    int i, dg, dh;
    
    if ( k == 1 ) { zp_poly_find_split_root (r[0], h, d, p); return; }
    assert (d > 0 && k > 0 && zp_poly_is_monic (h, d) );
    switch (d) {
    case 1: zp_neg (r[0], h[0], p); return;
    case 2: if ( zp_poly_quadratic_roots (h, h, p) != 2 ) { fprintf (stderr, "irreducible quadratic in zp_poly_split_roots!\n"); abort(); } for ( i = 0 ; i < k && i < 2 ; i++ ) mpz_set (r[i], h[i]); return;
    case 3: if ( zp_poly_cubic_roots (h, h, p) != 3 ) { fprintf (stderr, "irreducible cubic in zp_poly_split_roots!\n"); abort(); } for ( i = 0 ; i < k && i < 3 ; i++ ) mpz_set (r[i], h[i]); return;
    }
    
#if ZP_POLY_EXP_TIME_TRACE
    printf ("begin zp_poly_find_split_roots: roots=%d, deg=%d, pbits=%ld, algo=%d\n", k, d, zp_bits(p), zp_poly_algo);
    start = clock();
#endif
    zp_poly_mod_init (mod, h, d, p, zp_poly_algo);
#if ZP_POLY_EXP_TIME_TRACE
    end = clock();
    printf("zp_poly_mod_init deg %d took %.3f secs\n", d, (double) (end-start)/CLOCKS_PER_SEC);
#endif
    zp_init (a, p); mpz_init_modp (n, p->p);
    mpz_sub_ui(n, p->p, 1); mpz_div_2exp (n, n, 1);     // n = (p-1)/2
    i = 0;
    do {
        zp_random (a, p);
        zp_poly_mod_pow_xan (h, a, n, mod);  zp_dec (h[0], p); dh = zp_poly_degree (h, d-1);
        if ( dh > 0 ) {
#if ZP_POLY_EXP_TIME_TRACE
            start = clock();
            printf ("gcd(deg %d,deg %d) took ", mod->d, dh);
#endif
            zp_poly_gcd (h, &dh, mod->g, mod->d, h, dh, p);
#if ZP_POLY_EXP_TIME_TRACE
            end = clock();
            printf ("%.3f secs\n", (double) (end-start)/CLOCKS_PER_SEC);
#endif
        } else {
            i += d;
            if ( i > 100 ) { printf ("Unable to split degree %d poly in zp_poly_find_split_root, are you sure the input poly really splits?\n", d); zp_poly_print (h, d); abort(); }
        }
    } while ( dh <= 0 );
    
    // if we can get all the roots we want from the smaller factor, do so
    if ( k < dh && k < d-dh ) {
        if ( dh > d-dh ) { zp_poly_div (h, mod->g, mod->d , h, dh, p); dh = d-dh; }
        zp_poly_mod_clear (mod);
        mpz_clear (a);  mpz_clear (n);
        zp_poly_find_split_roots (r, k, h, dh, p);
        return;
    }
    // if we can get all the roots we want from the bigger factor, do so
    if ( k < dh || k < d-dh ) {
        if ( dh < d-dh ) { zp_poly_div (h, mod->g, mod->d , h, dh, p); dh = d-dh; }
        zp_poly_mod_clear (mod);
        mpz_clear (a);  mpz_clear (n);
        zp_poly_find_split_roots (r, k, h, dh, p);
        return;
    }
    g = zp_poly_alloc (d-dh, p);
    zp_poly_div (g, mod->g, mod->d , h, dh, p); dg = d-dh;
    zp_poly_mod_clear (mod);
    mpz_clear (a);  mpz_clear (n);
    zp_poly_find_split_roots (r, dg, g, dg, p);
    zp_poly_free (g, d-dh);
    zp_poly_find_split_roots (r+dg, k-dg, h, dh, p);
#if ZP_POLY_EXP_TIME_TRACE
    printf ("end zp_poly_find_split_roots: roots=%d, deg=%d, pbits=%ld, algo=%d\n", k, d, zp_bits(p), zp_poly_algo);
#endif
}

// attempts to find a root r of a monic poly f of deg d > 0
// returns the total number of *distinct* roots
int zp_poly_find_root (zp_t r, zp_t f[], int d, zp_mod_t p)
{
    zp_t *h;
    int dh;

    assert ( zp_poly_is_monic(f, d) && d > 0 );
    h = zp_poly_alloc (d, p);
    switch (d) {
    case 1: zp_neg (r, f[0], p); dh=1; break;
    case 2:
        dh = zp_poly_quadratic_roots (h, f, p); 
        if ( dh ) { zp_set (r, h[0]); if ( zp_equal (h[0],h[1]) ) dh = 1; } // adjust count if roots are not distinct
        break;
    case 3:
        dh = zp_poly_cubic_roots (h, f, p);
        if ( dh ) { zp_set (r, h[0]); if ( dh == 3 ) { if ( zp_equal (h[0],h[1]) ) dh--;  if ( zp_equal(h[0],h[2]) ) dh--;  if ( dh ==3 && zp_equal (h[1],h[2]) ) dh--; } } // adjust count if roots are not distinct
        break;
    default:
        zp_poly_linear_factor_product (h, &dh, f, d, p);
        if ( dh > 0 ) zp_poly_find_split_root (r,  h, dh, p);
    }
    zp_poly_free (h, d);
    return dh;
}

// attempts to find up to k roots of a monic poly f of deg d > 0
// returns the total number of roots, with multiplicity (which may be greater than k, but only k roots will be placed in r[])
int zp_poly_find_roots (zp_t r[], int k, zp_t f[], int d, zp_mod_t p)
{
    zp_t *g, *h, *w;
    int i, n, dg, dw;

    assert ( zp_poly_is_monic(f, d) && d > 0 );
    switch (d) {
    case 1: if ( k >= 1 ) zp_neg (r[0], f[0], p); return 1;
    case 2:
        if ( k >= 2 ) return zp_poly_quadratic_roots (r, f, p);
        h = zp_poly_alloc (d, p);
        n = zp_poly_quadratic_roots (h, f, p); 
        if ( n && k >= 1) zp_set (r[0], h[0]);
        break;
    case 3:
        if ( k >= 3 ) return zp_poly_cubic_roots (r, f, p);
        h = zp_poly_alloc (d, p);
        n = zp_poly_cubic_roots (h, f, p);
        for ( i = 0 ; i < n && i < k ; i++ ) zp_set (r[i], h[i]);
        break;
    default:
        h = zp_poly_alloc (d, p);
        zp_poly_linear_factor_product (h, &n, f, d, p);
        if ( !n ) break;
        if ( n==d ) { zp_poly_find_split_roots (r,  k, h, n, p); break; }
        // save a copy of h since it will be destroyed by zp_poly_find_split_roots
        w = zp_poly_alloc (d, p); zp_poly_copy (w, h, n); dw = n;
        zp_poly_find_split_roots (r,  k, h, n, p);
        // we have more than 0 and less than d roots, we need to check for multiple roots
        g = zp_poly_alloc (d-1, p);
        zp_poly_derivative (g, f, d, p);  dg = d-1;
        for (;;) {
            zp_poly_gcd (w, &dw, g, dg, w, dw, p);
            if ( !dw ) break;
            if ( n < k ) zp_poly_find_split_roots (r+n, k-n, w, dw, p);
            n += dw;
            zp_poly_derivative (g, g, dg, p);  dg = dg-1;
        }
        zp_poly_free (g, d-1);  zp_poly_free (w, d);
    }
    zp_poly_free (h, d);
    return n;
}

// completely factors a square-free monic polynomial h that is the product of d/k irreducible polynomials of degree k, using Cantor-Zassenhaus (Algs 14.8 and 14.10 in GG)
// r must have space for d=(d/k)*k coefficients:  a concatenated list of d/k implicitly monic polys of degree k will be stored in w, each using just d/k coeffs (monic leading coeff is omitted).
// the polynomial h is destroyed in the process and cannot overlap r
void zp_poly_equal_degree_factor (zp_t r[], zp_t h[], int d, int k, zp_mod_t p)
{
    zp_poly_mod_t mod;
    zp_t *g;
    mpz_t n;
    int dg, dh;

    assert ( k > 0 && !(d%k) && zp_poly_is_monic (h, d) );
    if ( d == k ) { zp_poly_copy (r, h, d-1); return; }

    g = zp_poly_alloc (d, p);

    zp_poly_mod_init (mod, h, d, p, zp_poly_algo);  
    mpz_init2 (n, k*zp_bits(p));  mpz_pow_ui (n, p->p, k);  mpz_sub_ui(n, n, 1); mpz_div_2exp (n, n, 1);            // n = (p^k-1)/2
    do {
        zp_poly_randomize (g, d-1, p);
        zp_poly_mod_pow (h, g, n, mod);  zp_dec (h[0], p);  dh = zp_poly_degree (h, mod->d-1);              // h = g^n-1
        if ( dh > 0 ) zp_poly_gcd (h, &dh, mod->g, mod->d, h, dh, p);
    } while ( dh <= 0 );
    mpz_clear (n);
    
    assert ( !(dh%k) );
    
    zp_poly_div (g, mod->g, d, h, dh, p);  dg = d-dh;
    zp_poly_mod_clear (mod);
    zp_poly_equal_degree_factor (r, g, dg, k, p);
    zp_poly_free (g, d);
    zp_poly_equal_degree_factor (r+dg, h, dh, k, p);
}

// finds a single factor of a square-free monic polynomial h that is the product of d/k irreducible polynomials of degree k, using Cantor-Zassenhaus (Algs 14.8 and 14.10 in GG)
// the polynomial h is destroyed in the process and cannot overlap r
void zp_poly_find_equal_degree_factor (zp_t r[], zp_t h[], int d, int k, zp_mod_t p)
{
    zp_poly_mod_t mod;
    zp_t *g;
    mpz_t n;
    int dh;

    assert ( k > 0 && !(d%k) && zp_poly_is_monic (h, d) );
    if ( d == k ) { zp_poly_copy (r, h, d-1); return; }

    g = zp_poly_alloc (d, p);

    zp_poly_mod_init (mod, h, d, p, zp_poly_algo);  
    mpz_init2 (n, k*zp_bits(p));  mpz_pow_ui (n, p->p, k);  mpz_sub_ui(n, n, 1); mpz_div_2exp (n, n, 1);            // n = (p^k-1)/2
    do {
        zp_poly_randomize (g, d-1, p);
        zp_poly_mod_pow (h, g, n, mod);  zp_dec (h[0], p);  dh = zp_poly_degree (h, mod->d-1);              // h = g^n-1
        if ( dh > 0 ) zp_poly_gcd (h, &dh, mod->g, mod->d, h, dh, p);
    } while ( dh <= 0 );
    mpz_clear (n);
    
    assert ( !(dh%k) );
    
    // make sure we use the lower degree factor
    if ( dh > d-dh ) {
        zp_poly_div (g, mod->g, d, h, dh, p);
        zp_poly_copy (h, g, d-dh);  dh = d-dh;
    }
    zp_poly_mod_clear (mod);
    zp_poly_free (g, d);
    zp_poly_find_equal_degree_factor (r, h, dh, k, p);
}


// Determines the irreducible factors of a monic polynomial f, using the standard large characteristic algorithm using Cantor-Zassenhaus for equal-degree splitting  -- see Algorithm 14.13 in GG.
// TODO: optimize this algorithm, it is currently quite slow

// r should have space for d coefficients to hold a concatenated list of f's irreducible factors, stored as implicitly monic polys (so leading monic coeff is omitted), r and f may coincide
// n[i] is the degree of the ith factor,  and e[i] is its multiplicity.  return value is the number of distinct factors
// factors are returned in order of increasing degree, with repeated factors listed in succession. 
int zp_poly_factor (zp_t r[], int n[], int e[], zp_t f[], int d, zp_mod_t p)
{
    zp_poly_mod_t mod;
    zp_t *g, *h, *v, *w, *x;
    int i, j, k, m, dg, dv;
        
    assert ( zp_poly_is_monic (f, d) && d > 0 );
    if ( d == 1 ) { zp_poly_copy (r, f, d);  n[0] = 1; return 1; }

    g = zp_poly_alloc (d, p);  h = zp_poly_alloc (d, p);  v = zp_poly_alloc (d, p);  w = zp_poly_alloc (d/2, p);  x = zp_poly_alloc (d/2, p);
    zp_poly_mod_init (mod, f, d, p, zp_poly_algo);
    
    // initialize v to f, we will remove factors of f from v as we go
    zp_poly_copy (v, f, d);  dv = d;
    
    // j indexes the factors, k is the degree of the factors we are currently working on
    m = 0;
    for ( j = 0, k = 1 ; k <= dv/2 ; k++ ) {
        if ( k == 1 ) zp_poly_mod_pow_xn (g, p->p, mod); else zp_poly_mod_pow (g, h, p->p, mod);
        zp_poly_mod_copy (h, g, mod);
        zp_dec (g[1], p);           // g = x^(p^k)-x mod f
        dg = zp_poly_degree (g, d-1);
        zp_poly_gcd (g, &dg, g, dg, v, dv, p);
        if ( dg > 0 ) {
            assert (!(dg%k));
            zp_poly_equal_degree_factor (r+m, g, dg, k, p);
            for ( i = 0 ; i < dg/k ; i++, j++ ) {
                n[j] = k;  e[j] = 0;
                zp_poly_copy (w, r+m, k-1);  zp_set_one (w[k]);  m += k;
                // remove all factors of h from v, computing multiplicity of h as we go
                for ( e[j] = 0 ; dv ; e[j]++ ) {
                    zp_poly_div (g, v, dv, w, k, p);
                    if ( zp_poly_rem (x, v, dv, w, k, g, p) >= 0 ) break;
                    zp_poly_copy (v, g, dv-k);  dv -= k;
                }
                assert ( e[j] );
            }
        }       
    }
    
    // If anything is left of v, it must be irreducible
    if ( dv ) {
        assert ( dv >= k );
        zp_poly_copy (r+m, v, dv-1); n[j] = dv;  e[j] = 1; j++;
    }
    // verify degrees match up -- this can be removed eventually
    for ( i = k = 0 ; i < j ; i++ ) k += n[i]*e[i];
    assert ( k == d );
    zp_poly_mod_clear (mod);
    zp_poly_free (g, d);  zp_poly_free (h, d);  zp_poly_free (v, d);  zp_poly_free (w, d/2);   zp_poly_free (x, d/2); 
    return j;
}

// Finds a factor a monic poly f of degree d that has minimal positive degree.  Returns the degree
int zp_poly_min_degree_factor (zp_t r[], zp_t f[], int d, zp_mod_t p)
{
    zp_poly_mod_t mod;
    zp_t *g, *h;
    int k, dg;
        
    assert ( zp_poly_is_monic (f, d) && d > 0 );
    if ( d == 1 ) { zp_poly_copy (r, f, d);  return 1; }

    g = zp_poly_alloc (d, p);  h = zp_poly_alloc (d, p);
    zp_poly_mod_init (mod, f, d, p, zp_poly_algo);
    
    for ( k = 1 ; k <= d/2 ; k++ ) {
        if ( k == 1 ) zp_poly_mod_pow_xn (g, p->p, mod); else zp_poly_mod_pow (g, h, p->p, mod);
        zp_poly_mod_copy (h, g, mod);
        zp_dec (g[1], p);           // g = x^(p^k)-x mod f
        dg = zp_poly_degree (g, d-1);
        zp_poly_gcd (g, &dg, g, dg, f, d, p);
        if ( dg > 0 ) {
            assert (!(dg%k));           
            zp_poly_find_equal_degree_factor (r, g, dg, k, p);
            zp_set_one (r[k]);
            goto done;
        }       
    }
    
    // if we get here then f is irreducible
    zp_poly_copy (r, f, d);  k = d;
    
done:
    zp_poly_mod_clear (mod);
    zp_poly_free (g, d);  zp_poly_free (h, d);
    return k;
}

// computes the Jacobi symbol of two polynomials as defined in Bach-Shallit p.141, using Alg "Poly-Jacobi" on p. 142
int zp_poly_jacobi_symbol (zp_t u[], int du, zp_t v[], int dv, zp_p p)
{
    zp_t *w;
    zp_t c;
    int s, dw;

    assert ( dv > 0 && zp_poly_is_monic (v, dv) );
    w = zp_poly_alloc (dv-1, p);
    dw = zp_poly_mod (w, u, du, v, dv, p);
    if ( dw < 0 ) { zp_poly_free (w, dv-1); return 0; }
    if ( !(dv&1) ) s = 1; else s = zp_jacobi_symbol (w[dw], p);             // s = jacobi(a/p)^dv
    if ( ! dw ) { zp_poly_free (w, dv-1); return s; }
    zp_init (c, p);
    zp_poly_make_monic (w, dw, c, p);
    zp_clear (c);
    if ( mpz_congruent_ui_p (p, 3,4) && (dw&1) && (dv&1) ) s = -s;
    s *= zp_poly_jacobi_symbol (v, dv, w, dw, p);
    zp_poly_free (w, dv-1);
    return s;
}


#define ZP_POLY_PRODUCT_TREE_MAX_LEVELS     18
#define ZP_POLY_PRODUCT_TREE_BASE_CUTOFF    8
#define ZP_POLY_PRODUCT_TREE_MAX_DEGREE (1<<20)

struct product_tree_struct {
    int *nodes;
    int *level[ZP_POLY_PRODUCT_TREE_MAX_LEVELS];
    int level_mind[ZP_POLY_PRODUCT_TREE_MAX_LEVELS];
    int h;
    int d;
};

void zp_poly_build_product_tree (struct product_tree_struct *tree, int deg)
{
    register int i, j;

    assert ( deg <= ZP_POLY_PRODUCT_TREE_MAX_DEGREE );
    // TODO: only allocate the memeory needed to handle deg
    tree->nodes = zp_malloc ((1<<ZP_POLY_PRODUCT_TREE_MAX_LEVELS)*sizeof(*tree->nodes));
    for ( tree->level[0] = tree->nodes, i = 0 ; tree->level[i] + (1<<(i+2)) + (1<<(i+3)) <= tree->nodes+(1<<ZP_POLY_PRODUCT_TREE_MAX_LEVELS) ;  i++ ) tree->level[i+1] = tree->level[i] + (1<<(i+2));
    tree->level[0][0] = (deg+1)/2;  tree->level[0][1] = deg - tree->level[0][0];
    tree->level_mind[0] = tree->level[0][1];
    for ( i = 1 ; tree->level[i-1][0] > ZP_POLY_PRODUCT_TREE_BASE_CUTOFF ; i++ ) {
        tree->level_mind[i] = deg;
        for ( j = 0 ; j < (1<<i) ; j++ ) {
            tree->level[i][2*j] = (tree->level[i-1][j]+1)/2;  tree->level[i][2*j+1] = tree->level[i-1][j] - tree->level[i][2*j];
            if ( tree->level[i][2*j+1] < tree->level_mind[i] ) tree->level_mind[i] = tree->level[i][2*j+1];
        }
    }
    tree->h = i;
    tree->d = deg;
}

void zp_poly_clear_product_tree (struct product_tree_struct *tree)
    { zp_free (tree->nodes, (1<<ZP_POLY_PRODUCT_TREE_MAX_LEVELS)*sizeof(*tree->nodes)); }

// f = \prod_i  (x-r[i]), f and r may point to the same location
void zp_poly_from_roots (zp_t f[], zp_t r[], int d, zp_p p)
{
    struct product_tree_struct tree[1];
    mpz_t *s;
    register int i, j, k, m, n, d1, d2;

    if ( d <= ZP_POLY_PRODUCT_TREE_BASE_CUTOFF ) { zp_poly_from_roots_naive (f, r, d, p); return; }
    zp_poly_build_product_tree (tree, d);

    // we could do this all in place (i.e. in f) if we knew f was big enough
    m = d+1+(2*d)/ZP_POLY_PRODUCT_TREE_BASE_CUTOFF;
    s = zp_array_alloc (m, zp_bits(p));
    i = tree->h-1;
    
    // note that we don't take advantage of the fact that we are dealing with monic polys
    for ( j = 0, k = 0 ; j < (1<<(i+1)) ; k += n, j++ )  {
        n = tree->level[i][j];
        zp_poly_from_roots_naive (s+k+j,r+k,n,p);
    }
    for ( ; i >= 0 ; i-- ) {
        n = tree->level_mind[i];
        for ( j = 0, k=0 ; j < (1<<i) ; k += d1+d2, j++ ) {
            d1 = tree->level[i][2*j];  d2 = tree->level[i][2*j+1];
            zp_array_mul_mod_xn (s+k+j, s+k+2*j, d1+1, s+k+2*j+d1+1, d2+1, 0, d1+d2+2, p);
        }
    }
    zp_poly_copy (f, s, d);
    zp_array_free (s, m);
    zp_poly_clear_product_tree (tree);
}

// attempts to compute h = sqrt(f) mod g using slow Tonelli-Shanks, where g is assumed to be irreducible
// h and f may coincide
int zp_poly_mod_sqrt (zp_t h[], zp_t f[], zp_poly_mod_t mod)
{
    zp_t *a, *ai, *g;
    mpz_t m, n, e, w;
    int i, s;

//printf ("computing sqrt of "); zp_poly_print (f, mod->d-1); printf (" mod "); zp_poly_print (mod->g, mod->d); puts ("");  
    if ( zp_poly_mod_is_zero (f, mod) ) { zp_poly_mod_set_zero (h, mod); return 1; }
    if ( h == f ) g = zp_poly_mod_alloc (mod); else g = h;
    mpz_init2 (n, mod->d*mpz_sizeinbase(mod->p,2));  mpz_init2 (m, mod->d*mpz_sizeinbase(mod->p,2));
    mpz_pow_ui (m, mod->p, mod->d);  mpz_sub_ui (m, m, 1);      // m = p^d-1 is the order of the multiplicative group of Fp[x]/(mod)
    mpz_div_2exp (n, m, 1);                                 // n = (p^d-1)/2
    
    // Determine whether f has a square root or not by computing f^((p^d-1)/2)
    zp_poly_mod_pow (g, f, n, mod);
    if ( ! zp_poly_mod_is_one (g, mod) ) { if ( h == f ) zp_poly_mod_free (g, mod); mpz_clear(m);  mpz_clear(n); return 0; }
    
    // Get a random non-residue a
    a = zp_poly_mod_alloc (mod);  ai = zp_poly_mod_alloc (mod);
    do {
        do { zp_poly_mod_randomize (a, mod); } while ( zp_poly_mod_is_zero (a, mod) );
        zp_poly_mod_pow (g, a, n, mod);
    } while ( zp_poly_mod_is_one (g, mod) );
//printf ("a = "); zp_poly_print (a, mod->d-1); printf (" is a random non-residue mod "); zp_poly_print (mod->g, mod->d); puts ("");
    if ( ! zp_poly_mod_inv (ai, a, mod) ) { printf ("Unable to compute inverse in zp_poly_mod_sqrt, modulus poly not irreducible!\n");  zp_poly_print (mod->g, mod->d); puts ("");  abort(); }
    
    for ( s = 1 ; !mpz_tstbit (n,0) ; s++ ) mpz_div_2exp (n,n,1);       // m = n*2^s with n odd

    mpz_init2 (e,s);  mpz_init2 (w, mpz_sizeinbase(n,2));
    for ( i = 2 ; i <= s ; i++ ) {                              // this loop computes e such that f=a^e mod the 2-Sylow
        zp_poly_mod_pow (g, ai, e, mod);
        zp_poly_mod_mul (g, g, f, mod);
        mpz_div_2exp (w, m, i);                             // w = (p^d-1)/2^i
        zp_poly_mod_pow (g, g, w, mod);
        if ( ! zp_poly_mod_is_one (g, mod) ) mpz_setbit (e, i-1);
    }
    zp_poly_mod_pow (g, ai, e, mod);
    zp_poly_mod_mul (g, g, f, mod);                         // g = f*a^-e
    mpz_add_ui (n, n, 1);  mpz_div_2exp (n, n, 1);
    zp_poly_mod_pow (g, g, n, mod);                         // g = f*a^(-e(n+1)/2) is the square-root of the odd part of f
    mpz_div_2exp (e, e, 1);
    zp_poly_mod_pow (ai, a, e, mod);                            // ai = a^(e/2) is the square root of the even part of f
    zp_poly_mod_mul (g, g, ai, mod);                            // g = f*a^(-e(n+1)/2)*a^(e/2) = f*a^(-en) is the square root of f
    zp_poly_mod_mul (ai, g, g, mod);                            // verify that g^2 = f
    assert ( zp_poly_mod_equal (ai, f, mod) );
    zp_poly_mod_free (a, mod); zp_poly_mod_free (ai, mod);
    mpz_clear (m); mpz_clear (n); mpz_clear (e); mpz_clear (w);
    if ( h == f ) { zp_poly_mod_copy (h, g, mod); zp_poly_mod_free (g, mod); }
    return 1;
}

#ifdef ZP_POLY_MOD_SQRT_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_mod_t p;
    zp_poly_mod_t mod;
    mpz_t P, *r, *f, *g;
    int *ds, *es;
    int d, i, k, n, mind, maxd, delta, pbits;
    
    if ( argc < 4 ) { printf ("zp_poly_sqrt_mod_test mindeg maxdeg pbits [iterations]\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) n = atoi(argv[4]); else n = 1;
    pbits = atoi(argv[3]);
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);

    ds = malloc (maxd*sizeof(*ds));  es = malloc (maxd*sizeof(*es));
    mpz_init2 (P,pbits);
    zp_random_prime (P, pbits);
    gmp_printf ("Using random prime %Zd\n", P);
    zp_mod_init (p, P);
    mpz_clear (P);
    f = zp_poly_alloc (maxd, p);  r = zp_poly_alloc (maxd, p);  g = zp_poly_alloc (maxd, p); 
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    if ( mind < 2 ) mind = 2;
    for ( d = mind ; d <= maxd  ; d+= delta ) {
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            // generate random monic irred poly of degree d -- this may take a while if d is big
            do {
                zp_poly_randomize (f, d-1,p);  zp_set_one (f[d]);
                k = zp_poly_factor (r, ds, es, f, d, p);
            } while ( k > 1 || es[0] > 1 );
            zp_poly_mod_init (mod, f, d, p, zp_poly_algo);
            zp_poly_mod_randomize (g, mod);
            zp_poly_mod_mul (f, g, g, mod);
            start = clock();
            if ( ! zp_poly_mod_sqrt (g, f, mod) ) { printf ("zp_poly_mod_sqrt failed to find sqrt of known square "); zp_poly_print (f, d-1); printf (" mod "); zp_poly_print (mod->g, d); puts (""); abort(); }
            end = clock();
            zp_poly_mod_mul (g, g, g, mod);
            if ( ! zp_poly_mod_equal (f, g, mod) ) { printf ("zp_poly_mod_sqrt returned invalid sqrt of "); zp_poly_print (f, d-1); printf (" mod "); zp_poly_print (mod->g, d); puts (""); abort(); }
            zp_poly_mod_clear (mod);
            tot += (double)(end-start)/CLOCKS_PER_SEC;
        }
        printf ("Computed sqrt of random poly modulo a %d degree poly in %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
    }
    free (ds); free (es);
    zp_mod_clear (p);
    zp_poly_free (f,maxd);  zp_poly_free (r, maxd);  zp_poly_free (g,maxd);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif

#ifdef ZP_POLY_MIN_DEGREE_FACTOR_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_mod_t p;
    mpz_t P, *r, *f, *g;
    int d, i, k, n, totk, mind, maxd, delta, pbits;
        
    if ( argc < 4 ) { printf ("zp_poly_min_degree_factor_test mindeg maxdeg pbits [iterations]\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) n = atoi(argv[4]); else n = 1;
    pbits = atoi(argv[3]);
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);

    mpz_init2 (P,pbits);
    zp_random_prime (P, pbits);
    gmp_printf ("Using random prime %Zd\n", P);
    zp_mod_init (p, P);
    mpz_clear (P);
    f = zp_poly_alloc (maxd, p);  r = zp_poly_alloc (maxd, p);  g = zp_poly_alloc (maxd, p); 
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    if ( mind < 2 ) mind = 2;
    for ( d = mind ; d <= maxd  ; d+= delta ) {
        tot = 0.0;
        for ( i = 0, totk=0 ; i < n ; i++ ) {
            zp_poly_randomize (f, d-1,p);  zp_set_one (f[d]);
            start = clock();
            k = zp_poly_min_degree_factor (r, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            totk += k;
            zp_poly_div (g, f, d, r, k, p);
            zp_poly_mul (g, g, d-k, r, k, p);
            if ( ! zp_poly_equal (f, d, g, d) ) { printf ("zp_poly_min_degree_factor failed to correctly factor: ");  zp_poly_print (f, d); abort(); }
        }
        printf ("Found min degree %.2f factor of degree %d random poly factorization in %.3f secs", (double)totk/n, d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
    }
    zp_mod_clear (p);
    zp_poly_free (f,maxd);  zp_poly_free (r, maxd);  zp_poly_free (g,maxd);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif


#ifdef ZP_POLY_FACTOR_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_mod_t p;
    mpz_t P, *r, *f, *g, *h;
    int *ds, *es;
    int d, i, j, k, m, n, t, mind, maxd, totd, delta, pbits, dg;
    
    if ( argc < 4 ) { printf ("zp_poly_factor_test mindeg maxdeg pbits [iterations]\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) n = atoi(argv[4]); else n = 1;
    pbits = atoi(argv[3]);
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);

    ds = malloc (maxd*sizeof(*ds));  es = malloc (maxd*sizeof(*es));
    mpz_init2 (P,pbits);
    zp_random_prime (P, pbits);
    gmp_printf ("Using random prime %Zd\n", P);
    zp_mod_init (p, P);
    mpz_clear (P);
    f = zp_poly_alloc (maxd, p);  r = zp_poly_alloc (maxd, p);  g = zp_poly_alloc (maxd, p);  h = zp_poly_alloc (maxd, p);
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    if ( mind < 2 ) mind = 2;
    for ( d = mind ; d <= maxd  ; d+= delta ) {
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            zp_poly_randomize (f, d-1,p);  zp_set_one (f[d]);
            start = clock();
            k = zp_poly_factor (r, ds, es, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            for ( j = 0, totd = 0 ; j < k ; j++ ) totd += ds[j]*es[j];
            if ( totd != d ) { printf ("Bad degree list returned by zp_poly_factor\n"); abort(); }
            zp_set_one (g[0]);  dg = 0;
            for ( j = m = 0 ; j < k ; m += ds[j++] ) {
if ( es[j] == 1 ) printf ("%d ", ds[j]); else printf ("%d^%d ", ds[j], es[j]);
                zp_poly_copy (h, r+m, ds[j]-1);  zp_set_one (h[ds[j]]);
                for ( t = 0 ; t < es[j] ; t++ ) { zp_poly_mul (g, g, dg, h, ds[j], p);  dg += ds[j]; }
            }
puts ("");
            if ( ! zp_poly_equal (f, d, g, dg) ) { printf ("zp_poly_factor failed to correctly factor: ");  zp_poly_print (f, d); abort(); }
        }
        printf ("Degree %d random poly factorization took %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
    }
    free (ds); free (es);
    zp_mod_clear (p);
    zp_poly_free (f,maxd);  zp_poly_free (r, maxd);  zp_poly_free (g,maxd);  zp_poly_free (h, maxd);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif


#ifdef ZP_POLY_FIND_ROOT_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_mod_t p;
    mpz_t r, *f;
    int d, i, k, n, mind, maxd, delta, pbits;
    
    if ( argc < 4 ) { printf ("zp_poly_find_root_test mindeg maxdeg pbits [iterations]\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) n = atoi(argv[4]); else n = 1;
    pbits = atoi(argv[3]);
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);

    mpz_init2(r,pbits);
    zp_random_prime (r, pbits);
    gmp_printf ("Using random prime %Zd\n", r);
    zp_mod_init (p, r);
    f = zp_poly_alloc (maxd+1, p);
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    if ( mind < 2 ) mind = 2;
    for ( d = mind ; d <= maxd  ; d+= delta ) {
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            zp_poly_randomize (f, d-2,p);  zp_set_one (f[d-1]); zp_random (r, p);  zp_poly_add_root (f, f, d-1, r, p);
            start = clock();
            k = zp_poly_find_root (r, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            if ( ! k ) { printf ("zp_poly_find_root failed for 1-root poly with d=%d, no root found!\n", d); abort(); }
            zp_poly_eval (r, f, d, r, p);
            if ( ! zp_is_zero(r) ) { printf ("zp_poly_find_root failed for 1-root poly, returned a non-root!\n"); abort(); }
        }
        printf ("Degree %d 1-root poly took %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            zp_poly_randomize (f, d-3,p);  mpz_set_one (f[d-2]); zp_random (r, p); zp_poly_add_root (f, f, d-2, r, p);  zp_random (r, p); zp_poly_add_root (f, f, d-1, r, p);
            start = clock();
            k = zp_poly_find_root (r, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            if ( ! k ) { printf ("zp_poly_find_root failed for 2-root poly with d=%d, found %d < 2 roots!\n", d, k); abort(); }
            zp_poly_eval (r, f, d, r, p);
            if ( ! zp_is_zero(r) ) { printf ("zp_poly_find_root failed for 2-root poly, returned a non-root!\n"); abort(); }
        }
        printf ("Degree %d 2-root poly took %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            zp_set_one(f[0]);
            for ( k = 0 ; k < d; k++ ) { zp_random (r, p); zp_poly_add_root (f, f, k, r, p); }
            start = clock();
            k = zp_poly_find_root (r, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            if ( !k ) { printf ("zp_poly_find_root failed for split poly with d=%d, found %d roots!\n", d, k); abort(); }
            zp_poly_eval (r, f, d, r, p);
            if ( ! zp_is_zero(r) ) { printf ("zp_poly_find_root failed, returned a non-root!\n"); abort(); }
        }
        printf ("Degree %d d-root poly took %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
    }
    zp_mod_clear (p);
    zp_poly_free (f,maxd+1);
    mpz_clear (r);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif


#ifdef ZP_POLY_FIND_ROOTS_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_mod_t p;
    mpz_t P, *r, *f;
    int d, i, j, k, n, mind, maxd, delta, pbits;
    
    if ( argc < 4 ) { printf ("zp_poly_find_roots_test mindeg maxdeg pbits [iterations]\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) n = atoi(argv[4]); else n = 1;
    pbits = atoi(argv[3]);
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);

    mpz_init2(P,pbits);
    zp_random_prime (P, pbits);
    gmp_printf ("Using random prime %Zd\n", P);
    zp_mod_init (p, P);
    f = zp_poly_alloc (maxd+1, p);  r = zp_poly_alloc (maxd, p);
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    if ( mind < 2 ) mind = 2;
    for ( d = mind ; d <= maxd  ; d+= delta ) {
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            zp_poly_randomize (f, d-2,p);  zp_set_one (f[d-1]); zp_random (P, p);  zp_poly_add_root (f, f, d-1, P, p);
            start = clock();
            k = zp_poly_find_roots (r, d, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            if ( ! k ) { printf ("zp_poly_find_root failed for 1-root poly with d=%d, no root found!\n", d); abort(); }
            for ( j = 0 ; j < k ; j++ ) {
                zp_poly_eval (r[j], f, d, r[j], p);
                if ( ! zp_is_zero(r[j]) ) { printf ("zp_poly_find_root failed for 1-root poly, returned a non-root!\n"); abort(); }
            }
        }
        printf ("Degree %d 1-root poly took %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            zp_poly_randomize (f, d-3,p);  mpz_set_one (f[d-2]); zp_random (P, p); zp_poly_add_root (f, f, d-2, P, p);  zp_random (P, p); zp_poly_add_root (f, f, d-1, P, p);
            start = clock();
            k = zp_poly_find_roots (r, d, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            if ( k < 2 ) { printf ("zp_poly_find_root failed for 2-root poly with d=%d, found %d < 2 roots!\n", d, k); abort(); }
            for ( j = 0 ; j < k ; j++ ) {
                zp_poly_eval (r[j], f, d, r[j], p);
                if ( ! zp_is_zero(r[j]) ) { printf ("zp_poly_find_root failed for 2-root poly, returned a non-root!\n"); abort(); }
            }
        }
        printf ("Degree %d 2-root poly took %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
        tot = 0.0;
        for ( i = 0 ; i < n ; i++ ) {
            zp_set_one(f[0]);
            for ( k = 0 ; k < d; k++ ) { if ( k < (3*d)/4 ) zp_random (P, p); zp_poly_add_root (f, f, k, P, p); }       // make 1/4 of the roots duplicates
            start = clock();
            k = zp_poly_find_roots (r, d, f, d, p);
            end = clock();
            tot += (double)(end-start)/CLOCKS_PER_SEC;
            if ( k != d ) { printf ("zp_poly_find_root failed for split poly with d=%d, found %d roots!\n", d, k); abort(); }
            for ( j = 0 ; j < k ; j++ ) {
                zp_poly_eval (r[j], f, d, r[j], p);
                if ( ! zp_is_zero(r[j]) ) { printf ("zp_poly_find_root failed for %d-root poly, returned a non-root!\n", d); abort(); }
            }
        }
        printf ("Degree %d d-root poly took %.3f secs", d, tot/n);  if ( n > 1 ) printf (" (average over %d iterations)\n", n); else puts ("");
    }
    zp_mod_clear (p);
    zp_poly_free (f,maxd+1); zp_poly_free (r,maxd);
    mpz_clear (P);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif


#ifdef ZP_POLY_FIND_SPLIT_ROOT_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_mod_t p;
    mpz_t r, *f, *g;
    int d, i, k, mind, maxd, delta, pbits;
    
    if ( argc < 4 ) { printf ("zp_poly_find_split_root_test mindeg maxdeg pbits [algo]\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) zp_poly_algo = atoi(argv[4]);
    pbits = atoi(argv[3]);  mpz_init2(r,pbits);
    zp_random_prime (r, pbits);
    zp_mod_init (p, r);
    f = zp_poly_alloc (maxd+1, p);  g = zp_poly_alloc (maxd+1, p);
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    if ( mind < 2 ) mind = 2;
    for ( d = mind ; d <= maxd  ; d+= delta ) {
        tot = 0.0;
        mpz_set_one(f[0]);
        for ( k = 0 ; k < d ; k++ ) zp_random (f[k], p);        // we assume these are distinct
        zp_poly_from_roots (f, f, d, p);
        zp_poly_copy (g, f, d);
        start = clock();
        zp_poly_find_split_root (r, f, d, p);
        end = clock();
        tot += (double)(end-start)/CLOCKS_PER_SEC;
        zp_poly_eval (r, g, d, r, p);
        if ( ! mpz_is_zero(r) ) { printf ("zp_poly_find_split_root failed, returned a non-root!\n"); abort(); }
        printf ("Degree %d split poly took %.3f secs", d, tot);
    }
    zp_poly_free (f,maxd+1);  zp_poly_free (g,maxd+1);
    zp_mod_clear (p);
    mpz_clear (r);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif

#ifdef ZP_POLY_FIND_SPLIT_ROOTS_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_mod_t p;
    mpz_t P, *r, *f, *g;
    int d, i, mind, maxd, maxk, delta, pbits;
    
    if ( argc < 4 ) { printf ("zp_poly_find_split_roots_test mindeg maxdeg pbits [maxk algo]\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) maxk = atoi (argv[4]); else maxk = 0;
    if ( argc > 5 ) zp_poly_algo = atoi(argv[5]);
    pbits = atoi(argv[3]);  mpz_init2(P,pbits);
    zp_random_prime (P, pbits);
    zp_mod_init (p, P);
    r = zp_poly_alloc (maxd, p); f = zp_poly_alloc (maxd+1, p);  g = zp_poly_alloc (maxd+1, p);
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    if ( mind < 2 ) mind = 2;
    for ( d = mind ; d <= maxd  ; d+= delta ) {
        tot = 0.0;
        mpz_set_one(f[0]);
        for ( i = 0 ; i < d ; i++ ) zp_random (f[i], p);        // we assume these are distinct
        zp_poly_from_roots (f, f, d, p);
        zp_poly_copy (g, f, d);
        start = clock();
        zp_poly_find_split_roots (r, (maxk?maxk:d), f, d, p);
        end = clock();
        tot += (double)(end-start)/CLOCKS_PER_SEC;
        for ( i = 0 ; i < maxk && i < d ; i++ ) {
            zp_poly_eval (r[i], g, d, r[i], p);
            if ( ! mpz_is_zero(r[i]) ) { printf ("zp_poly_find_split_roots failed, returned a non-root!\n"); abort(); }
        }
        printf ("Degree %d split poly took %.3f secs", d, tot);
    }
    zp_poly_free (f,maxd+1);  zp_poly_free (g,maxd+1); zp_poly_free (r, maxd);
    zp_mod_clear (p);
    mpz_clear (P);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif

#ifdef ZP_POLY_XN_MOD_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double phi_sqr, ker_sqr, phi, ker, sea, tot;
    zp_poly_mod_t mod;
    mpz_t ell, p, *f, *g, *h;
    int pbits, k, m, n, mind, maxd;

    if ( argc < 4 ) { printf ("zp_poly_xn_mod_test mindeg maxdeg pbits\n"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    pbits = atoi(argv[3]);  mpz_init2(p,pbits);
    zp_random_prime (p, pbits);
    f = zp_poly_alloc (maxd, p);  g = zp_poly_alloc (maxd+1, p); h = zp_poly_alloc(maxd, p);
    tot = 0.0;
    if ( mind <= 2 ) mind = 2;
    mpz_init (ell); mpz_set_ui (ell, mind);
    for ( mpz_nextprime(ell,ell) ; mpz_cmp_ui(ell,maxd) <= 0 ; mpz_nextprime(ell,ell) ) {
        n = mpz_get_ui(ell);
        zp_poly_randomize (f,n,p); zp_poly_randomize (g,n,p); mpz_set_ui(g[n+1],1);
        zp_poly_mod_init (mod, g, n+1, p, zp_poly_algo);
        start = clock();
        for ( k = 1 ; ; k++ ) {
            zp_poly_mod_mul (h,f,f,mod);
            end = clock();
            if ( delta_secs(start,end) >= 5 ) break;
        }
        zp_poly_mod_clear (mod);
        phi_sqr = (double)(end-start) / ((double)k*CLOCKS_PER_SEC);
        phi = pbits*phi_sqr;
        m = (n-1)/2;
        zp_poly_randomize (f,m,p); zp_poly_randomize (g,m,p); mpz_set_ui(g[m+1],1);
        zp_poly_mod_init (mod, g, m+1, p, zp_poly_algo);
        start = clock();
        for ( k = 1 ; ; k++ ) {
            zp_poly_mod_mul (h,f,f,mod);
            end = clock();
            if ( delta_secs(start,end) >= 5 ) break;
        }
        zp_poly_mod_clear (mod);
        ker_sqr = (double)(end-start) / ((double)k*CLOCKS_PER_SEC);
        ker = pbits*ker_sqr;
        sea = phi + ker;  tot += sea;
        printf ("deg %5d/%5d: %5.3f/%5.3f secs per sqrmod, x^p mod time %6.0f/%6.0f secs, total %6.0f secs, grand total %9.0f secs\n", n, m, phi_sqr, ker_sqr, phi, ker, sea, tot);
    }
    zp_poly_free (f,maxd); zp_poly_free (g, maxd+1); zp_poly_free (h, maxd);
    mpz_clear (ell);  mpz_clear (p);
    zp_mem_report(1);
}
#endif

#ifdef ZP_POLY_FIND_SPLIT_ROOT_COST

int main (int argc, char *argv[])
{
    clock_t start, end;
    double tot;
    zp_poly_mod_t mod;
    mpz_t a, p, n, *h;
    int d, maxd, maxd0, maxd1, pbits, bits, mpbits;
    
    if ( argc < 3 ) { printf ("zp_poly_find_split_root_cost deg pbits\n"); return 0; }
    maxd = atoi(argv[1]);
    if ( argc == 2 && maxd ) { printf ("please specify pbits or use deg=0 for full scan\n"); return 0; }
    pbits = atoi(argv[2]);
    if ( ! maxd ) { maxd0 = 64;  maxd1 = 65536; } else maxd0 = maxd1 = maxd;
    if ( ! pbits ) { pbits = 64; mpbits = 32768; } else mpbits = pbits;
for ( ; pbits <= mpbits ; pbits *= 2 ) {
for ( maxd = maxd0 ; maxd <= maxd1 && maxd*pbits < (1<<29) ; maxd *= 2 ) {
    mpz_init2(p,pbits);
    zp_random_prime (p, pbits);
    mpz_init_modp (a, p);
    mpz_init_set_ui (n,0xAAAAAAAAAAAAAAAAUL);       // 64 bits, half of which are 1s
    h = zp_poly_alloc (maxd+1, p);
    tot = 0.0;  bits = 0;
    do {
        for ( d = maxd ; d > 2 ; d -= d/2 + (int)sqrt(d/M_PI_2) ) {     // on average, we expect to do slightly better than a 50/50 split
            // to speed up testing, use a random poly, not necessarily a split one
            zp_poly_randomize (h, d-1, p); mpz_set_ui(h[d],1);
            zp_poly_mod_init (mod, h, d, p, ZP_POLY_MOD_ALGO_KRONECKER);
            zp_random (a, p);
            start = clock();
            zp_poly_mod_pow_xan (h, a, n, mod);
            end = clock();
            tot += (double)(end-start) / CLOCKS_PER_SEC;
            zp_poly_mod_clear (mod);
        }
        bits += 64;
    } while ( tot < 5.0 );
    printf ("Estimated cost for find split root using KRONECKER with deg=%d and pbits=%d is %.2f secs\n", maxd, pbits, (tot*pbits)/bits);
    tot = 0.0;  bits = 0;
    do {
        for ( d = maxd ; d > 2 ; d -= d/2 + (int)sqrt(d/M_PI_2) ) {     // on average, we expect to do slightly better than a 50/50 split
            // to speed up testing, use a random poly, not necessarily a split one
            zp_poly_randomize (h, d-1, p); mpz_set_ui(h[d],1);
            zp_poly_mod_init (mod, h, d, p, ZP_POLY_MOD_ALGO_BARRETT);
            zp_random (a, p);
            start = clock();
            zp_poly_mod_pow_xan (h, a, n, mod);
            end = clock();
            tot += (double)(end-start) / CLOCKS_PER_SEC;
            zp_poly_mod_clear (mod);
        }
        bits += 64;
    } while ( tot < 5.0 );
    printf ("Estimated cost for find split root using BARRETT with deg=%d and pbits=%d is %.2f secs\n", maxd, pbits, (tot*pbits)/bits);
    zp_poly_free (h,maxd+1);
    mpz_clear (a);  mpz_clear (p); mpz_clear (n);
}}
    zp_cleanup ();
    zp_mem_report (0);
}

#endif
