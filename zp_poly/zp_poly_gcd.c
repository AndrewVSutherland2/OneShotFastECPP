/*

zp_poly_gcd.c: polynomial gcd and inversion functions for zp_poly

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

// TODO: switch mpz_t to zp_t

#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"
#include "zp_poly_pair.h"

// TODO: base-case cutoff should depend on the size of p
#define XGCD_HALF_BASE      18  // call base case if maxd <= XGCD_HALF_BASE
#define XGCD_BASE               18  // call base case if maxd <= XGCD_BASE
#define INV_BASE                64  // use Euclidean part GCD if d <= INV_BASE
#define XGCD_MAX_DEPTH      32
#define XGCD_STACK_MULT     14  // allocate XGCD_STACK_MULT*(d+1)
#define XGCD_STACK_FUDGE        256  //          + XGCD_STACK_FUDGE entries
#define XGCD_SCOL_ONLY      1   // only compute S column of Q (used for inverses)
#define XGCD_NOQ                2   // don't update Q at all, just compute the gcd
#define XGCD_NOQ_IFZERO     4   // don't update Q on the last Euclidean step

void zp_poly_xgcd_euclidean (mpz_t r[], int *dr, mpz_t s[], int *ds, mpz_t t[], int *dt, mpz_t f[], int df, mpz_t g[], int dg, mpz_t p);
void zp_poly_inv_mod_verify (mpz_t h[], mpz_t f[], mpz_t g[], int d, mpz_t p);

void zp_poly_gcd (zp_t r[], int *dr, zp_t f[], int df, zp_t g[], int dg, zp_p p)
    { zp_poly_xgcd (r,dr,0,0,0,0,f,df,g,dg,p); }

// h = 1/f mod g, where f and h are filled to degree d-1 (doesn't require g to be monic but g[dg] cannot be zero)
int zp_poly_inv_mod (zp_t h[], zp_t f[], zp_t g[], int dg, zp_p p)
{
    int dr, df, dh;

    assert ( mpz_sgn(g[dg]) );
    df = zp_poly_degree(f,dg-1);
    zp_poly_xgcd (0, &dr, h, &dh, 0, 0, f, df, g, dg, p);
    zp_poly_zero_pad (h, dh+1, dg-1);
    return  ( dr == 0 ? 1 : 0 );
}

// h[i] = 1/f[i] mod g, for i in [0,n-1], where f and h are zero-filled to degree d-1
int zp_poly_mod_array_inv (zp_t *h[], zp_t *f[], int n, zp_poly_mod_t g)
{
    mpz_t **w;
    int i, e, sts;

    e = g->d-1;
    assert ( n >= 0 );
    if ( !n ) return 1;
    if ( n==1 ) return zp_poly_inv_mod (h[0], f[0], g->g, g->d, g->p);
    if ( h == f ) { w = zp_malloc (n * sizeof(mpz_t *)); for ( i = 0 ; i < n ; i++ ) w[i] = zp_poly_alloc (e, g->p); } else w = h;
    zp_poly_mod_mul (w[1], f[0], f[1], g);
    for ( i = 2 ; i < n ; i++ ) zp_poly_mod_mul (w[i], w[i-1], f[i], g);
    sts = zp_poly_inv_mod (w[0], w[n-1], g->g, g->d, g->p);
    if ( sts ) {
        for ( i = n-1 ; i > 1 ; i-- ) {
            zp_poly_mod_mul (w[i], w[i-1], w[0], g);
            zp_poly_mod_mul (w[0], w[0], f[i], g);
            zp_poly_mod_copy (h[i], w[i], g);           // trivial when w==h
        }
        zp_poly_mod_mul (w[1], f[0], w[0], g);
        zp_poly_mod_mul (w[0], w[0], f[1], g);
        zp_poly_mod_copy (h[1], w[1], g);
        zp_poly_mod_copy (h[0], w[0], g);
    }
    if ( w != h ) { for ( i = 0 ; i < n ; i++ ) zp_poly_free (w[i], e); zp_free (w, n* sizeof(mpz_t *)); }
    return sts;
}

// replaces r by r'=Mr with M diagonal and r' monic (or zero) and replaces Q by MQ (so Qx=r implies Q'x=r' for all x)
void xgcd_make_monic (zp_poly_pair_t r, zp_poly_quad_t Q, int i, mpz_t p)
{
    mpz_t c;
    
    mpz_init_modp(c, p);
    // note that zp_poly_make_monic and zp_poly_smul check for 1 (and do nothing in this case)
    if ( ! zp_poly_make_monic (r->f[i], r->d[i], c, p) ) assert(0);
    zp_poly_smul (Q[0]->f[i], Q[0]->f[i], Q[0]->d[i], c, p);    // col 0 row i
    zp_poly_smul (Q[1]->f[i], Q[1]->f[i], Q[1]->d[i], c, p);    // col 1 row i
    mpz_clear (c);
}

// replaces max(r) with max(r) mod min(r), and outputs Q s.t. Q*r_in = r_out, subject to cols
// flags specifies which columns of Q, if any, are to be updated.   w[2] are work variables
// This function does not assume the input r is monic, and even if it is, the output r need not be.
int xgcd_euclidean_step (zp_poly_pair_t r, zp_poly_quad_t Q, int flags, mpz_t p, zp_poly_stack_t stack)
{
    mpz_t *a, *b, *q, *w, *mark;
    int i, da, db;

    // for readability, set a to the larger poly and b to the smaller
    i = zp_poly_pair_maxi(r);  a= r->f[i]; da = r->d[i];  b = r->f[1-i];  db = r->d[1-i];
    if ( db < 0 ) return 0;
    mark = zp_poly_stack_mark (stack);  q = zp_poly_stack_alloc (da-db, stack);
    zp_poly_div (q, a, da, b, db, p);
    zp_poly_rem (a, a, da, b, db, q, p);  r->d[i] = zp_poly_degree(a,db-1);
    if ( (flags&XGCD_SCOL_ONLY) ) {
        w = zp_poly_stack_alloc (Q[0]->d[1-i]+da-db,stack);
        zp_poly_mul (w, q, da-db, Q[0]->f[1-i], Q[0]->d[1-i], p);
        Q[0]->d[i] = zp_poly_sub (Q[0]->f[i], Q[0]->f[i], Q[0]->d[i], w, da-db+Q[0]->d[1-i], p);
    } else if ( ! (flags&XGCD_NOQ) && ( r->d[i] >= 0 || !(flags&XGCD_NOQ_IFZERO) ) ) {
        zp_poly_quad_mul_lefte (Q, q, da-db, i, p, stack);      // this is trivial (and fast) if Q=id
    }
    zp_poly_stack_pop (stack, mark);
    return ( r->d[i] >= 0 ? 1 : 0 );
}

// compute r,s,t such that r =sf+tg with r a monic divisor of f and g and ds<max(dg,0) and dt<max(df,0)
void zp_poly_xgcd_euclidean (mpz_t r[], int *dr, mpz_t s[], int *ds, mpz_t t[], int *dt, mpz_t f[], int df, mpz_t g[], int dg, mpz_t p)
{
    zp_poly_quad_t Q;
    zp_poly_pair_t u;
    zp_poly_stack_t stack;
    int i, d, k, flags;

    // setup
    flags = XGCD_NOQ_IFZERO | ( s && !t ? XGCD_SCOL_ONLY : 0) | ( !s && !t ? XGCD_NOQ : 0 );
    df = zp_poly_degree(f,df);  dg = zp_poly_degree(g,dg);  d = ( df > dg ? df : dg );
    if ( flags&XGCD_NOQ ) k = 4; else if ( flags&XGCD_SCOL_ONLY ) k=6; else k = 8;
    zp_poly_stack_init (stack, k*(d+1)+XGCD_STACK_FUDGE, p);
    zp_poly_pair_init_on_stack (u,d,stack);
    zp_poly_copy(u->f[0],f,u->d[0]=df); zp_poly_copy(u->f[1],g,u->d[1]=dg);
    zp_poly_pair_init_on_stack (Q[0], (s?d:0), stack); zp_poly_pair_init_on_stack (Q[1], (t?d:0), stack); 
    zp_poly_quad_id (Q);
    
    // step and normalize
    while ( xgcd_euclidean_step (u, Q, flags, p, stack) );
    i = zp_poly_pair_maxi(u);  xgcd_make_monic (u, Q, i, p);

    // copy results to output
    if ( r ) zp_poly_copy (r, u->f[i], u->d[i]);
    if ( dr ) *dr = u->d[i];
    if ( s ) zp_poly_copy (s, Q[0]->f[i], Q[0]->d[i]);
    if ( ds ) *ds = Q[0]->d[i];
    if ( t ) zp_poly_copy (t, Q[1]->f[i], Q[1]->d[i]);
    if ( dt ) *dt = Q[1]->d[i];
    zp_poly_stack_clear (stack);
}

// input r need not be monic, output is Q, r is unmodified
int xgcd_half_base (zp_poly_quad_t Q, zp_poly_pair_t r, mpz_t p, zp_poly_stack_t stack)
{
    mpz_t *mark;
    zp_poly_pair_t s;
    int d, m;

    d = zp_poly_pair_maxd(r);
    m = (d+1) / 2;
//printf("half_xgcd_base start (%d,%d) m=%d\n", zp_poly_pair_mind(r), zp_poly_pair_maxd(r), m);
    zp_poly_quad_id(Q);
    if ( zp_poly_pair_mind(r) < m ) return 0;
    mark = zp_poly_stack_mark (stack);
    zp_poly_pair_copy_on_stack (s, r, stack);                           // copy r to avoid modifying it
    while ( zp_poly_pair_mind(s) >= m && xgcd_euclidean_step (s, Q, 0, p, stack) );
//printf("half_xgcd_base end (%d,%d) -> (%d,%d)\n", zp_poly_pair_mind(r), zp_poly_pair_maxd(r), zp_poly_pair_mind(s), zp_poly_pair_maxd(s));
    d -= zp_poly_pair_maxd(s);
    zp_poly_stack_pop (stack, mark);
    return d;
}

// input r need not be monic, output is Q, r is unmodified
int xgcd_half (zp_poly_quad_t Q, zp_poly_pair_t r, mpz_t p, zp_poly_stack_t stack)
{
    zp_poly_pair_t s;
    zp_poly_quad_t S;
    mpz_t *mark;
    int d, m, maxd, delta;
    
    maxd = zp_poly_pair_maxd(r);
    if ( maxd <= XGCD_HALF_BASE ) { return xgcd_half_base (Q, r, p, stack); }
    m = (maxd+1) / 2;
    if ( zp_poly_pair_mind(r) < m ) { zp_poly_quad_id(Q); return 0; }
//printf ("half_xgcd begin (%d,%d)  m=%d\n", zp_poly_pair_mind(r), zp_poly_pair_maxd(r), m);
    mark = zp_poly_stack_mark (stack);
    zp_poly_pair_init_on_stack (s, maxd, stack);    // more space than we need (3/4maxd should suffice), but unused entries cost very little
    zp_poly_pair_shift_right (s, r, m);
    delta = xgcd_half (Q, s, p, stack);
    zp_poly_quad_pair_mul_mod_xn (s, Q, r, maxd-delta+1, p, stack);
    d = zp_poly_pair_mind(s);
    if ( d < m ) goto done;
    delta += zp_poly_pair_maxd(s);
    xgcd_euclidean_step (s, Q, 0, p, stack);
    delta -= zp_poly_pair_maxd(s);
    m = 2*m - d;
    zp_poly_pair_shift_right (s, s, m);
    zp_poly_quad_init_on_stack (S, zp_poly_pair_maxd(s), stack);
    delta += xgcd_half (S, s, p, stack);
    zp_poly_quad_mul (Q, S, Q, p, stack);
done:
//printf ("half_xgcd end (%d,%d)\n", zp_poly_pair_mind(r), zp_poly_pair_maxd(r));
    zp_poly_stack_pop (stack, mark);
    return delta;
}

void zp_poly_xgcd  (mpz_t r[], int *dr, mpz_t s[], int *ds, mpz_t t[], int *dt, mpz_t f[], int df, mpz_t g[], int dg, zp_p _p)
{
    zp_p_to_mpz p=_p;   
//  clock_t start, start1, end;
    zp_poly_pair_t u;
    zp_poly_quad_t S[2*XGCD_MAX_DEPTH];
    zp_poly_stack_t stack;
    mpz_t x[3];
    int i, d, maxd, delta;
    long mem;

    if ( ! t && (df <= INV_BASE || dg <= INV_BASE) ) { zp_poly_xgcd_euclidean (r, dr, s, ds, t, dt, f, df, g, dg, p); return; }
    mem = zp_mem_used;
    // allocate and setup
    d = ( df > dg ? df : dg );

    // TODO: analyze memory usage, we shouldn't need to use this much
    zp_poly_stack_init (stack, XGCD_STACK_MULT*(d+1)+XGCD_STACK_FUDGE, p);
    zp_poly_pair_init_on_stack (u, d, stack);
    zp_poly_copy (u->f[0], f, df); u->d[0] = df;  zp_poly_copy (u->f[1], g, dg); u->d[1] = dg;
    for ( i = 0 ; i < 3 ; i++ ) mpz_init2 (x[i],2*zp_bits(p)+ZP_PBITS_FUDGE);

    // compute half_xgcd
    for ( i = 0 ; zp_poly_pair_mind(u) > XGCD_BASE ; i++ ) {
        maxd = zp_poly_pair_maxd(u);
        zp_poly_quad_init_on_stack(S[i], maxd,stack);
        delta = xgcd_half (S[i], u, p, stack);
        zp_poly_quad_pair_mul_mod_xn (u, S[i], u, maxd-delta+1, p, stack);
        if ( zp_poly_pair_gap(u) > 2 ) xgcd_euclidean_step (u, S[i], 0, p, stack);
    }
    zp_poly_quad_init_on_stack(S[i],zp_poly_pair_maxd(u),stack);
    zp_poly_quad_id (S[i]);

    // finalize xgcd 
    while ( zp_poly_pair_mind(u) >= 0 && xgcd_euclidean_step (u, S[i], 0, p, stack) );
    for ( i-- ; i >= 0 ; i-- ) zp_poly_quad_mul (S[i], S[i+1], S[i], p, stack);

    // copy out results
    i = zp_poly_pair_maxi(u);  if ( u->d[i] < 0 ) i = 1-i;
    xgcd_make_monic (u, S[0], i, p);        // normalize result
    if ( r ) zp_poly_copy (r,u->f[i],u->d[i]);
    if ( dr ) *dr = u->d[i];
    if ( s ) zp_poly_copy (s,S[0][0]->f[i],S[0][0]->d[i]);
    if ( ds ) *ds = S[0][0]->d[i];
    if ( t ) zp_poly_copy (t,S[0][1]->f[i],S[0][1]->d[i]);
    if ( dt ) *dt = S[0][1]->d[i];
    for ( i = 0 ; i < 3 ; i++ ) mpz_clear(x[i]);
    zp_poly_stack_clear (stack);
    if ( zp_mem_used != mem ) { puts ("Memory leak in xgcd!"); zp_mem_report (1); abort(); }
}

void zp_poly_xgcd_verify (mpz_t r[], int dr, mpz_t s[], int ds, mpz_t t[], int dt, mpz_t f[], int df, mpz_t g[], int dg, mpz_t p)
{
    mpz_t *w1, *w2;
    int d1;
    
    if ( dr < 0 ) { printf ("Error, xgcd returned r=0\n"); goto failure; }
    w1 = zp_poly_alloc(df+dg,p);  w2 = zp_poly_alloc(df+dg,p);
    if ( ds && ds >= dg ) { printf ("ds=%d >= dg=%d\n", ds, dg); goto failure; }
    if ( dt && dt >= df ) { printf ("dt= %d >= dg=%d\n", dt, dg); goto failure; }
    zp_poly_mul2 (w1, &d1, f, df, s, ds, g, dg, t, dt, p);
    if ( ! zp_poly_equal (r, dr, w1, d1) ) { printf ("s*f + t*g = "); zp_poly_print(w1,d1); puts ("");  goto failure; }
    zp_poly_div (w1,f,df,r,dr,p);  d1=df-dr;
    zp_poly_mul (w2,r,dr,w1,d1,p);
    if ( ! zp_poly_equal(w2,df,f,df) ) { printf ("r is not a divisor of f!\n");  goto failure; }
    zp_poly_div (w1,g,dg,r,dr,p);  d1=dg-dr;
    zp_poly_mul (w2,r,dr,w1,d1,p);
    if ( ! zp_poly_equal(w2,dg,g,dg) ) { printf ("r is not a divisor of g!\n"); goto failure; }
    zp_poly_free(w1,df+dg); zp_poly_free(w2,df+dg);
    return;
failure:
    puts ("*** zp_poly_xgcd_verification failed ***");
    printf ("f = "); zp_poly_print(f,df); puts ("");
    printf ("g = "); zp_poly_print(g,dg); puts ("");
    printf ("r = "); zp_poly_print(r,dr); puts ("");
    printf ("s = "); zp_poly_print(s,ds); puts ("");
    printf ("t = "); zp_poly_print(t,dt); puts ("");
    abort();
}

void zp_poly_gcd_verify (mpz_t r[], int dr, mpz_t f[], int df, mpz_t g[], int dg, mpz_t p)
{
    mpz_t *w1, *w2;
    int d1;

    w1 = zp_poly_alloc(df+dg,p);  w2 = zp_poly_alloc(df+dg,p);
    zp_poly_div (w1,f,df,r,dr,p);  d1=df-dr;
    zp_poly_mul (w2,r,dr,w1,d1,p);
    if ( ! zp_poly_equal(w2,df,f,df) ) { printf ("r is not a divisor of f!\n");  goto failure; }
    zp_poly_div (w1,g,dg,r,dr,p);  d1=dg-dr;
    zp_poly_mul (w2,r,dr,w1,d1,p);
    if ( ! zp_poly_equal(w2,dg,g,dg) ) { printf ("r is not a divisor of g!\n"); goto failure; }
    zp_poly_free(w1,df+dg); zp_poly_free(w2,df+dg);
    return;
failure:
    puts ("*** zp_poly_gcd verification failed ***");
    printf ("f = "); zp_poly_print(f,df); puts ("");
    printf ("g = "); zp_poly_print(g,dg); puts ("");
    printf ("r = "); zp_poly_print(r,dr); puts ("");
    abort();
}

// assumes g is monic and f and h are filled to degree d-1
void zp_poly_inv_mod_verify (mpz_t h[], mpz_t f[], mpz_t g[], int d, mpz_t p)
{
    zp_poly_mod_t mod;
    mpz_t *w1, *w2;
    
    w1 = zp_poly_alloc(d-1,p);  w2 = zp_poly_alloc(d-1,p);
    zp_poly_mod_init (mod, g, d, p, 0);
    zp_poly_mod_mul (w1, f, h, mod);
    if ( ! zp_poly_mod_is_one(w1,mod) ) goto failure;
    zp_poly_mod_clear (mod);
    zp_poly_free (w1,d-1);  zp_poly_free (w2, d-1);
    return;
failure:
    puts ("*** zp_poly_inv_mod verification failed ***");
    printf ("f = "); zp_poly_print(f,d-1); puts ("");
    printf ("g = "); zp_poly_print(g,d); puts ("");
    printf ("h = "); zp_poly_print(h,d-1); puts ("");
    printf ("f*h = "); zp_poly_print(w1,d-1); puts ("");
    abort();
}

#ifdef ZP_POLY_GCD_TEST
#include <math.h>
#define LOG2            0.69314718055994530941723212145818
#define log2(x)     (log(x)/LOG2)

static char *opstr[6] = { "xgcd", "inv", "gcd", "xgcd euclidean", "gcd euclidean", "inv4" };

int main(int argc, char *argv[])
{
    clock_t start,end;
    zp_mod_t P;
    zp_poly_mod_t mod;
    unsigned long usecs, mtime;
    mpz_t *r, *s, *t, *f, *g, *h, *v[4], *w[4], p;
    int d, dr, ds, dt, mind, maxd, pbits, n, op, sts, delta;
    register int i, j;

    if ( argc < 4 ) { puts ("zp_poly_xgcd mind maxd pbits [iterations op]\n    op=0 for xgcd, op=1 for inverse, op=2 for gcd, op=3 quad xgcd, op=4 euclidean xgcd. op=5 inv4"); return 0; }
    mind = atoi(argv[1]);  maxd = atoi(argv[2]);
    if ( argc > 4 ) n = atoi(argv[4]); else n = 1;
    if ( argc > 5 ) op = atoi(argv[5]); else op = 0;
    if ( op < 0 || op > 5 ) { puts ("Invalid op specified in arg5, must be in [0,5]"); return 0; }
    pbits = atoi(argv[3]);
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);
    mpz_init2(p,pbits);
    zp_random_prime (p,pbits);  pbits = mpz_sizeinbase (p,2);
    zp_mod_init (P,p);
    
    f = zp_poly_alloc(maxd,p);  g = zp_poly_alloc(maxd,p);  r = zp_poly_alloc(maxd,p);  s = zp_poly_alloc(maxd,p);  t = zp_poly_alloc(maxd,p);  h = zp_poly_alloc(2*maxd,p);
    if ( op == 5 ) for ( i = 0 ; i < 4 ; i++ ) { v[i] = zp_poly_alloc (maxd,p); w[i] = zp_poly_alloc (maxd,p); }
    
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD
    for ( d = mind ; d <= maxd ; d+= delta ) {
        zp_poly_randomize (f, d, p);  zp_poly_randomize (g, d, p);
        start = clock();
        if ( d >= 500 ) {
            start = clock();
            for ( i = 1 ; ; i++ ) {
                zp_poly_mul (h,f,d,g,d,p);
                if ( clock()-start > 500000 ) break;
            }
        }
        mtime = (clock()-start)/i;
        usecs = 0;
        if ( op == 5 ) { mpz_set_one(g[d]);  zp_poly_mod_init (mod, g, d, p, 0); }
        start = clock();
        for ( j = 0 ; j < n ; j++ ) {
            if ( j ) { zp_poly_randomize (f, d, p); zp_poly_randomize (g, d, p); }
            if ( op==5 ) for ( i = 0 ; i < 4 ; i++ ) zp_poly_randomize (v[i], d-1, p);
            if ( op==1 ) mpz_set_ui(g[d],1);
            sts = 0;    // shut up compiler
            start = clock();
            switch (op) {
            case 0: zp_poly_xgcd (r, &dr, s, &ds, t, &dt, f, d, g, d, p); break;
            case 1: sts = zp_poly_inv_mod (r, f,  g, d, p); break;
            case 2: zp_poly_gcd (r, &dr, f, d, g, d, P); break;
            case 3: zp_poly_xgcd_euclidean (r, &dr, s, &ds, t, &dt, f, d, g, d, p); break;
            case 4: zp_poly_xgcd_euclidean (r, &dr, 0, 0, 0, 0, f, d, g, d, p); break;
            case 5: sts = zp_poly_mod_array_inv (w, v, 4, mod); break;
            }
            end = clock();
            usecs += end-start;
            switch (op) {
            case 0: case 3: zp_poly_xgcd_verify (r, dr, s, ds, t, dt, f, d, g, d, p); break;
            case 1: if ( sts ) zp_poly_inv_mod_verify (r, f, g, d, p); break;   // no verification on error
            case 2: case 4: zp_poly_gcd_verify (r, dr, f, d, g, d, p); break;
            case 5: if ( sts ) for ( i = 0 ; i < 4 ; i++ ) zp_poly_inv_mod_verify (w[i], v[i], g, d, p); break;
            }
        }
        end = clock();
        printf ("Verified %d degree %d calls to %s, each call took %.3f secs", n, d, opstr[op], usecs/(1000000.0*n));
        if ( d >= 500 ) printf (" = %.3f M log2(d) with M=%.3f secs\n", (double)usecs / (mtime * log2(d)), (double)mtime/1000000.0); else printf ("\n");
        zp_mem_report (0);
    }
    zp_poly_free(f,maxd);  zp_poly_free(g,maxd);  zp_poly_free(r,maxd);  zp_poly_free(s,maxd);  zp_poly_free(t,maxd); zp_poly_free (h,2*maxd);
    if ( op == 5 ) { zp_poly_mod_clear (mod); for ( i = 0 ; i < 4 ; i++ ) { zp_poly_free (v[i],maxd); zp_poly_free (w[i],maxd); } }
    mpz_clear (p);
    zp_cleanup ();
    zp_mem_report (1);
}

#endif
