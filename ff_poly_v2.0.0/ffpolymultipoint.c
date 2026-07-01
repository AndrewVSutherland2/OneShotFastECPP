#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ff.h"
#include "ffpoly.h"
#include "ffpolyalloc.h"
#include "ffpolybig.h"
#include "ffpolymultipoint.h"
#include "cstd.h"

/*
    Copyright 2008-2017 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// compute c[i]=1/delta(i,d)=1/prod_{j=0..d, j!=i}(i-j)
// uses 6d+O(1) mults
void bgs_idelta (ff_t c[], int d)
{
    register ff_t x;
    register int i;
    
    _ff_set_one(c[0]); _ff_set_zero(x);
    for ( i = 1 ; i <= d ; i++ ) { _ff_inc(x); _ff_multby(c[0],x); _ff_set(c[i],x); }
    ff_parallel_invert (c, c, d+1);
    if ( d&1 ) ff_negate(c[0]);
    _ff_inc(x); ff_negate(x); // x = -d-1
    for ( i = 1 ; i <= d ; i++ ) { _ff_inc(x); _ff_multby(c[i],x); _ff_multby(c[i],c[i-1]); }
}

// given a, d, and S[i]=1/(a+d+i)
// compute c[i]=Delta(a,i,d)=prod_{j=0..d}(a+i-j)
// Not currently used, integrated into translate_values below
void bgs_Delta (ff_t c[], int d, ff_t a, ff_t *S)
{
    register ff_t x, w;
    register int i;
    
    // compute Delta(a,0,d)=a*(a-1)*...*(a-d)
    ff_falling_factorial(c,a,d+1);
    _ff_set(x,a); _ff_set(w,c[0]);
    // compute w=Delta(a,i,d)=(a+i)*S[i-1]*Delta(a,i-1,d)
    for ( i = 1 ; i <= d ; i++ ) { _ff_inc(x); _ff_multby(w,x); _ff_multby(w,S[i-1]); _ff_set(c[i],w); }
}

// given y=[f(0),f(m),f(2m),...,f(dm)] with deg f <= d and a=s/m
// computes z=[f(s),f(s+m),f(s+2m),...,f(s+dm)]
// complexity is 2*M(d log p)+O(d log p)
// requires precomputed idelta[i] = 1/delta(i,d), S[i] = 1/(a+i-d) for i=0..2d, and D0 = Delta(a,0,d) = a(a-1)...(a-d)
void ff_poly_translate_values(ff_t z[], ff_t y[], int d, ff_t a, ff_t *idelta, ff_t *S)
{
    register ff_t *P,*R,*Delta;
    register int i;
    
    // _ff_set_ui(x,m); _ff_invert(x,x); _ff_mult(a,x,s);  // a = s/m

    P = ff_poly_stack_alloc(d); R = ff_poly_stack_alloc(d); Delta = ff_poly_stack_alloc(d);
    bgs_Delta (Delta, d, a, S);
    
    // P[i] = y[i]/delta(i,d) for i in [0..d]
    for ( i = 0 ; i <= d ; i++ ) _ff_mult(P[i],idelta[i],y[i]);

    // S[i] = 1/(a+i-d) for i in [0..2d]
    
    // Q = P*S (this is where all the time is spent)
    ff_poly_middle_product (R,S,2*d,P,d);
    
    // z[i] = Q[d+i]*Delta(a,i,d), using w=Delta(a,i,d)=(a+i)*S[i-1]*Delta(a,i-1,d)
    for ( i = 0 ; i <= d ; i++ ) _ff_mult(z[i],R[i],Delta[i]);
    ff_poly_stack_pop (P);
}

typedef struct {
    void *mpctx;
    ff_t *S;
    ff_t a, *P, *Q, *D;
    int d;
} ff_poly_translate_values_ctx_t[1];

// setup for translating [f(0),f(m),f(2m),...,f(dm)] with deg f <= d by s, given a=s/m and S[i]=1/(a+i-d)
void ff_poly_translate_values_setup (ff_poly_translate_values_ctx_t ctx, ff_t a, ff_t *S, int d)
{
    ctx->d = d;
    ctx->a = a;
    ctx->P = malloc(3*(d+1)*sizeof(ff_t));
    ctx->Q = ctx->P+d+1; ctx->D = ctx->Q+d+1;
    bgs_Delta(ctx->D,d,a,S);
    if ( d <= FF_POLY_MULMID_SMALL_DEGREE ) {
        ctx->S = malloc((2*d+1)*sizeof(*S));
        memcpy (ctx->S, S, (2*d+1)*sizeof(*S));
        ctx->mpctx = 0;
    } else {
        ctx->mpctx = ff_poly_middle_product_setup1 (S,2*d,d);
        ctx->S = 0;
    }
}

// given y=[f(0),f(m),f(2m),...,f(dm)] with deg f <= d computes z=[f(s),f(s+m),f(s+2m),...,f(s+dm)]
// using idelta[i] = 1/delta(i,d) computed by bgs_delta and ctx computed by ff_poly_translate_values_setup
// note that a=s/m, and d are stored in ctx
void ff_poly_translate_values_ctx (ff_t z[], ff_t y[], ff_t *idelta, ff_poly_translate_values_ctx_t ctx)
{
    register int i;
    
    // P[i] = y[i]/delta(i,d) for i in [0..d]
    for ( i = 0 ; i <= ctx->d ; i++ ) _ff_mult(ctx->P[i],idelta[i],y[i]);

    // Q = P*S (this is where all the time is spent)
    if ( ctx->mpctx ) ff_poly_middle_product_execute1 (ctx->Q, ctx->P, ctx->mpctx);
    else ff_poly_middle_product (ctx->Q,ctx->S,2*ctx->d,ctx->P,ctx->d);
    
    // z[i] = Q[i]*Delta(a,i,d), where Q[i] is the coefficient of x^(d+i) in P*S (computed using middle product algorithm)
    for ( i = 0 ; i <= ctx->d ; i++ ) ff_mult(z[i],ctx->Q[i],ctx->D[i]);
}

void ff_poly_translate_values_clear (ff_poly_translate_values_ctx_t ctx)
{
    free (ctx->P);
    if ( ctx->mpctx ) ff_poly_middle_product_clear1 (ctx->mpctx);
    if ( ctx->S ) free (ctx->S);
}

void ff_poly_matrix_translate_values(ff_t z[], ff_t y[], int r, int d, ff_t a, ff_t *idelta, ff_t *S)
{
    ff_poly_translate_values_ctx_t ctx;
    ff_t *v, *w;
    register int i,j,k;

    v = ff_poly_stack_alloc(d); w = ff_poly_stack_alloc(d);
    ff_poly_translate_values_setup (ctx, a, S, d);
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < r ; j++ ) {
        for ( k = 0 ; k <= d ; k++ ) _ff_set(v[k],y[k*r*r+i*r+j]);
        ff_poly_translate_values_ctx (w, v, idelta, ctx);
        for ( k = 0 ; k <= d ; k++ ) _ff_set(z[k*r*r+i*r+j],w[k]);
    }
    ff_poly_translate_values_clear (ctx);
    ff_poly_stack_pop(v);
}

void ff_poly_matrix_translate_values_ctx (ff_t z[], ff_t y[], int r, ff_t *idelta, ff_poly_translate_values_ctx_t ctx)
{
    ff_t *v, *w;
    register int i,j,k;

    v = ff_poly_stack_alloc(ctx->d); w = ff_poly_stack_alloc(ctx->d);
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < r ; j++ ) {
        for ( k = 0 ; k <= ctx->d ; k++ ) _ff_set(v[k],y[k*r*r+i*r+j]);
        ff_poly_translate_values_ctx (w, v, idelta, ctx);
        for ( k = 0 ; k <= ctx->d ; k++ ) _ff_set(z[k*r*r+i*r+j],w[k]);
    }
    ff_poly_stack_pop (v);
}

void ff_poly_multipoint_product_naive (ff_t y[1], ff_t f[], int d, long n)
{
    ff_t t, x, z;
    register long i;
    
    _ff_set_one(t); _ff_set_zero(x);
    for ( i = 0 ; i < n ; i++) { ff_poly_eval (&z, f, d, &x); _ff_multby(t,z); _ff_inc(x); }
    _ff_set(y[0],t);
}

// Given r-by-r matrix f(x) of polys of degree d returns f(0)f(1)...f(n-1)
void ff_poly_matrix_multipoint_product_naive (ff_t y[], ff_t f[], int r, int d, long n)
{
    ff_t x, *t,*z;
    register long i;
    
    t = ff_poly_stack_alloc(r*r-1); z = ff_poly_stack_alloc(r*r-1);
    ff_matrix_set_one(y,r); _ff_set_zero(x);
    for ( i = 0 ; i < n ; i++) { ff_poly_matrix_eval (z,f,r,d,&x); ff_matrix_mult (t,y,z,r); ff_matrix_set(y,t,r); _ff_inc(x); }
    ff_poly_stack_pop(t);
}

// Given r-by-r matrix f(x) of polys of degree d returns f(0)f(1)...f(n-1)
// Scales f(x) to degree 2^k*d (here k is a simply a parameter one can optimize, it does not change the output)
void ff_poly_matrix_multipoint_product_scale (ff_t y[], ff_t f[], int r, int d, long n, int k)
{
    ff_t a, t, *g, *h, *w;
    int dd;
    
    if ( k+ui_len(d) > 31 ) k = 31-ui_len(d);
    if ( k <= 0 ) { ff_poly_matrix_multipoint_product (y,f,r,d,n); return; }
    while ( !(n>>k) ) k--;
    g = malloc(r*r*((d<<k)+1)*sizeof(*g));  h = malloc(r*r*((d<<k)+1)*sizeof(*h)); w = malloc(r*r*((d<<k)+1)*sizeof(*h));
    memcpy (h,f,r*r*(d+1)*sizeof(*f)); dd = d;
    _ff_set_one(a); _ff_set_one(t); _ff_x2(t);
    for ( int i = 0 ; i < k ; i++ ) {
        ff_poly_matrix_translate(g,h,r,dd,a);
        ff_poly_matrix_mult (w,h,g,r,dd); dd <<= 1;
        ff_poly_matrix_scale (h,w,r,dd,t);
    }
    long m = n >> k;
    ff_poly_matrix_multipoint_product (y,h,r,dd,m);
    m <<= k;
    if ( m < n ) {
        _ff_set_ui(a,m);
        ff_poly_matrix_translate (g,f,r,d,a);
        ff_poly_matrix_multipoint_product_naive(w,g,r,d,n-m);
        ff_matrix_mult (g,y,w,r);
        ff_matrix_set (y,g,r);
    }
    free (g); free (h); free (w);
}


// Given f(x) and n computes the product f(0)f(1)...f(n-1) using Bostan-Gaudry-Schost O(M(sqrt(n))) algorithm
void ff_poly_multipoint_product (ff_t y[1], ff_t f[], int d, long n)
{
    ff_poly_translate_values_ctx_t ctx;
    ff_t g[d+1];
    ff_t a, t, z, *v, *w, *idelta, *S, *T, ee, mm;
    register int i,e,k;
    register long m;
    int dh;
    
    if ( d < 0 ) { if ( n ) _ff_set_zero(y[0]); else _ff_set_one(y[0]); return; }
    if ( !d ) { ff_exp_ui(y,f,n); return; }
    if ( n*(d+1) < 5500 || n <= (d+1) ) { ff_poly_multipoint_product_naive(y,f,d,n); return; }
    if ( d <= 16 ) {
        ff_t g[33], h[33];
        ff_poly_copy(h,&dh,f,d);
        _ff_set_one(a); _ff_set_one(t); _ff_x2(t);
        for ( e = 1 ; dh <= 16 ; e*= 2 ) {
            ff_poly_translate(g,h,dh,a);
            ff_poly_mult (h,&dh,h,dh,g,dh);
            ff_poly_scale (h,h,dh,t);
        }
        m = n/e;
        ff_poly_multipoint_product(y,h,dh,m);
        if ( e*m < n ) {
            _ff_set_ui(a,e*m);
            ff_poly_translate(g,f,d,a);
            ff_poly_multipoint_product_naive(&t,g,d,n-e*m);
            _ff_multby(y[0],t);
        }
        return;
    }
    e = d+1;
    k = (int)floor(log2(sqrt((double)n/e)));
    m = 1<<k;
    v = malloc(8*e*m*sizeof(ff_t)); w = v + 4*e*m;
    _ff_set_ui (mm,m); _ff_set_zero(z);
    for ( i = 0 ; i < e ; i++ ) { ff_poly_eval(v+i,f,d,&z); _ff_addto(z,mm); }
    
    // precompute data used/reused by ff_poly_translate
    idelta = malloc(6*e*m*sizeof(ff_t)); S = idelta+2*e*m;  T = S+2*e*m;
    _ff_set_one(t);
    ff_fraction_sequence (S,t,2*e*m);

    // a runs through 1/m,2/m,4/m,...,m/m=1, where m=2^k
    // e runs through e,2e,4e,...,m*e, mirrored by ff_t ee
    _ff_invert(a,mm); _ff_set_ui(ee,e);
    while ( ! _ff_one(a) ) {
        bgs_idelta(idelta,e-1);
        
        // compute w=[f(s),f(s+m),...,f(s+(e-1)m)] where s=ma is running through 1,2,4,8,...
        _ff_sub(t,a,ee); _ff_inc(t);
        ff_fraction_sequence(T,t,2*e-1);
        ff_poly_translate_values(w,v,e-1,a,idelta,T);
        _ff_x2(a);
        
        // extend v=[f(0),f(m),...,f((e-1)m)] to [f(0).f(m),...,f((2e-1)m)] by shifting by e*m
        // extend w=[f(s),f(s+m),...,f(s+(e-1)m)] to [f(s).f(s+m),...,f(s+(2e-1)m)] by shifting by e*m
        ff_poly_translate_values_setup (ctx,ee,S,e-1);
        ff_poly_translate_values_ctx (v+e,v,idelta,ctx);
        ff_poly_translate_values_ctx (w+e,w,idelta,ctx);
        ff_poly_translate_values_clear (ctx);
        e += e; _ff_x2(ee);

        // replace v with pointwise product of v and w (i.e. replace f(x) with f(x)f(x+a)
        for ( i = 0 ; i < e ; i++ ) _ff_multby(v[i],w[i]);
    }
    if ( e*m < n ) {
        bgs_idelta(idelta,e-1);
        ff_poly_translate_values_setup (ctx,ee,S,e-1);
        for ( i = 1 ; i*e*m < n ; i++ ) ff_poly_translate_values_ctx (v+i*e,v+(i-1)*e,idelta,ctx);
        ff_poly_translate_values_clear (ctx);
    }
    free(idelta);
    e = n/m;
    ff_product(&t,v,e);
    free(v);
    if ( m*e == n ) { _ff_set(y[0],t); return; }
    _ff_set_ui(a,m*e);
    ff_poly_translate(g,f,d,a);
    ff_poly_multipoint_product (y, g, d, n-m*e);
    _ff_multby(y[0],t);
    return;
}


// Given r-by-r matrix f(x) of polynomials of degree <= d computes f(0)f(1)...f(n-1) using Bostan-Gaudry-Schost O(M(sqrt(n))) algorithm
void ff_poly_matrix_multipoint_product (ff_t y[], ff_t f[], int r, int d, long n)
{
    ff_poly_translate_values_ctx_t ctx;
    ff_t a, t, z, *v, *w, *idelta, *A, *B, *S, *T, ee, mm;
    register int k,e;
    register long m;

    if ( d < 0 ) { if ( n ) ff_matrix_set_zero(y,r); else ff_matrix_set_one(y,r); return; }
    e = d+1;
    if ( n < 2*(d+1) ) { ff_poly_matrix_multipoint_product_naive (y,f,r,d,n); return; }
    k = (int)floor(log2(sqrt((double)n/e)));
    m = 1<<k;
    assert (e*m <= n);
    v = malloc((8*e*m+2)*r*r*sizeof(ff_t)); w = v + 4*e*m*r*r;
    A = w + 4*e*m*r*r; B = A+r*r;
    _ff_set_ui (mm,m); _ff_set_zero(z);
    for ( int i = 0 ; i < e ; i++ ) { ff_poly_matrix_eval(v+i*r*r,f,r,d,&z); _ff_addto(z,mm); }

    // precompute data used/reused by ff_poly_translate
    idelta = malloc(6*e*m*sizeof(ff_t));
    S = idelta + 2*e*m;  T = S + 2*e*m;
    _ff_set_one(t);
    ff_fraction_sequence (S,t,2*e*m);
    // a runs through 1/m,2/m,4/m,...,m/m=1, where m=2^k
    // e runs through e,2e,4e,...,m*e, mirrored by ff_t ee
    _ff_invert(a,mm); _ff_set_ui(ee,e);
    while ( ! _ff_one(a) ) {
        bgs_idelta(idelta,e-1);

        // compute w=[f(s),f(s+m),...,f(s+(e-1)m)] where s=m*a
        _ff_sub(t,a,ee); _ff_inc(t);
        ff_fraction_sequence(T,t,2*e-1);
        ff_poly_matrix_translate_values(w,v,r,e-1,a,idelta,T);
        _ff_x2(a);

        // extend v=[f(0),f(m),...,f((e-1)m)] to [f(0).f(m),...,f((2e-1)m)] by shifting by e*m
        ff_poly_translate_values_setup (ctx,ee,S,e-1);
        ff_poly_matrix_translate_values_ctx (v+e*r*r,v,r,idelta,ctx);
        ff_poly_matrix_translate_values_ctx (w+e*r*r,w,r,idelta,ctx);
        ff_poly_translate_values_clear (ctx);
        e += e; _ff_x2(ee);

        // replace v with pointwise product of v and w (i.e. replace f(x) with f(x)f(x+a)
        for ( int i = 0 ; i < e ; i++ ) { ff_matrix_mult(A,v+i*r*r,w+i*r*r,r); ff_matrix_set(v+i*r*r,A,r); }
    }
    if ( e*m < n ) {
        bgs_idelta(idelta,e-1);
        ff_poly_translate_values_setup (ctx,ee,S,e-1);
        for ( int i = 1 ; i*e*m < n ; i++ ) ff_poly_matrix_translate_values_ctx (v+i*e*r*r,v+(i-1)*e*r*r,r,idelta,ctx);
        ff_poly_translate_values_clear (ctx);
    }
    free(idelta);
    e = n/m;
    ff_matrix_product(A,v,r,e);
    if ( m*e == n ) {
        ff_matrix_set(y,A,r);
    } else {
        _ff_set_ui(a,m*e);
        ff_poly_matrix_translate(v,f,r,d,a);
        ff_poly_matrix_multipoint_product (B,v,r,d,n-m*e);
        ff_matrix_mult(y,A,B,r);
    }
    free(v);
    return;
}


// TODO: save modulus ctx in product tree to speed up repeated interpolations at the same points (should give a 2x speedup)

void ffpt_create_ctx (ffpt_ctx_t ctx, ff_t *r, int n)
{
    register ff_t *f, *g;
    register int *p;
    register int i, j;
    int d;
    
    ctx->n = n;
    ctx->r = malloc(n*sizeof(*r));
    ctx->s = malloc(2*(n+1)*sizeof(*ctx->s));
    ctx->t = ctx->s+n+1;
    memcpy (ctx->r, r, n*sizeof(*r));
    ctx->p[0] = p = malloc((2*n+FF_POLY_TREE_LEVELS)*sizeof(*p));
    *p++ = n; *p++ = 0;
    for ( i = 0 ; ; i++ ) {
        if ( *(p-2) <= FF_POLY_TREE_BASE ) break;
        assert (i+1 < FF_POLY_TREE_LEVELS);
        ctx->p[i+1] = p;
        for ( j = 0 ; j < (1<<i) ; j++ ) { *p++ = (ctx->p[i][j])/2; *p++ = (ctx->p[i][j]+1)/2; }
        *p++ = 0;
    }
    ctx->k = i; // index of bottom level, not number of levels
    i = ctx->k;
    ctx->f[i] = f = malloc((n*(i+1)+(1<<(i+1)))*sizeof(*f));
    for ( j = 0, p = ctx->p[i] ; *p ; p++ ) {
        ff_poly_from_roots (f, r+j, *p);
        f+=*p+1; j+=*p;
    }
    for ( i-- ; i >= 0 ; i-- ) {
        ctx->f[i] = f;
        for ( p = ctx->p[i+1], g = ctx->f[i+1] ; *p ; p+=2 ) {
            ff_poly_mult (f, &d, g, *p, g+*p+1, *(p+1));
            f+=*p+*(p+1)+1; g+=*p+*(p+1)+2;
        }
    }
    ff_poly_from_roots(ctx->t, ctx->r, ctx->n);
    assert ( ff_poly_equal (ctx->f[0],n,ctx->t,ctx->n) );
}

void ffpt_destroy_ctx (ffpt_ctx_t ctx)
{
    free (ctx->r);
    free (ctx->s);
    free (ctx->p[0]);
    free (ctx->f[ctx->k]);
}

void ffpt_eval (ff_t y[], ff_t h[], int d, ffpt_ctx_t ctx)
{
    register ff_t *f, *g, *r;
    register int dg, *p;
    register int i;
    
    ff_poly_mod_zpad (ctx->t, h, d, ctx->f[0], ctx->n);
    
    for ( i = 1 ; i <= ctx->k ; i++ ) {
        for ( g = ctx->t, p = ctx->p[i], f = ctx->f[i] ; *p ; p += 2 ) {
            dg = *p+*(p+1)-1;
            ff_poly_mod_zpad (ctx->s, g, dg, f, *p);  f += *p+1;
            ff_poly_mod_zpad (ctx->s+*p, g, dg, f, *(p+1));  f += *(p+1)+1;
            memcpy (g, ctx->s, (*p+*(p+1))*sizeof(*g));
            g += *p+*(p+1);
        }
    }
    for ( p = ctx->p[ctx->k], r = ctx->r, f = ctx->t ; *p ; p++ ) {
        ff_poly_multipoint_eval_small(y,f,*p-1,r,*p);
        y += *p;  r += *p;  f += *p;
    }
}

void ff_poly_multipoint_eval_small (ff_t y[], ff_t f[], int d, ff_t x[], int n)
    { for ( register int i = 0 ; i < n ; i++ ) ff_poly_eval(y+i,f,d,x+i); }

void ff_poly_multipoint_eval_big (ff_t y[], ff_t f[], int d, ff_t x[], int n)
{
    ffpt_ctx_t ctx;
    
    ffpt_create_ctx (ctx, x, n);
    ffpt_eval (y, f, d, ctx);
    ffpt_destroy_ctx (ctx);
}

void ff_poly_multipoint_eval (ff_t y[], ff_t f[], int d, ff_t x[], int n)
{
    if ( d <= 20 || n <= 1 ) { ff_poly_multipoint_eval_small (y,f,d,x,n); return; }
    if ( d <= 40 && n <= 2 ) { ff_poly_multipoint_eval_small (y,f,d,x,n); return; }
    ff_poly_multipoint_eval_big (y,f,d,x,n);
}


void ffpi_alloc_ctx (ffpi_ctx_t ctx, int n)
{
    int *cntp;
    register int i, j;

    // allocate one extra entry in each work variable, just in case we want to store a poly of degree n
    ctx->c = malloc((6*n+10)*sizeof(ff_t));
    ctx->r = ctx->c + n+1;
    ctx->s = ctx->r + n+1;
    ctx->t = ctx->s + n+1;
    ctx->mbase = malloc (n*FFPI_MAX_BASE_N*sizeof(ff_t));
    ctx->cntbase = malloc(n*sizeof(int));
    ctx->n = n;
    cntp = ctx->cntbase;
    for ( i = 0, j = n ; j > FFPI_MAX_BASE_N ; i++, j = j-j/2 );
    if ( i > FF_POLY_TREE_LEVELS ) { err_printf ("Exceeded FF_POLY_TREE_LEVELS=%d with n=%d in ffpi_alloc_ctx\n", FF_POLY_TREE_LEVELS, n); abort(); }
    ctx->levels = i;  ctx->cnts[i] = cntp++;  ctx->lens[i] = 1;  ctx->cnts[i][0] = n;
    for ( i-- ; i >= 0 ; i-- ) {
        ctx->cnts[i] = cntp;  ctx->lens[i] = 2*ctx->lens[i+1]; cntp += ctx->lens[i];
        for ( j = 0 ; j < ctx->lens[i+1] ; j++ ) {
            ctx->cnts[i][2*j] = ctx->cnts[i+1][j]/2;
            ctx->cnts[i][2*j+1] = ctx->cnts[i+1][j] - ctx->cnts[i+1][j]/2;
        }
    }
    // unroll the first iteration of the loop below to work around a bizarre compiler bug
    ctx->mtree[0] = malloc ((n+ctx->lens[0])*sizeof(ff_t));
    for ( i = 1 ; i < ctx->levels ; i++ ) ctx->mtree[i] = malloc ((n+ctx->lens[i])*sizeof(ff_t));
}

void ffpi_free_ctx (ffpi_ctx_t ctx)
{
    free (ctx->c);
    free (ctx->mbase);
    free (ctx->cntbase);
    for ( int i = 0 ; i < ctx->levels ; i++ ) free(ctx->mtree[i]);
}

// naively computes s[i] = 1 / prod_{j!=i} (x[i]-x[j]}, where the x[i] must be distinct (not verified)  s and x cannot overlap!
void ffpi_compute_si_naive (ff_t s[], ff_t x[], int n)
{
    register ff_t t0, t1;
    register int i, j;
    
    for ( i = 0 ; i < n ; i++ ) {
        _ff_set_one(t1);
        for ( j = 0 ; j < i ; j++ ) { _ff_sub(t0,x[i],x[j]);  ff_mult(t1,t1,t0); }
        for ( j++ ; j < n ; j++ ) { _ff_sub(t0,x[i],x[j]);  ff_mult(t1,t1,t0); }
        _ff_set(s[i],t1);
    }
    ff_parallel_invert(s,s,n);
}

/*
    This is a simple non-recursive version of Steps 1 and 2 of Alg. 10.11 in von zur Gathen and Gerhard.
    For the moment we just use a naive O(n^2) approach, since we expect ffpi_interpolate to be called much more often (e.g. O(n) times)
    We may want to change this in the future if it will be used in places where only one or a few interpolations use the same abcissa
*/
void ffpi_setup_ctx (ffpi_ctx_t ctx, ff_t x[], int n, int combo_only)
{
    register ff_t *b, *m, *xx;
    register int i, j, k, u, v;
    int d;

    if ( n != ctx->n ) { printf ("Error, inconsistent n=%d != %d in ffpi_ctx\n", n, ctx->n); abort(); }
    xx = x;  m = ctx->mtree[0]; b = ctx->mbase;
    for ( u = 0 ; u < ctx->lens[0] ; u++ ) {
        k = ctx->cnts[0][u];
        // compute m = prod_i (x-x[i])
        ff_poly_from_roots (m,xx,k);
        /*
            Compute polys m_i=m/(x-[i]) and transpose coefficients so that the first n entries in mbase are the constant coefficients, the next n entries the x coefficients, etc...
            Store temporarily in ctx->c.  This will make forming linear combinations of the m_i via dot-products quick and cache-friendly
        */
        for ( i = 0 ; i < k ; i++ ) {
            ff_poly_remove_root(ctx->c,m,k,xx++);
            for ( j = 0 ; j < k ; j++ ) _ff_set(b[j*k+i],ctx->c[j]);
        }
        m+=(k+1); b += k*k;
    }
    for ( v = 1 ; v < ctx->levels ; v++ ) {
        m = ctx->mtree[v];  b = ctx->mtree[v-1];
        for ( u = 0 ; u < ctx->lens[v] ; u++ ) {
            j = ctx->cnts[v-1][2*u];  k = ctx->cnts[v-1][2*u+1];
            if ( ctx->cnts[v][u] != j+k ) { printf("Inconsistent level cnts in ffpi_ctx %d != %d!\n", ctx->cnts[v][u], j+k); abort(); }
            ff_poly_mult (m, &d, b, j, b+j+1, k);
            b += j+k+2;
            m += j+k+1;
        }
    }
    if ( combo_only ) return;
    if ( ctx->levels ) {
        ff_poly_from_roots (ctx->s, x, n);
        ff_poly_derivative (ctx->c, 0, ctx->s, n);
    } else {
        ff_poly_derivative (ctx->c, 0, ctx->mtree[0], n);
    }
    ff_poly_multipoint_eval (ctx->s, ctx->c, n-1, x, n);
    ff_parallel_invert (ctx->s, ctx->s, n);
}

static inline int imin(int a, int b) { return (a < b ? a : b); }

// o specifies the number of output coefficients, which may be less than the number of input coordinates (yielding a slight speedup)
void _ffpi_interpolate (ff_t f[], ff_t y[], int n, int o, int combo, ffpi_ctx_t ctx)
{
    ff_t *c, *r, *b, *m;
    register int i,j,k,u,v;
    int d;
    
    if ( n != ctx->n ) { printf ("Error, inconsistent n in ffpi_ctx\n"); abort(); }
    if ( combo ) {
        memcpy (ctx->c, y, n*sizeof(y[0]));
    } else {
        for ( i = 0 ; i < n ; i++ ) _ff_mult(ctx->c[i],ctx->s[i],y[i]);
    }

    c = ctx->c;  r = ctx->r;  b = ctx->mbase;
    for ( u = 0 ; u < ctx->lens[0] ; u++ ) {
        k = ctx->cnts[0][u];
        for ( i = 0 ; i < k ; i++, b+=k ) ff_dot_product(r+i,b,c,k);
        c += k; r+= k;
    }
    if ( ctx->levels == 0 ) {
        memcpy (f, ctx->r, o*sizeof(*f));
        return;
    }
    for ( v = 1 ; v <= ctx->levels ; v++ ) {
        m = ctx->mtree[v-1];  r = ctx->r;
        for ( u = 0 ; u < ctx->lens[v] ; u++ ) {
            j = ctx->cnts[v-1][2*u];  k = ctx->cnts[v-1][2*u+1];
            ff_poly_mult (ctx->t, &d, m+j+1, imin(k,o), r, imin(j-1,o));
            ff_poly_mult (ctx->t+j+k, &d, m, imin(j,o), r+j, imin(k-1,o));
            for ( i = 0 ; i < j+k ; i++ ) _ff_add(r[i],ctx->t[i],ctx->t[j+k+i]);
            m += j+k+2;
            r += j+k;
        }
    }
    memcpy (f, ctx->r, o*sizeof(*f));
}

void ff_poly_interpolate (ff_t f[], int *d, ff_t x[], ff_t y[], int n)
{
    ffpi_ctx_t ctx;
    
    ffpi_alloc_ctx (ctx, n);
    ffpi_setup_ctx (ctx, x, n, 0);
    ffpi_interpolate (f, y, n, ctx);
    ffpi_free_ctx (ctx);
    if ( d ) *d = ff_poly_degree(f,n-1);
}

