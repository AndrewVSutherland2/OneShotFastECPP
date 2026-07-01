#ifndef _FFPOLYMULTIPOINT_INCLUDE_
#define _FFPOLYMULTIPOINT_INCLUDE_

/*
    Copyright 2007-2016 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "ff.h"
#include "ffpoly.h"
#include "cstd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FF_POLY_TREE_LEVELS    21
#define FF_POLY_TREE_BASE       8
#define FFPI_MAX_BASE_N        96   // crossover tested on an AMD Phenom II 3.0 GHz with FF_BIG_PRIME=0 and zn_poly version 0.9 (11/3/2009)

typedef struct {
    ff_t *r,*s,*t;
    int n;                          // number of roots in r
    int k;                          // number of levels
    int *p[FF_POLY_TREE_LEVELS];    // p[i] is a partition of n into 2^i parts whose sizes differ by at most 1 (p[i][0] is min and p[i][2^i-1] is max), 0 terminated
    ff_t *f[FF_POLY_TREE_LEVELS];   // f[i] is a list of 2^i polys f_0,...f_{2^i-1} with deg f_j = p[i][j] and roots of f_j are corresponding elements of r
} ffpt_ctx_t[1];

void ffpt_create_ctx (ffpt_ctx_t ctx, ff_t *x, int n);
void ffpt_destroy_ctx (ffpt_ctx_t ctx);
void ffpt_eval (ff_t y[], ff_t f[], int d, ffpt_ctx_t ctx);

void ff_poly_multipoint_eval_small (ff_t y[], ff_t f[], int d, ff_t x[], int n);
void ff_poly_multipoint_eval_big (ff_t y[], ff_t f[], int d, ff_t x[], int n);
void ff_poly_multipoint_eval (ff_t y[], ff_t f[], int d, ff_t x[], int n);

// Given f(x) and n computes the product f(0)f(1)...f(n-1) using Bostan-Gaudry-Schost O(M(sqrt(n))) algorithm
void ff_poly_multipoint_product (ff_t y[1], ff_t f[], int d, long n);
void ff_poly_multipoint_product_naive (ff_t y[1], ff_t f[], int d, long n);

// Convert a sequence of n r-by-y matrices to an r-by-r matrix of sequences of length n
static inline void ff_matrix_sequence_to_matrix (ff_t G[], ff_t F[], int r, int n)
{
    register int i, j, k;
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < r ; j++ ) for ( k = 0 ; k < n ; k++ ) _ff_set(G[(i*r+j)*n+k],F[k*r*r+i*r+j]);
}

// Given r-by-r matrix F(x) of polynomials of degree <= d computes F(0)F(1)...F(n-1) using Bostan-Gaudry-Schost O(M(sqrt(n))) algorithm
void ff_poly_matrix_multipoint_product (ff_t y[], ff_t f[], int r, int d, long n);
void ff_poly_matrix_multipoint_product_naive (ff_t y[], ff_t f[], int r, int d, long n);
void ff_poly_matrix_multipoint_product_scale (ff_t y[], ff_t f[], int r, int d, long n, int k); // scales degree by 2^k (k can be tuned, output is unchanged)

typedef struct {
    ff_t *c, *r, *s, *t;
    int n;
    int levels;
    int lens[FF_POLY_TREE_LEVELS+1];
    int *cnts[FF_POLY_TREE_LEVELS+1];
    int *cntbase;
    ff_t *mtree[FF_POLY_TREE_LEVELS+1];
    ff_t *mbase;
} ffpi_ctx_t[1];

void ffpi_compute_si_naive (ff_t s[], ff_t x[], int n);                                 // computes s[i] = 1 / prod_{j!=i} (x[i]-x[j]), requires x[i] distinct and s and x cannot overlap
void ffpi_alloc_ctx (ffpi_ctx_t ctx, int n);
void ffpi_free_ctx (ffpi_ctx_t ctx);
void ffpi_setup_ctx (ffpi_ctx_t ctx, ff_t x[], int n, int combo_only);                  // if combo_only flag is set, only ffpi_combo may be called subsequently (no interpolation)
void _ffpi_interpolate (ff_t f[], ff_t y[], int n, int o, int combo, ffpi_ctx_t ctx);   // if combo is set, computes f(X) = sum y_i g(X)/(X-x_i), otherwise interpolates f(X) s.t. f(x_i)=y_i, computing the first o coefficients

// computes f(X) = sum c_i g(X)/(X-x_i)  where g(X) = prod (X-x_i) (the x_i are specified in ffpi_setup_ctx)
static inline void ffpi_combo (ff_t f[], ff_t c[], int n, ffpi_ctx_t ctx)
    { _ffpi_interpolate (f,c,n,n,1,ctx); }

// computes f(X) of degree less than n such that f(x_i)=y_i, overlap is ok
static inline void ffpi_interpolate(ff_t f[], ff_t y[], int n, ffpi_ctx_t ctx)
    { _ffpi_interpolate (f,y,n,n,0,ctx); }

// computes f(X) of degree less than n such that f(x_i)=y_i, overlap is ok

// computes f(X) of degree less than n such that f(x_i)=y_i, overlap is ok
void ff_poly_interpolate_r (ff_t f[], int *d, ff_t x[], ff_t y[], int n);
void ff_poly_interpolate (ff_t f[], int *d, ff_t x[], ff_t y[], int n);



#ifdef __cplusplus
}
#endif

#endif
