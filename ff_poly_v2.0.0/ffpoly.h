#ifndef _FFPOLY_INCLUDE
#define _FFPOLY_INCLUDE

/*
    Copyright 2007-2017 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <assert.h>
#include <gmp.h>
#include "ff.h"
#include "cstd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FF_POLY_VERSION_STRING      "2.0.0"

typedef struct ff_poly_modulus_struct {
    ff_t *g, *rgi, *w, *w1;                 // rgi is the inverse of rev(g) mod x^e
    int d, e;                               // d+e is the max degree of a poly that may be reduced
} ff_poly_modulus_t[1];

// These crossovers are quite crude, the optimal value really depends on the size of p (smaller p means smaller crossover, because Kronecker wins sooner)

#define FF_POLY_MULT_SMALL_DEGREE       20      // crossover degree for big mults
#define FF_POLY_COMPOSE_SMALL_DEGREE    100     // crossover degree for modular composition, actually depends on p (this values is for p ~ 2^30
#define FF_POLY_MULMID_SMALL_DEGREE     30      // crossover degree for big middle products
#define FF_POLY_SQUARE_SMALL_DEGREE     30      // crossover degree for big squares
#define FF_POLY_MOD_SMALL_DEGREE        780     // crossover degree for big mods (setup+reduce)
#define FF_POLY_TRANSLATE_SMALL_DEGREE  30      // crossover degree to switch from naive to recursive poly translates
#define FF_POLY_TRANSLATE_BIG_DEGREE    84      // crossover degree to switch from recursive to convolution translates
#define FF_POLY_FROMROOTS_SMALL_DEGREE  64      // crossover degree to switch from inlines to recursive
#define FF_POLY_XNMOD_SMALL_DEGREE      255     // crossover degree to switch to using modulus when computing x^n mod g(x)

#define FF_POLY_SPLIT                   1       // indicates that poly is a product of linear factors over Fp
#define FF_POLY_ONE_ROOT                2       // asks for just one root (even if there are more)
#define FF_POLY_EXACTLY_ONE_ROOT        4       // used in ff_poly_roots_d3 to indicate a root should be returned only when there is exactly one (but the correct count is always returned)

// handy macros - these require integer variables d_f be declared for each poly f
#define _poly_set(a,b)              ff_poly_copy(a,&d_ ## a,b,d_ ## b)
#define _poly_print(a)              ff_poly_print(a,d_ ## a)
#define _poly_neg(c,a)              ff_poly_neg(c,&d_ ## c,a,d_ ## a)
#define _poly_monic(c,a)            ff_poly_monic(c,&d_ ## c,a,d_ ## a)
#define _poly_add(c,a,b)            ff_poly_add(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_addto(c,a)            ff_poly_addto(c,&d_ ## c,a,d_ ## a)
#define _poly_sub(c,a,b)            ff_poly_sub(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_subfrom(c,a)          ff_poly_subfrom(c,&d_ ## c,a,d_ ## a)
#define _poly_mult(c,a,b)           ff_poly_mult(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_div(q,r,a,b)          ff_poly_div(q,&d_ ## q,r,&d_ ## r,a,d_ ## a,b,d_ ## b)
#define _poly_mod(c,a,b)            ff_poly_mod(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_gcd(c,a,b)            ff_poly_gcd(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_gcdext(c,u,v,a,b)     ff_poly_gcdext(c,&d_ ## c,u,&d_ ## u, v,&d_ ## v,a,d_ ## a,b,d_ ## b)
#define _poly_expmod(c,a,e,b)       ff_poly_exp_mod(c,&d_ ## c,a,&d_ ## a,e,b,&d_ ## b)

// these are for general use, see ffpolyalloc.h for fast safe-stack allocation functions
static inline ff_t *ff_poly_alloc (int d) { return malloc ((d+1)*sizeof(ff_t)); }
static inline void ff_poly_free (ff_t *f) { free (f); }

static inline void ff_poly_randomize (ff_t f[], int d) { register int i; _ff_random_nz(f[d]); for ( i = 0 ; i < d; i++ ) _ff_random(f[i]); }

static inline void ff_poly_set_mpz (ff_t f[], mpz_t F[], int d)
    { register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_set_mpz(f[i],F[i]); }

static inline ff_t *ff_poly_copy (ff_t b[], int *pd_b, ff_t a[], int d_a)
    { if ( b!=a ) memmove (b, a, (d_a+1)*sizeof(ff_t)); if ( pd_b ) *pd_b = d_a; return b; }
    
static inline ff_t *ff_poly_zpad (ff_t f[], int dold, int dnew)
    { memset (f+dold+1, 0, (dnew-dold)*sizeof(ff_t)); return f; }


// finds a rational point (x0:y0:1) on conic c0*x^2+c1*x*y+c2*x*y+c3*y^2+c4*y*z+c5*z^2 (lex monomial order)
void ff_poly_conic_rational_point (ff_t *x0, ff_t *y0, ff_t c[6]);

// given plane quartic f(x,y,z)=f0*x^4+f1*x^3*y+...+f14*z^4 and conic c(x,y,z)=c0*x^2+c1*x*y...+c5*z^2 (lex monomial order)
// that define a conic double cover [c(x,y,z), w^2*z^2-q(x,y,z)] in P^3, returns genus 3 hyperelliptic g^2=f(x) (so deg(g) = 7 or 8).
// MOVE THIS CODE TO SMALLJAC ?
int ff_poly_hyperelliptic_curve_from_conic_cover (ff_t g[9], ff_t c[6], ff_t f[15]);

void ff_poly_med_weierstrass (ff_t f[4], ff_t W[5]);
void ff_poly_short_weierstrass (ff_t f[4], ff_t W[5]);

int ff_poly_parse (ff_t f[], int maxd, char *expr);
void ff_poly_print (ff_t f[], int d_f);
int ff_poly_sprint (char *s, ff_t f[], int d_f);

void ff_poly_twist (ff_t g[], ff_t f[], int d);
static inline void ff_poly_twist_mpz (ff_t f[], mpz_t F[], int d)
{
    ff_poly_set_mpz (f, F, d);
    ff_poly_twist (f, f, d);
}

static inline int ff_poly_is_zero (ff_t f[], int d_f) { return d_f == -1 ? 1 : 0; }
static inline int ff_poly_is_one (ff_t f[], int d_f ) { return (!d_f && _ff_one(f[0])) ? 1 : 0; }

int ff_poly_irreducible (ff_t f[], int d_f, int *pnroots);              // pnroots is the number of distinct roots (optional)
int ff_poly_factorization_pattern_and_root (int counts[], ff_t f[], int d_f, ff_t *r);      // counts[i] = # of irred factors of deg i.  returns total number of irred factors, if r is non-null sets r to a root (if any exist)
static inline int ff_poly_factorization_pattern (int counts[], ff_t f[], int d_f)
    { return ff_poly_factorization_pattern_and_root (counts, f, d_f, 0); }
static inline int ff_poly_count_factors (ff_t f[], int d_f)             // total number of  irreducible factors (including multiplicity) of f, which need not be monic
    { return ff_poly_factorization_pattern_and_root (0, f, d_f, 0); }
int ff_poly_factors (ff_t r[], int n[], ff_t f[], int d_f);             // determines the irreducible monic factors of f (which need not be monic).  Returned poly's are implicitly monic (leading 1 omitted), and concatenated in r[], with degrees in n[]
int ff_poly_roots (ff_t *r, ff_t f[], int d_f);                         // f must be monic. returns all Fp-roots (not necessarily distinct)
int _ff_poly_distinct_roots (ff_t *r, ff_t f[], int d_f, int flags);    // f must be monic.  FF_POLY_SPLIT flag => f splits into linear factors.  FF_POLY_ONE_ROOT flag => only one root will be computed (if any).  Returns count of distinct roots
static inline int ff_poly_find_root (ff_t *r, ff_t f[], int d_f)        // finds one Fp-root and returns count of distinct roots, or returns 0 if f has no Fp-roots
    { return _ff_poly_distinct_roots (r,f,d_f,FF_POLY_ONE_ROOT); }
int ff_poly_distinct_roots (ff_t *r, ff_t f[], int d_f);                // f must be monic
int ff_poly_count_roots (ff_t f[], int d_f);                            // returns the number of Fp-roots of f, with multiplicity
int ff_poly_count_distinct_roots (ff_t f[], int d_f);                   // returns the number of distinct Fp-roots of f, does not require f to be monic
int ff_poly_count_distinct_factors (int counts[], ff_t f[], int d_f, int n);  // counts[i] = # distinct irred factors of deg i, for i in [1..n]
int ff_poly_gcd_xpn (ff_t g[], int *pd_g, int n, ff_t f[], int d_f, int comp);    // returns GCD(f,x^(p^n)-x)

void ff_poly_resultant (ff_t D[1], ff_t f[], int df, ff_t g[], int dg);
void ff_poly_disc (ff_t D[1], ff_t f[], int d);
void ff_poly_disc_small (ff_t D[1], ff_t f[], int d_f);                 // monic f of degree at most 4
int ff_poly_disc_legendre_small (ff_t f[], int d_f);                    // returns Legendre symbol of discriminant for monic f of degree at most 4
static inline int ff_poly_disc_legendre (ff_t f[], int d)
    { ff_t D; if ( d <= 4 && _ff_one(f[d]) ) return ff_poly_disc_legendre_small(f,d); ff_poly_disc(&D,f,d); return ( _ff_zero(D) ? 0 : ( ff_residue(D) ? 1 : -1 ) ); }
int ff_poly_squarefree (ff_t f[], int d);                               // returns 1 if discriminant is nonzero, 0 otherwise (uses gcd with derivative, does not compute the discriminant)
int ff_poly_gcd (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b);        // output g is not made monic.  returns deg
int ff_poly_gcd_reduce (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b); // output g is not made monic.  returns deg
int ff_poly_gcdext (ff_t g[], int *pd_g, ff_t u[], int *pd_u, ff_t v[], int *pd_v, ff_t a[], int d_a, ff_t b[], int d_b);   // computes monic extended gcd, returns deg
int ff_poly_inv_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g);
static inline int ff_poly_inv_modulus (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_poly_modulus_t mod)
    { return ff_poly_inv_mod (h, pd_h, f, d_f, mod->g, mod->d); }
int ff_poly_trivial_gcd (ff_t f[], int d_f, ff_t g[], int d_g);             // very fast code to test whether f and g have a common factor
void ff_poly_div (ff_t q[], int *pd_q, ff_t r[], int *pd_r, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_div_exact (ff_t q[], int *pd_q, ff_t a[], int d_a, ff_t b[], int d_b);
void  ff_poly_inverse_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n);  // given f with f(0)!=0, computes g(x) = 1/f(x) mod x^n
void  ff_poly_log_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n);      // given f with f(0)=1, computes log(f(x)) mod x^n
void  ff_poly_exp_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n);      // given f with f(0)!=0, computes exp(f(x)) mod x^n
void ff_poly_antiderivative (ff_t *g, int *d_g, ff_t *f, int d_f);
void ff_poly_fold_mod_xn (ff_t *f, int *pd_f, ff_t *a, int d_a, ff_t *b, int d_b, ff_t *c, int d_c, int n);

void ff_poly_mult_small (ff_t o[], ff_t f[], ff_t g[], int n);                  // polys must be the same size (padded to degree n, meaning n+1 coefficients, o must have space for poly of degree 2n
void ff_poly_mult_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);         // uses zn_poly FFT implmenetation
void ff_poly_mult (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b);  // uses hard-wired schoolbook for very small degrees, then ff_poly_mult_small, then ff_poly_mult_big based on FF_POLY_MULT_SMALL_DEGREE

// g must be large enough to hold a poly of degree 2d_f
void ff_poly_square_small (ff_t g[], ff_t f[], int d_f);
void ff_poly_square_big (ff_t g[], ff_t f[], int d_f);                          // uses zn_poly FFT implementation
static inline void ff_poly_square (ff_t g[], ff_t f[], int d_f)
{
    if ( d_f <= FF_POLY_SQUARE_SMALL_DEGREE ) ff_poly_square_small (g,f,d_f);
    else ff_poly_square_big (g,f,d_f);
}

// computes middle product c[i] = sum_j a[d_b+i-j]*b[j] for i in [0,d_a-d_b]
static inline void ff_poly_middle_product_small (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
    { for ( register int i = 0 ; i <= d_a-d_b ; i++ ) ff_conv (c+i, a+i, b, d_b+1); }
void ff_poly_middle_product_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);   // uses zn_poly FFT implementation
static inline void ff_poly_middle_product (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
{
    if ( d_b <= FF_POLY_MULMID_SMALL_DEGREE ) ff_poly_middle_product_small (c,a,d_a,b,d_b);
    else ff_poly_middle_product_big(c,a,d_a,b,d_b);
}

void ff_poly_from_roots_small (ff_t f[], ff_t r[], int d);              // computes f(x) = prod_i (x-r[i]) where i ranges from o to d-1 (so deg f = d)
void ff_poly_from_roots (ff_t f[], ff_t r[], int d);                    // recursive divide-and-conquer down to ff_poly_from_roots_small based on FF_POLY_FROMROOTS_SMALL_DEGREE

// y[i] = f(x[i]) for i in [0,n-1], uses O(M(n)log(n)) algorithm
void ff_poly_multipoint_eval_small (ff_t y[], ff_t f[], int d, ff_t x[], int n);
void ff_poly_multipoint_eval_big (ff_t y[], ff_t f[], int d, ff_t x[], int n);
void ff_poly_multipoint_eval (ff_t y[], ff_t f[], int d, ff_t x[], int n);

// g(x)=f(x+c), works in place, uses O(M(n) algorithm
void ff_poly_translate_small (ff_t g[], ff_t f[], int d_f, ff_t c);
void ff_poly_translate_big (ff_t g[], ff_t f[], int d_f, ff_t c);
void ff_poly_translate (ff_t g[], ff_t f[], int d_f, ff_t c);

// O(M(n)log(n)) algorithm
void ff_poly_interpolate (ff_t f[], int *d, ff_t x[], ff_t y[], int n);

// Given f(x) and n computes the product f(0)f(1)...f(n-1) using Bostan-Gaudry-Schost O(M(sqrt(n))) algorithm
void ff_poly_multipoint_product (ff_t y[1], ff_t f[], int d, long n);

// Given r-by-r matrix F(x) of polynomials of degree <= d computes F(0)F(1)...F(n-1) using Bostan-Gaudry-Schost O(M(sqrt(n))) algorithm
void ff_poly_matrix_multipoint_product (ff_t y[], ff_t f[], int r, int d, long n);
void ff_poly_matrix_multipoint_product_scale (ff_t y[], ff_t f[], int r, int d, long n, int k);

void ff_poly_mod_setup_max (ff_poly_modulus_t mod, ff_t *g, int d_g, int max);
static inline void ff_poly_mod_setup (ff_poly_modulus_t mod, ff_t *g, int d_g)
    { ff_poly_mod_setup_max (mod, g, d_g, 0); }
void ff_poly_mod_clear (ff_poly_modulus_t mod) ;
void ff_poly_mod_reduce (ff_t *h, int *d_h, ff_t *f, int d_f, ff_poly_modulus_t mod);

// computes g = a mod f
void ff_poly_mod_small (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t f[], int d_f);
void ff_poly_mod_small_inplace (ff_t *g, int d_g, ff_t *f, int d_f); // replaces g with g mod f, zero padded to degree d_f-1
void ff_poly_mod_big (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t f[], int d_f);
static inline void ff_poly_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t f[], int d_f)
{
    if ( d_f <= FF_POLY_MOD_SMALL_DEGREE )  { ff_poly_mod_small (g,pd_g,a,d_a,f,d_f); return; }
    ff_poly_mod_big (g,pd_g,a,d_a,f,d_f);
}
static inline void ff_poly_mod_zpad (ff_t h[], ff_t *f, int d_f, ff_t *g, int d_g)
    { int d_h; ff_poly_mod (h, &d_h, f, d_f, g, d_g); ff_poly_zpad (h, d_h, d_g-1); }

// g = h1(h2(x)) mod f
void ff_poly_compose_mod_small (ff_t g[], int *pd_g, ff_t h1[], int d_h1, ff_t h2[], int d_h2, ff_t f[], int d_f);
void ff_poly_compose_mod_big (ff_t g[], int *pd_g, ff_t h1[], int d_h1, ff_t h2[], int d_h2, ff_t f[], int d_f);
static inline void ff_poly_compose_mod (ff_t g[], int *pd_g, ff_t h1[], int d_h1, ff_t h2[], int d_h2, ff_t f[], int d_f)
{
    if ( d_f <= 20 || 30*d_f <= (FF_POLY_COMPOSE_SMALL_DEGREE * ui_len(_ff_p)) ) { // this cutoff is sub-optimal in many cases
        ff_poly_compose_mod_small (g, pd_g, h1, d_h1, h2, d_h2, f, d_f);
    } else {
        ff_poly_compose_mod_big (g, pd_g, h1, d_h1, h2, d_h2, f, d_f);
    }
}

void ff_poly_xn_mod (ff_t g[], int *pd_g, unsigned long e, ff_t f[], int d_f);
void ff_poly_xn_mod_mpz (ff_t g[], int *pd_g, mpz_t e, ff_t f[], int d_f);
void ff_poly_xn_modulus  (ff_t g[], int *pd_g, unsigned long e, ff_poly_modulus_t mod);
void ff_poly_xn_modulus_mpz (ff_t g[], int *pd_g, mpz_t e, ff_poly_modulus_t mod);
void ff_poly_xpan_mod (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_t f[], int d_f);
void ff_poly_xpan_modulus  (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_poly_modulus_t mod);
void ff_poly_pow_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e, ff_t f[], int d_f);
void ff_poly_pow_mod_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_t f[], int d_f);
void ff_poly_pow_modulus (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e, ff_poly_modulus_t mod);
void ff_poly_pow_modulus_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_poly_modulus_t mod);

int ff_poly_sqrt_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g);
int ff_poly_sqrt_modulus (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_poly_modulus_t mod);

void _ff_poly_depress_monic_inplace (ff_t s[1], ff_t f[], int d_f);
static inline void ff_poly_depress_monic_inplace (ff_t s[1], ff_t f[], int d_f)
    { if ( d_f > 0 && _ff_zero(f[d_f-1]) ) { if ( s ) _ff_set_zero(s[0]); return; } _ff_poly_depress_monic_inplace (s,f,d_f); }
static inline ff_t *ff_poly_depress_monic (ff_t s[1], ff_t g[], ff_t f[], int d_f)
    { ff_poly_depress_monic_inplace(s,ff_poly_copy(g,0,f,d_f),d_f); return g; }

int _ff_poly_roots_d3 (ff_t r[3], ff_t f[4], ff_t *pd, int flags);      // f must be of the form x^3+ax+b, r may be null if only a count is desired, pd is an optional pointer to sqrt(-D/3) (used for 3-torsion)
int _ff_poly_roots_d3_mod3 (ff_t r[3], ff_t f[4]);              // used only when _ff_p==3.  f must be monic
int _ff_poly_roots_d4 (ff_t r[4], ff_t f[5], ff_t *pd, int flags);      // f must be of the form x^3+ax^2+bx+c, r is required, pd is optional
int _ff_poly_roots_d4_mod3 (ff_t r[4], ff_t f[5]);              // used only when _ff_p==3.  f must be monic

// all the ff_poly_g1 functions have been moved to ecurve in smalljac
void ff_poly_xn_mod_d3 (ff_t g[3], unsigned long n, ff_t f[2]); // computes x^n mod f=x^3+ax+b
void ff_poly_xn_mod_d4 (ff_t g[4], unsigned long n, ff_t f[3]);
void ff_poly_xpan_mod_d3 (ff_t g[3], ff_t a, unsigned long n, ff_t f[2]);
void ff_poly_xpan_mod_d4 (ff_t g[4], ff_t a, unsigned long n, ff_t f[3]);

void ff_poly_xpan_mod_d2 (ff_t g[2], ff_t a, unsigned long n, ff_t f[1]);
void ff_poly_xpan_mod_d2_o (ff_t g[2], ff_t a, unsigned long n, ff_t f[1]);

static inline int ff_poly_degree (ff_t f[], int d_f) { while ( d_f>=0 && _ff_zero(f[d_f]) ) d_f--;  return d_f; }

// note that f may have leading zeroes after reversing inplace -- use ff_poly_degree to get the correct degree
static inline ff_t *ff_poly_reverse_inplace (ff_t f[], int d)
{
    register ff_t t;
    register int i,j;
    
    if ( d <= 0 ) return f;
    j = (d>>1) + (d&1);     // j = ceil(d_f/2)
    for ( i = 0 ; i < j ; i++ ) { _ff_set(t,f[i]);  _ff_set(f[i],f[d-i]);  _ff_set (f[d-i],t); }
    return f;
}

static inline ff_t *ff_poly_reverse (ff_t g[], int *pd_g, ff_t f[], int d_f)
{
    if ( g == f ) ff_poly_reverse_inplace (f, d_f); else for ( int i = 0 ; i <= d_f ; i++ ) _ff_set (g[i],f[d_f-i]);
    if ( pd_g ) *pd_g = ff_poly_degree (g, d_f);
    return f;
}

static inline int ff_poly_equal (ff_t a[], int d_a, ff_t b[], int d_b)
    { return ( d_a == d_b && ! memcmp(a,b,d_a*sizeof(d_a)) ); }

static inline int ff_poly_set_i (ff_t f[], long F[], int d) { register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_set_i(f[i],F[i]); return ff_poly_degree (f, d); }
static inline int ff_poly_set_ui (ff_t f[], unsigned long F[], int d) { register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_set_ui(f[i],F[i]); return ff_poly_degree(f, d); }

// overlap ok
static inline void ff_poly_double (ff_t b[], int *d_b, ff_t a[], int d_a)
{
    for ( int i = 0 ; i <= d_a ; i++ ) _ff_add(b[i],a[i],a[i]);
    if ( d_b ) *d_b = d_a;
}

// overlap ok
static inline void ff_poly_triple (ff_t b[], int *d_b, ff_t a[], int d_a)
{
    register ff_t t0;
    
    for ( int i = 0 ; i <= d_a ; i++ ) { _ff_add(t0,a[i],a[i]); _ff_add(b[i],a[i],t0); }
    if ( d_b ) *d_b = d_a;
}


// overlap ok
static inline void ff_poly_scalar_mult (ff_t c[], int *d_c, ff_t a, ff_t b[], int d_b)
{
    for ( int i = 0 ; i <= d_b ; i++ ) ff_mult(c[i],a,b[i]);
    if ( d_c ) *d_c = d_b;
}

// overlap ok
static inline void ff_poly_mult_xpa (ff_t g[], ff_t f[], int d, ff_t a)
{
    _ff_set(g[d+1],f[d]);
    for ( int i = d ; ; i-- ) { _ff_mult(g[i],f[i],a); if ( ! i ) return; _ff_addto(g[i],f[i-1]); } 
}

// overlap not ok, use ff_poly_addto
static inline void ff_poly_add (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b)
{
    register int i,d_c;

    assert (c!=a && c != b);
    if ( d_a < 0 ) { ff_poly_copy (c, pd_c, b, d_b);  return; }
    if ( d_b < 0 ) { ff_poly_copy (c, pd_c, a, d_a);  return; }
    
    d_c = d_a;
    if ( d_b > d_a ) d_c = d_b;
    for ( i = 0 ; i <= d_c ; i++ ) {
        _ff_set_zero (c[i]);
        if ( i <= d_a ) _ff_addto (c[i], a[i]);
        if ( i <= d_b ) _ff_addto (c[i], b[i]);
    }
    *pd_c = ff_poly_degree(c,d_c);
}

// overlap ok
static inline void ff_poly_addto (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
    register int i;

    if ( d_b < 0 ) return;
    for ( i = 0 ; i <= d_b ; i++ ) {
        if ( i > *pd_c ) {
            _ff_set (c[i], b[i]);
            *pd_c = i;
        } else {
            _ff_addto (c[i], b[i]);
        }
    }
    *pd_c = ff_poly_degree(c,*pd_c);
}

// overlap ok
static inline ff_t *ff_poly_neg (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
    if ( pd_c ) *pd_c = d_b;
    if ( c == b ) { for ( int i = 0 ; i <= d_b ; i++ ) ff_negate(c[i]); }
    else { for ( int i = 0 ; i <= d_b ; i++ ) _ff_neg(c[i], b[i]); }
    return c;
}

// overlap not ok - use poly_subfrom instead
static inline void ff_poly_sub (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b)
{
    int d_c;
    int i;

    assert (c!=a && c !=b);
    if ( d_a < 0 ) { ff_poly_neg (c, pd_c, b, d_b);  return; }
    if ( d_b < 0 ) { ff_poly_copy (c, pd_c, a, d_a);  return; }
    
    d_c = d_a;
    if ( d_b > d_a ) d_c = d_b;
    for ( i = 0 ; i <= d_c ; i++ ) {
        _ff_set_zero (c[i]);
        if ( i <= d_a ) _ff_addto (c[i], a[i]);
        if ( i <= d_b ) _ff_subfrom (c[i], b[i]);
    }
    *pd_c = ff_poly_degree(c,d_c);
}

// overlap ok
static inline void ff_poly_subfrom (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
    register int i;

    if ( d_b < 0 ) return;
    for ( i = 0 ; i <= d_b ; i++ ) {
        if ( i > *pd_c ) {
            _ff_neg (c[i], b[i]);
            *pd_c = i;
        } else {
            _ff_subfrom (c[i], b[i]);
        }
    }
    *pd_c = ff_poly_degree(c,*pd_c);
}

// overlap ok
static inline void ff_poly_shift_up (ff_t b[], int *pd_b, ff_t a[], int d_a, int n)
{
    register int i;
    
    for ( i = d_a+n ; i >= n ; i-- ) _ff_set(b[i],a[i-n]);
    while ( i >= 0 ) { _ff_set_zero(b[i]); i--; }
    *pd_b = d_a+n;
}

// overlap ok
static inline void ff_poly_shift_down (ff_t b[], int *pd_b, ff_t a[], int d_a, int n)
{
    register int i;
    
    *pd_b = d_a-n;
    for ( i = 0 ; i <= *pd_b ; i++ ) _ff_set(b[i],a[i+n]);
}

// overlap ok -- note that if b is already monic and c=b then nothing will be done
static inline ff_t *ff_poly_monic (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
    register ff_t z;
    register int i;
    
    if ( d_b < 0 || _ff_one (b[d_b]) ) return ff_poly_copy (c, pd_c, b, d_b);
    _ff_invert (z, b[d_b]);
    if ( pd_c ) *pd_c = d_b;
    for ( i = 0 ; i < d_b ; i++ ) ff_mult(c[i],z,b[i]);
    _ff_set_one(c[i]);
    return c;
}

// discriminant and root-finding code specialized for polys of degree d <=4
    
static inline int ff_poly_x3axb_disc (ff_t d[1], ff_t f[2])                 // f(x)=x^3+ax+b, returns D =-4a^3-27b^2
{
    register ff_t t0, t1, t2;
    
    _ff_square(t0,f[0]);  _ff_set_ui(t1,27);  ff_mult(t2,t0,t1);                // t2 = 27f0^2
    _ff_square(t0,f[1]);  _ff_mult(t1,t0,f[1]);  _ff_x2(t1); _ff_x2(t1);        // t1 = 4f1^3
    _ff_addto(t1,t2); _ff_neg(d[0],t1);
    return ! _ff_zero(d[0]);
}
   
// assumes f monic, char > 3 (not checked), depresses in place and sets t to translation value
static inline void ff_poly_depress_cubic (ff_t t[1], ff_t f[3])
{
    register ff_t t1, t2, t3, t4;
    
    _ff_mult(t1,_ff_third,f[2]);        // f2/3
    _ff_set(t[0],t1);
    _ff_square(t2,t1);                  // f2^2/9
    _ff_add(t3,t2,t2);                  // 2f2^2/9
    _ff_sub(t4,t3,f[1]);
    _ff_addto(t2,t3);
    _ff_subfrom(f[1],t2);
    _ff_mult(t2,t1,t4);
    _ff_addto(f[0],t2);
    _ff_set_zero(f[2]);
    // 3M+5A
}

// assumes f monic (does not check) depresses in place and sets t to translation value
static inline void ff_poly_depress_quartic (ff_t t[1], ff_t f[4])
{
    register ff_t t0, t1, t2, t3, t4;

    _ff_mult(t1,_ff_fourth,f[3]);
    _ff_set(t[0],t1);
    _ff_neg(t0,t1);                         // t0 = c
    _ff_square(t2,t0);                      // t2 = c^2
    _ff_triple(t3,t2);                      // t3 = 3c^2
    _ff_sub (t4,f[2],t3);                   // t4 = f2-3c^2
    _ff_sum_2_mults(t1,t0,t2,t4,f[1]);      // t1 = cf1+c^2(f2-3c^2)
    _ff_addto(f[0],t1);                     // f0 += cf1+c^2(f2-3c^2)
    _ff_subfrom(t4,t2);                     // t4 = f2-4c^2
    _ff_mult(t1,t0,t4);                     // t1 = c(f2-4c^2)
    _ff_add(t0,t1,t1);                      // t0 = 2c(f2-4c^2)
    _ff_addto(f[1],t0);                     // f1 += 2c(f2-4c^2)
    _ff_add(t0,t3,t3);                      // t0 = 6c^2
    _ff_subfrom(f[2],t0);                   // f2 -= 6c^2
    _ff_set_zero(f[3]);                     // f3 = 0
    // 5M+10A
}

// computes cubic resolvent of f=x^4+f2x^2+f1x+f0 and depresses it
static inline void ff_poly_depressed_cubic_resolvent (ff_t t[1], ff_t g[3], ff_t f[3])
{
    register ff_t t1,t2;

    _ff_set_one(g[3]);
    _ff_add(t1,f[2],f[2]);
    _ff_neg(g[2],t1);                                                           // g2 = -2f2
    _ff_add(t1,f[0],f[0]); _ff_x2(t1);
    _ff_square(t2,f[2]);
    _ff_sub(g[1],t2,t1);                                                        // g1 = f2^2-4f0
    _ff_square(g[0],f[1]);                                                      // g0 = f1^2
    ff_poly_depress_cubic (t,g);
    // total 5M+9A
}


static inline int ff_poly_roots_d1 (ff_t r[1], ff_t f[2])               // f needn't be monic
{
    register ff_t t1;
    
    if ( _ff_one(f[1]) ) { _ff_set_one(t1); } else { _ff_invert(t1,f[1]); }
    ff_negate(t1); _ff_mult(r[0],t1,f[0]);  return 1;
}

// Compute roots of poly of degree 2, or less.  We don't assume leading coefficient is non-zero (or that f is monic))
static inline int ff_poly_roots_d2 (ff_t r[2], ff_t f[3], int d_f)
{
    ff_t t1,t2,D;
    
    switch (d_f) {
    case 0: return 0;
    case 1: 
        if ( _ff_one(f[1]) ) { _ff_set_one(t1); } else { if ( _ff_zero(f[1]) ) return 0; _ff_invert(t1,f[1]); }
        ff_negate(t1); _ff_mult(r[0],t1,f[0]);  return 1;
    case 2:
        _ff_square(t1,f[1]);
        _ff_mult(t2,f[2],f[0]);
        _ff_x2(t2);  _ff_x2(t2);
        _ff_subfrom(t1,t2);
        if ( ! ff_sqrt(&D,&t1) ) return 0;
        if ( _ff_one(f[2]) ) {                  
            _ff_set(t1,_ff_half);
        } else {
            if ( _ff_zero(f[2]) ) return ff_poly_roots_d1 (r,f);       // put zero-test here to optimize for the monic case
            _ff_add(t2,f[2],f[2]);
            _ff_invert(t1,t2);
        }
        _ff_sub(t2,D,f[1]);
        _ff_mult(r[0],t1,t2);
        ff_negate(D);
        _ff_sub(t2,D,f[1]);
        _ff_mult(r[1],t1,t2);
        return 2;
    default:
        printf ("Invalid degree %d in ff_poly_roots_d2\n", d_f); abort();
    }
}

static inline int ff_poly_roots_d3 (ff_t r[3], ff_t f[4])           // f must be monic, returns all roots with multiplicity (not nescarilly distinct)
{
    ff_t g[4],t;
    int i, k;

    assert ( _ff_one(f[3]) );
    if ( _ff_zero(f[2]) ) return _ff_poly_roots_d3 (r, f, 0, 0);
    if ( _ff_p == 3 ) return _ff_poly_roots_d3_mod3 (r, f);
    for ( i = 0 ; i <= 3 ; i++ ) _ff_set(g[i],f[i]);
    ff_poly_depress_cubic(&t,g);
    k = _ff_poly_roots_d3 (r, g, 0, 0);
    for ( i = 0 ; i < k ; i++ ) _ff_subfrom(r[i],t); 
    return k;
}

static inline int ff_poly_roots_d4 (ff_t r[4], ff_t f[5])           // f must be monic, returns all roots with multiplicity (not nescarilly distinct)
{
    ff_t g[5], t;
    int i, k;

    assert ( _ff_one(f[4]) );
    if ( _ff_zero(f[3]) ) return _ff_poly_roots_d4 (r, f, 0, 0);
    for ( i = 0 ; i <= 4 ; i++ ) _ff_set(g[i],f[i]);
    ff_poly_depress_quartic(&t,g);
    k = _ff_poly_roots_d4 (r, g, 0, 0);
    for ( i = 0 ; i < k ; i++ ) _ff_subfrom(r[i],t);
    return k;
}


// f must be monic
static inline void ff_poly_remove_root_d2 (ff_t g[], ff_t f[], ff_t r[1])
{
    _ff_add(g[0],f[1],*r);
    _ff_set_one(g[1]);
}

// computes g=f/(x-r), f is monic degree 3, f and g may point to the same place
static inline void ff_poly_remove_root_d3 (ff_t g[], ff_t f[], ff_t r[1])
{
    register ff_t t1,t2;
    
    _ff_add(t1,f[2],*r);
    ff_mult(t2,t1,*r);
    _ff_add (g[0],t2,f[1]);
    _ff_set(g[1],t1);
    _ff_set_one(g[2]);
}

// computes g=f/(x-r), f is monic degree 4, f and g may point to the same place
static inline void ff_poly_remove_root_d4 (ff_t g[], ff_t f[], ff_t r[1])
{
    register ff_t r0,t1,t2;

    _ff_set(r0,r[0]);
    _ff_add(t1,f[3],r0);
    _ff_set(t2,f[2]);
    _ff_set(g[2],t1);
    ff_mult(t1,t1,r0);
    _ff_addto (t1,t2);
    _ff_set(t2,f[1]);
    _ff_set(g[1],t1);
    ff_mult(t1,t1,r0);
    _ff_add (g[0],t1,t2);
    _ff_set_one(g[3]);
}

// computes g=f/(x-r), assuming f(r) = 0.  f and g may point to the same place
static inline void ff_poly_remove_root (ff_t g[], ff_t f[], int d, ff_t r[1])
{
    register ff_t r0,t1,t2;
    register int i;

    _ff_set(r0,r[0]);
    _ff_set(t1,f[d]);
    for ( i = d-1 ; i > 0 ; i-- ) {
        _ff_set(t2,f[i]);
        _ff_set(g[i],t1);
        ff_mult(t1,t1,r0);
        _ff_addto (t1,t2);
    }
    _ff_set(g[0],t1);
}

// computes g=f*(x-r).  f is assumed to be monic and of degree > 0.  f and g may point to the same place
static inline void ff_poly_add_root (ff_t g[], ff_t f[], int d, ff_t r[1])
{
    register ff_t r0,t0;
    register int i;

    _ff_neg(r0,r[0]);
    _ff_set(g[d+1],f[d]);
    for ( i = d ; i > 0 ; i-- ) {
        _ff_mult(t0,r0,f[i]);
        _ff_add(g[i],f[i-1],t0);
    }
    _ff_mult(g[0],f[0],r0);
}

// compute o=f(x)
static inline void ff_poly_eval (ff_t o[1], ff_t f[], int d, ff_t x[1])
{
    register ff_t t,y;
    register int i;
    
    if ( d < 0 ) { _ff_set_zero(o[0]); return; }
    _ff_set (y,f[d]);
    _ff_set (t,*x);
    for ( i = d-1 ; i >= 0 ; i-- ) ff_multadd(y,y,t,f[i]);
    _ff_set(o[0],y);
}

// compute g(x)=f(a*x)
static inline void ff_poly_scale (ff_t g[], ff_t f[], int d, ff_t a)
{
    register int i;
    register ff_t x;
    
    if ( d < 0 ) return;
    _ff_set(g[0],f[0]);
    if ( !d ) return;
    _ff_set(x,a); _ff_mult(g[1],f[1],x);
    for ( i = 2 ; i <= d ; i++ ) { _ff_multby(x,a); _ff_mult(g[i],f[i],x); }
}

// naive algorithm, uses (d+1)(d+2)/2 mults and d(d+1) adds, overlap not allowed
static inline void ff_poly_square_naive (ff_t o[], ff_t f[], int d)
{
    register ff_t t;
    register int i,j;
    
    for ( i = 0 ; i <= d ; i++ ) { _ff_square(o[2*i],f[i]); if ( i ) _ff_set_zero(o[2*i-1]); }
    for ( i = 0 ; i <= d ; i++ ) for ( j = 0 ; j < i ; j++ ) { _ff_mult(t,f[i],f[j]); _ff_x2(t);  _ff_addto(o[i+j],t); }    
}

// computes c[i] = sum_{j+k=d_b+i} a[j]*b[k] for i in [0,d_a-d_b]
static inline void ff_poly_middle_product_naive (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
{
    register ff_t t;
    register int i,j;
    
    for ( i = 0 ; i <= d_a-d_b ; i++ ) { _ff_set_zero(c[i]); for ( j = 0 ; j <= d_b ; j++ ) { _ff_mult(t,a[d_b+i-j],b[j]); _ff_addto(c[i],t); } }
}

// Given an r-by-r matrix F of polys of degree <= d returns F(x)
static inline void ff_poly_matrix_eval_naive (ff_t Y[], ff_t F[], int r, int d, ff_t x[1])
{
    register int i, j;
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < r ; j++ ) ff_poly_eval(Y+i*r+j,F+(i*r+j)*(d+1),d,x);
}

// Given an r-by-r matrix F of polys of degree <= d returns F(x)
static inline void ff_poly_matrix_eval (ff_t Y[], ff_t F[], int r, int d, ff_t x[1])
{
    ff_t f[d+1];
    register int i, j;
    
    if ( r == 1 ) { ff_poly_eval(Y,F,d,x); return; }
    if ( d <= 0 ) { if ( d ) ff_matrix_set_zero(Y,r); else ff_matrix_set(Y,F,r); return; }
    _ff_set_one(f[0]); _ff_set(f[1],x[0]);
    for ( i = 2 ; i <= d; i++ ) _ff_mult(f[i],f[i-1],x[0]);
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < r ; j++ ) ff_dot_product(Y+i*r+j,F+(i*r+j)*(d+1),f,d+1);
}

// Compute H(x)=F(x)G(x) where F and G are r-by-r matrices of polys of degree <= d
static inline void ff_poly_matrix_mult (ff_t H[], ff_t F[], ff_t G[], int r, int d)
{
    ff_t f[2*d+1];
    register int i, j, k, m;
    
    for ( i = 0 ; i < r ; i++ ) {
        for ( j = 0 ; j < r ; j++ ) {
            memset(H+(i*r+j)*(2*d+1),0,(2*d+1)*sizeof(*H));
            for ( k = 0 ; k < r ; k++ ) {
                ff_poly_mult(f,0,F+(i*r+k)*(d+1),d,G+(k*r+j)*(d+1),d);
                for ( m = 0 ; m <= 2*d ; m++ ) _ff_addto(H[(i*r+j)*(2*d+1)+m],f[m]);
            }
        }
    }
}

// Compute G(x)=F(x+a) where F is an r-by-r matrix of polys of degree <= d
static inline void ff_poly_matrix_translate (ff_t G[], ff_t F[], int r, int d, ff_t a)
{
    register int i,j;
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < r ; j++ ) ff_poly_translate(G+(i*r+j)*(d+1),F+(i*r+j)*(d+1),d,a);
}

// Compute G(x)=F(ax) where F is an r-by-r matrix of polys of degree <= d
static inline void ff_poly_matrix_scale (ff_t G[], ff_t F[], int r, int d, ff_t a)
{
    register int i,j;
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < r ; j++ ) ff_poly_scale(G+(i*r+j)*(d+1),F+(i*r+j)*(d+1),d,a);
}


void ff_poly_from_roots_small (ff_t f[], ff_t r[], int d);      // d cannot exceed 64

// naive implementation: d(d-1)/2 mults. Both ff_poly_from_roots_small and ff_poly_from_roots_big are much faster
static inline void ff_poly_from_roots_naive (ff_t f[], ff_t r[], int d)
{
    register ff_t t;
    register int i,j;

    _ff_set_one(f[d]);
    _ff_neg(f[d-1],r[d-1]);
    for ( i = d-2 ; i >= 0 ; i-- ) {
        _ff_mult(f[i],r[i],f[i+1]);
        ff_negate(f[i]);
        for ( j = i+1 ; j < d-1 ; j++ ) {
            _ff_mult(t,r[i],f[j+1]);
            _ff_subfrom(f[j],t);
        }
        _ff_subfrom (f[j],r[i]);
    }
}

// g=f is ok
static inline ff_t *ff_poly_derivative (ff_t g[], int *pd_g, ff_t f[], int d_f)
{
    register int i;
    register ff_t t;
    
    if ( d_f > 0 ) {
        _ff_set(g[0],f[1]);  _ff_set_one(t);
        for ( i = 1 ; i < d_f ; i++ ) {  _ff_inc(t); _ff_mult(g[i],t,f[i+1]); }
        if ( pd_g ) *pd_g = ff_poly_degree(g,d_f-1);
    } else {
        if ( pd_g ) *pd_g = -1;
    }
    return g;
}

static inline void ff_poly_mod_xn (ff_t f[], int *pd_f, int n)
{
    if ( *pd_f < n ) return;
    *pd_f = ff_poly_degree(f,n-1);
}


// replace f(x) with f(x)-x
static inline void ff_poly_sub_x_from (ff_t f[], int *pd_f)
{
    if ( *pd_f > 1 ) { _ff_dec(f[1]); return; }
    if ( *pd_f == 1 ) { _ff_dec(f[1]); *pd_f = ( _ff_zero(f[1]) ? (_ff_zero(f[0]) ? -1 : 0 ) : 1 ); return; }
    _ff_set_neg_one (f[1]);  *pd_f = 1;
    if ( *pd_f < 0 ) _ff_set_zero(f[0]);
}


static inline int ff_poly_random (ff_t f[], int d)
    { register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_random(f[i]); return ff_poly_degree(f,d); }
    
static inline void ff_poly_random_monic (ff_t f[], int d)   // d must be >= 0
    { register int i;  for ( i = 0 ; i < d ; i++ ) _ff_random(f[i]); _ff_set_one(f[d]); return; }

static inline void ff_poly_random_depressed_monic (ff_t f[], int d) // d must be >= 1
    { register int i;  for ( i = 0 ; i < d-1 ; i++ ) _ff_random(f[i]); _ff_set_one(f[d]); _ff_set_zero(f[d-1]); return; }

#ifdef __cplusplus
}
#endif

#endif
