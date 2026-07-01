#ifndef _POLYPARSE_INCLUDE_
#define _POLYPARSE_INCLUDE_

/*
    Copyright 2011-2018 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
    The generic parsing functions below support multiple data types via the helper functions setzero, addto, and iszero, which operate on  a single coefficient f[i].
    f points to an array of maxd+1 coefficients, which may be of any datatype.  Note that addto is passed a rational number, and may return 0 if it can't be converted to the desired datatype.

    Note that the addto function is allowed to trash c if it likes (so it can use it as workspace)

    This isn't superfast, but it saves repeating the same code many times, and parsing poly strings is generally not a limiting performance issue.
    
    Returns the degree of the poly, with degree -1 for the zero polynomial.
    Returns -2 if the poly degree exceeds maxd, and returns -3 if an error is encountered.
*/
// poly_parse expects a univariate polynomial
int poly_parse_array (void *f, int off, int max, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);
int poly_parse (void *f, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);

// nf_poly_parse expects expr to be of the form [(coeff-poly in a)*x^d + ... + (coeff-poly in a)]/(minpoly of a)
int nf_poly_parse (void *f, int *df, void *g, int *dg, int maxd, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);

int poly_parse_hyperelliptic_curve (void *f, int *df, void *h, int *dh, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);
int nf_poly_parse_hyperelliptic_curve (void *f, int *df, void *h, int *dh, void *g, int *dg, int maxd, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);

int poly_parse_superelliptic_curve (int *m, void *f, int *df, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);

// poly_parse_plane_curve expects a plane curve of the specified degree in affine or project coords, f should point to an array of binom(d+2,2) coefficients stored in monomial lex order (e.g. x^2, x*y, x*z, y^2, y*z, z^2)
// the indeterminates do not have to be x,y,z, but they must be alpha chars and will be sorted alphabetically (if only one appears, it will be  treated as y)
int poly_parse_plane_curve (void *f, int d ,char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);
int nf_poly_parse_plane_curve (void *f, int deg, void *g, int *dg, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);

static inline int poly_parse_plane_conic (void *f, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
    { return poly_parse_plane_curve (f, 2, expr, setzero, addto, iszero, arg); }
static inline int poly_parse_plane_cubic (void *f, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
    { return poly_parse_plane_curve (f, 3, expr, setzero, addto, iszero, arg); }
static inline int poly_parse_plane_quartic (void *f, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
    { return poly_parse_plane_curve (f, 4, expr, setzero, addto, iszero, arg); }

static inline int nf_poly_parse_plane_conic (void *f, void *g, int *dg, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
    { return nf_poly_parse_plane_curve (f, 2, g, dg, maxn, expr, setzero, addto, iszero, arg); }
static inline int nf_poly_parse_plane_cubic (void *f, void *g, int *dg, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
    { return nf_poly_parse_plane_curve (f, 3, g, dg, maxn, expr, setzero, addto, iszero, arg); }
static inline int nf_poly_parse_plane_quartic (void *f, void *g, int *dg, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
    { return nf_poly_parse_plane_curve (f, 4, g, dg, maxn, expr, setzero, addto, iszero, arg); }

int poly_parse_conic_cover (void *g, void *f, int d ,char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);
    

// returns the degree of the highest degree term listed the poly specified by expr (actual degree could be less due to cancellation).
// expr must either be a univariate poly or homogeneous of degree 3
int poly_parse_max_degree (char *expr);
// sets vars to polynomial variable names (single characters) appearing in polynomial expression, returns number of variables
int poly_parse_variables (char vars[3], char *expr);
// poly_parse_minpoly returns a pointer to the opening paren of the minpoly in [...]/(...) or null if poly-expr is not over a number field
char *poly_parse_minpoly (char *expr);
char poly_parse_minpoly_variable (char *expr);
static inline int poly_parse_minpoly_degree (char *expr) { char *s = poly_parse_minpoly (expr); if ( ! s ) return -1; return poly_parse_max_degree (s); }
    
// parses polynomial expression with rational coefficients into an array of unsigned longs reduced mod p -- returns 0 if any coefficient contains an univertible denominator -- implemented via poly_parse
int ui_poly_parse_mod_p (unsigned long f[], int maxd, char *expr, unsigned long p);
int ui_poly_parse_plane_curve_mod_p (unsigned long f[], int d, char *expr, unsigned long p);
static inline int ui_poly_parse_plane_conic_mod_p (unsigned long f[6], char *expr, unsigned long p)
    { return ui_poly_parse_plane_curve_mod_p (f, 2, expr, p); }
static inline int ui_poly_parse_plane_cubic_mod_p (unsigned long f[10], char *expr, unsigned long p)
    { return ui_poly_parse_plane_curve_mod_p (f, 3, expr, p); }
static inline int ui_poly_parse_plane_quartic_mod_p (unsigned long f[15], char *expr, unsigned long p)
    { return ui_poly_parse_plane_curve_mod_p (f, 4, expr, p); }

#define MAXPLANEDEG                 9

// parses polynomial expression with integer coefficients into an array of longs -- returns 0 if any coefficient is non-integral, or in the case of overflow
int i_poly_parse (long f[], int maxd, char *expr);
int i_nf_poly_parse (long f[], int *df, long g[], int *dg, int maxd, int maxn, char *expr) ;

// these functions return the genus of the curve (assuming it is not singular)
int i_poly_parse_hyperelliptic_curve (long f[], int *df, long h[], int *dh, int maxd, char *expr);
int i_nf_poly_parse_hyperelliptic_curve (long f[], int *df, long h[], int *dh, long g[], int *dg, int maxd, int maxn, char *expr);
int i_poly_parse_superelliptic_curve (int *m, long f[], int *df, int maxd, char *expr);

// coefficient orderding for plane curves is lex monomial x...x, x...xy, ..., z...z (reverse lex order on exponent vectors)
int i_poly_parse_plane_curve (long f[], int d, char *expr);
static inline int i_poly_parse_plane_conic (long f[], char *expr)
    { return i_poly_parse_plane_curve (f, 2, expr); }
static inline int i_poly_parse_plane_cubic (long f[], char *expr)
    { return i_poly_parse_plane_curve (f, 3, expr); }
static inline int i_poly_parse_plane_quartic (long f[], char *expr)
    { return i_poly_parse_plane_curve (f, 4, expr); }
void i_poly_print_plane_curve (long f[], int d);
int i_poly_sprint_plane_curve (char *s, long f[], int d);
static inline void i_poly_print_plane_conic (long f[6])
    { i_poly_print_plane_curve (f, 2); }
static inline int i_poly_sprint_plane_conic (char *s, long f[6])
    { return i_poly_sprint_plane_curve (s, f, 2); }
static inline void i_poly_print_plane_cubic (long f[10])
    { i_poly_print_plane_curve (f, 3); }
static inline int i_poly_sprint_plane_cubic (char *s, long f[10])
    { return i_poly_sprint_plane_curve (s, f, 3); }
static inline void i_poly_print_plane_quartic (long f[15])
    { i_poly_print_plane_curve (f, 4); }
static inline int i_poly_sprint_plane_quartic (char *s, long f[15])
    { return i_poly_sprint_plane_curve (s, f, 4); }

int i_poly_parse_conic_cover (long g[6], long f[], int d, char *expr);
void i_poly_print_conic_cover (long g[6], long f[], int d);
int i_poly_sprint_conic_cover (char *s, long g[6], long f[], int d);

int ui_poly_sprint (char *s, unsigned long f[], int df);
void ui_poly_print (unsigned long f[], int df);
int i_poly_sprint (char *s, long f[], int df);
void i_poly_print (long f[], int df);
int i_poly_sprint_coeffs (char *s, long f[], int df);
void i_poly_print_coeffs (long f[], int df);

// no overflow check
static inline long i_poly_eval (long f[], int df, long x)
{
    register long y, *c;
    
    if ( df < 0 ) return 0;
    c = f+df; y = *c;
    while ( c > f ) y = x*y+*(--c);
    return y;
}

// Provided x is reduced mod p < 2^31 and coeffs of f are bounded by LONG_MAX/2, overflow should not be a problem
static inline long i_poly_eval_mod_p (long f[], int df, long x, long p)
{
    register long y, *c;
    
    if ( df < 0 ) return 0;
    c = f+df; y = *c % p;
    while ( c > f ) y = (x*y+*(--c))%p;
    if ( y < 0 ) y += p;
    return y;
}

// parses hyperelliptic curve spec in the form [f] or [f,h] defining the curve y^2 + h(x)y = f(x).  Returns the genus, or 0 for failure
// will allow elliptic curves as well, may modify (and then unmodify) str
static inline int i_hyperelliptic_curve_parse (long f[], int *df, long h[], int *dh, int maxd, char *str)   // maxd bounds the degree of f, (maxd+1)/2 bounds the degree of h, returns 0/1 for failure/success
    { return i_poly_parse_hyperelliptic_curve (f, df, h, dh, maxd, str); }

static inline int i_superelliptic_curve_parse (int *m, long f[], int *df, int maxd, char *str)              // maxd bounds the degree of f, (maxd+1)/2 bounds the degree of h, returns 0/1 for failure/success
    { return i_poly_parse_superelliptic_curve (m, f, df, maxd, str); }

static inline int i_hyperelliptic_curve_normalize (long g[], long f[], int df, long h[], int dh)            // sets g to 4*f+h^2, or just to f if h is zero (g is allowed to alias f but not h)
{
    register int i, j;
    
    if ( dh < 0 ) { if ( g != f ) for ( i = 0 ; i <= df ; i++ ) g[i] = f[i];  return df; }
    for ( i = 0 ; i <= df ; i++ ) g[i] = 4*f[i];
    while ( i <= 2*dh ) g[i++] = 0;
    for ( i = 0 ; i <= dh ; i++ ) for ( j = 0 ; j <= dh ; j++ ) g[i+j] += h[i]*h[j];
    if ( 2*dh > df ) df = 2*dh;
    for ( i = 0 ; i <= df ; i++ ) if ( (f[i]&3) ) break;
    if ( i > df ) for ( i = 0 ; i <= df ; i++ ) f[i] /= 4;
    return df;
}

static inline int i_hyperelliptic_curve_sprint (char *s, long f[], int df, long h[], int dh)
{
    char *t = s + i_poly_sprint (s, f, df) ;
    char *u = t + i_poly_sprint (t, h, dh);
    *(t-1) = ','; *t = ' ';
    return u-s;
}
void i_hyperelliptic_curve_print (long f[], int df, long h[], int dh);

static inline int i_hyperelliptic_curve_sprint_coeffs (char *s, long f[], int df, long h[], int dh)
{
    *s = '[';
    char *t = s + i_poly_sprint_coeffs (s+1, f, df);
    *t = ',';
    char *u = t + i_poly_sprint_coeffs (t+1, h, dh);
    *u++ = ']'; *u = '\0';
    return u-s;
}
void i_hyperelliptic_curve_print_coeffs (long f[], int df, long h[], int dh);

static inline int i_superelliptic_curve_sprint (char *s, int m, long f[], int df)
{
    char *t = s + sprintf (s, "y^%d = ", m); 
    t += i_poly_sprint(t, f, df);
    return t-s;
}
void i_supperelliptic_curve_print (int m, long f[], int df);

#ifdef __cplusplus
}
#endif

#endif
