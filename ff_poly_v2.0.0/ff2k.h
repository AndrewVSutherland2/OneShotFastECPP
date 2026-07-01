#ifndef _FF2K_INCLUDE_
#define _FF2K_INCLUDE_

/*
    Copyright 2011-2020 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>
#include <stdio.h>
#include <wmmintrin.h>
#include <ctype.h>
#include "cstd.h"

#ifdef __cplusplus
extern "C" {
#endif

// rudimentary implementation of arithmetic in F_{2^k}
// this is not meant to be fast, just functional

#define FF2K_MAXK           31      // it would not be hard to increase this to 63

typedef unsigned long ff2k_t;       // this needs to be the same as ff_t

#define FF2K_NULL           ((ff2k_t)(1UL<<FF2K_MAXK))

struct ff2k_ctx_struct {
    ff2k_t xk, xkr, m, tr, f;
    gmp_randstate_t rs;
    int k, ri;
};

extern struct ff2k_ctx_struct _ff2k_ctx[1];

struct ff2k_ctx_struct *ff2k_setup (int k);
static inline struct ff2k_ctx_struct *ff2k_get_context (void)
    { return _ff2k_ctx; }
static inline void ff2k_set_context (struct ff2k_ctx_struct *ctx)
    { _ff2k_ctx[0] = *ctx; }

static inline int ff2k_degree(void) { return _ff2k_ctx->k; }
static inline long ff2k_card(void) { return 1L << _ff2k_ctx->k; }
static inline long ff2k_cardinality(void) { return 1L << _ff2k_ctx->k; }

// We represent FF_2^k using primitive polys taken from Hansen & Mullen Supplement in Math Comp Vol 59, No 200, 1992, pp. S47-S50
// Depending on the value of k these may or may not be Conway polynomials.  The functions below are provided for conversion.
ff2k_t _ff2k_to_conway (ff2k_t x, int k);
static inline ff2k_t ff2k_to_conway (ff2k_t x) { return _ff2k_to_conway (x, _ff2k_ctx->k); }
ff2k_t _ff2k_from_conway (ff2k_t x, int k);
static inline ff2k_t ff2k_from_conway (ff2k_t x) { return _ff2k_from_conway (x, _ff2k_ctx->k); }

// given element of FF_2^k in HM rep, returns it as an element of FF_2^(nk) in HM rep
ff2k_t _ff2k_lift (ff2k_t x, int k, int n);
static inline ff2k_t ff2k_lift (ff2k_t x, int n) { return _ff2k_lift (x, _ff2k_ctx->k, n); }

static inline ff2k_t ff2k_rand (void)
{
    if ( ! _ff2k_ctx->ri ) { gmp_randinit_default (_ff2k_ctx->rs); gmp_randseed_ui (_ff2k_ctx->rs, cstd_seed()); _ff2k_ctx->ri = 1; }
    return gmp_urandomb_ui (_ff2k_ctx->rs, _ff2k_ctx->k);
}

// use of these macros is deprecated
#define _ff2k_set(o,a)      ((o) = (a))
#define _ff2k_equal(a,b)    ((a) == (b))
#define _ff2k_add(o,a,b)    ((o) = (a)^(b))
#define _ff2k_addto(o,a)    ((o) ^= (a))
#define _ff2k_sub(o,a,b)    _ff2k_add(o,a,b)
#define _ff2k_subfrom(o,a)  _ff2k_addto(o,a)
#define _ff2k_inc(o)        ((o) ^= 1UL)
#define _ff2k_zero(a)       (!(a))
#define _ff2k_set_zero(o)   ((o)=0)
#define _ff2k_one(a)        ((a) == 1UL)
#define _ff2k_set_one(o)    ((o)=1UL)
#define _ff2k_negate(o)
#define _ff2k_next(o)       (++(o) == _ff2k_ctx->xk ? (o)=0 : (o))

// inlines consistent with fff.h interface
static inline ff2k_t ff2k_zero(void) { return 0; }
static inline ff2k_t ff2k_one(void) { return 1; }
static inline ff2k_t ff2k_is_zero (ff2k_t a) { return !a; }
static inline ff2k_t ff2k_is_one (ff2k_t a) { return a==1; }
static inline ff2k_t ff2k_add (ff2k_t a, ff2k_t b) { return a^b; }
static inline ff2k_t ff2k_next (ff2k_t a) { return (a+1 == _ff2k_ctx->xk) ? 0 : a+1; }
static inline ff2k_t ff2k_from_int (long a) { return a&1; }

static inline ff2k_t _ff2k_pclmul (ff2k_t a, ff2k_t b)
{
    register __m128i A, B, C;
    A[0]=a; B[0] = b;
    C = _mm_clmulepi64_si128 (A,B,0);
    return (ff2k_t)C[0];
}

//#define _ff2k_pclmultby(a,b) asm ("pclmullqlqdq %1, %0;" : "+x"(a) : "x"(b));   // use Intel carry-less multiplication to compute a = a*b

static inline ff2k_t _ff2k_reduce(ff2k_t a)
{
    register ff2k_t q, r;

    // use Barret reduction without precomp (our primitive polys all have the upper half coeffs clear)
    // see https://www.cosic.esat.kuleuven.be/publications/article-1115.pdf
    q = a >> _ff2k_ctx->k; q = _ff2k_pclmul(q,_ff2k_ctx->f) >> _ff2k_ctx->k;
    r = a & _ff2k_ctx->m;  q = _ff2k_pclmul(q,_ff2k_ctx->f) & _ff2k_ctx->m;
    return r^q;
}

static inline ff2k_t ff2k_mul (ff2k_t a, ff2k_t b)
    { return _ff2k_reduce(_ff2k_pclmul(a,b)); }
static inline ff2k_t ff2k_sqr (ff2k_t a) { return ff2k_mul(a,a); }             // squaring is no faster than multiplying using carry-less mults

static inline ff2k_t ff2k_multiply (ff2k_t a, ff2k_t b) { return ff2k_mul(a,b); }   // backward compatibility
static inline ff2k_t ff2k_square (ff2k_t a) { return ff2k_sqr(a); }                 // backward compatibility

static inline ff2k_t ff2k_sum_2_muls (ff2k_t a1, ff2k_t a2, ff2k_t b2, ff2k_t b1)   // returns a1*b1 + a2*b2
    { return _ff2k_reduce (_ff2k_pclmul(a1,b1)^_ff2k_pclmul(a2,b2)); }

static inline ff2k_t ff2k_inv (ff2k_t a)
{
    register ff2k_t u, v, g1, g2;
    if ( _ff2k_zero(a) ) a=1/a; // force division by zero
    // this could be made faster but we don't use it much
    u = a; v = _ff2k_ctx->f;
    g1 = 1; g2 = 0;
    while ( u!=1 && v!=1 ) {
        while ( !(u&1) ) { u>>=1;  g1 = (g1&1) ? ((g1^_ff2k_ctx->f)>>1) : g1>>1; }
        while ( !(v&1) ) { v>>=1;  g2 = (g2&1) ? ((g2^_ff2k_ctx->f)>>1) : g2>>1; }
        if ( u > v ) { u ^=v; g1 ^= g2; } else { v ^= u; g2 ^= g1; }
    }
    return (u==1) ? g1 : g2;
}
static inline ff2k_t ff2k_invert (ff2k_t a) { return ff2k_inv(a); } // backward compatibility

static inline ff2k_t ff2k_trace (ff2k_t a)
    { return __builtin_parityl(a&_ff2k_ctx->tr); }

// deprecated, for backward compatibilty only
#define _ff2k_mult(o,a,b)       ((o)=ff2k_mul(a,b))
#define _ff2k_multby(o,a)       ((o)=ff2k_mul(o,a))
#define _ff2k_square(o,a)       ((o)=ff2k_sqr(a))
#define _ff2k_invert(o,a)       ((o)=ff2k_inv(a))
#define _ff2k_trace(o,a)        ((o)=ff2k_trace(a))

// one can do this slightly more quickly, but this is fast enough
static inline ff2k_t ff2k_sqrt (ff2k_t a)
{
    register int i;
    ff2k_t x = a;
    for ( i = 0 ; i < _ff2k_ctx->k-1 ; i++ ) x = ff2k_square(x);
    return x;
}
    
static inline ff2k_t ff2k_primitive_root () { return (_ff2k_ctx->k == 1) ? 1UL : 2UL; }

static inline int ff2k_poly_degree (ff2k_t f[], int d)
    { for ( ; d >= 0 && !f[d] ; d-- ); return d; }

static inline int ff2k_poly_set_mpz (ff2k_t f[], mpz_t F[], int d)
    { for ( int i = 0 ; i <= d ; i++ ) f[i] = mpz_tstbit(F[i],0) ? 1 : 0;  return ff2k_poly_degree(f,d); }

static inline int ff2k_poly_set_i (ff2k_t f[], long F[], int d)
    { for ( register int i = 0 ; i <= d ; i++ ) f[i] = F[i] >= 0 ? ((F[i]&1) ? 1 : 0) : (((-F[i])&1) ? 1 : 0 );  return ff2k_poly_degree(f,d); }
    
// computes o = f(x)
static inline ff2k_t ff2k_poly_eval (ff2k_t f[], int d, ff2k_t x)
{
    register ff2k_t y;
    register int i;
    
    if ( d < 0 ) { return 0; }
    for ( y = f[d], i = d-1 ; i >= 0 ; i-- ) { y = ff2k_mul(y, x);  y = ff2k_add (y,f[i]); }
    return y;
}

// computes x[0]*y[n-1]+x[1]*y[n-2]+...+x[n-2]*y[1]+x[n-1]*y[0]
static inline ff2k_t ff2k_conv (ff2k_t x[], ff2k_t y[], int n)
    { ff2k_t z = 0; for ( register int i = 0 ; i <= n ; i++ ) { z = ff2k_add(z,ff2k_mul (x[i],y[n-i])); } return z; }

// computes x[0]*y[n-1]+x[1]*y[n-2]+...+x[m-1]*y[n-m+1]+x[m]*y[n-m] with 0 <= m <= n
static inline ff2k_t ff2k_aconv (ff2k_t x[], int m, ff2k_t y[], int n)
    { ff2k_t z = 0; for ( register int i = 0 ; i <= m ; i++ ) { z = ff2k_add(z,ff2k_mul (x[i],y[n-i])); } return z; }

// computes x[0]*x[n-1]+x[1]*x[n-2]+...+x[n-2]*x[1]+x[n-1]*x[0]
static inline ff2k_t ff2k_uconv (ff2k_t x[], int n)
    { ff2k_t z = 0; for ( register int i = 0 ; i <= n ; i++ ) { z = ff2k_add(z,ff2k_mul (x[i],x[n-i])); } return z; }

static inline ff2k_t *ff2k_poly_mult (ff2k_t o[], ff2k_t f[], int df, ff2k_t g[], int dg)
{
    if ( df > dg ) return ff2k_poly_mult (o, g, dg, f, df);
    for ( register int i = df ; i >=0 ; i-- ) o[dg+i] = ff2k_aconv(f+i,df-i+1,g+i,dg-i+1);
    for ( register int i = dg-1 ; i >= df ; i-- ) o[i] = ff2k_aconv(f,df,g,i+1);
    for ( register int i = df-1 ; i >= 0 ; i-- ) o[i] = ff2k_conv(f,g,i+1);
}

static inline ff2k_t *ff2k_poly_square (ff2k_t o[], ff2k_t f[], int d)
{
    for ( register int i = d ; i >=0 ; i-- ) o[d+i] = ff2k_uconv(f+i,d-i+1);
    for ( register int i = d-1 ; i >= 0 ; i-- ) o[i] = ff2k_uconv(f,i+1);
    return o;
}

// brute-force root-finding, only use if k is small
static inline int ff2k_poly_distinct_roots (ff2k_t r[], ff2k_t f[], int d)
{
    register ff2k_t x, t;
    register int i;
    
    x = ff2k_zero();
    i = 0;
    do {
        t = ff2k_poly_eval (f, d, x);
        if ( ff2k_is_zero(t) ) r[i++] = x;
        x = ff2k_next(x);
    } while ( ! ff2k_is_zero(x) );
    return i;
}

// brute-force root-finding, only use if k is small
static inline int ff2k_poly_count_distinct_roots (ff2k_t f[], int d)
{
    register ff2k_t x, t;
    register int i;
    
    x = ff2k_zero();
    i = 0;
    do {
        t = ff2k_poly_eval (f, d, x);
        if ( ff2k_is_zero(t) ) i++;
        x = ff2k_next(x);
    } while ( ! ff2k_is_zero(x) );
    return i;
}

// naive point-counting for elliptic curves in WS form over F_2^k, returns -1 if curve is singular
long ff2k_WS_pointcount (ff2k_t w[]);

// naive point-counting for hyperelliptic curves over F_2^k, assumes curve is not singular but does not check
long ff2k_hyperelliptic_pointcount (ff2k_t f[], int df, ff2k_t h[], int dh);
long ff2k_hyperelliptic_pointcount_at_infinity (ff2k_t f[], int df, ff2k_t h[], int dh); // returns # of points at infinity (assuming good reduction), 0, 1, 2
static inline char *ff2k_sprint (char buf[], ff2k_t x, int conway)
{
    register char *s=buf;

    if ( conway ) x=ff2k_to_conway(x);
    *s++ = '['; *s++ = '0'+(x&1);
    for ( x|=(((ff2k_t)1)<<_ff2k_ctx->k), x >>= 1 ; x!=1 ; x >>= 1 ) { *s++ = ','; *s++ = '0'+(x&1); };
    *s++ = ']'; *s++ = '\0';
    return buf;
}
        
static inline void ff2k_print (ff2k_t x, int conway)
    { char buf[2*FF2K_MAXK+3]; printf("%s",ff2k_sprint(buf,x,conway)); }


// The root counting code below is based on Zinoviev's "On the solutions of equations of degree" RR-2829, 1996, inria-00073862
// https://hal.inria.fr/inria-00073862
// returns number of disinct roots of x^2+bx+c
static inline int ff2k_poly_count_distinct_roots_monic_d2 (ff2k_t b, ff2k_t c)
    {   return ff2k_is_zero(b) ? 1 : (ff2k_is_zero(ff2k_trace(ff2k_mul(c,ff2k_inv(ff2k_sqr(b))))) ? 2 : 0); }

static inline int ff2k_poly_count_distinct_roots_d2 (ff2k_t f[3])
{
    if ( ! ff2k_is_one(f[2]) ) {
        register ff2k_t t = ff2k_invert(f[2]);
        return ff2k_poly_count_distinct_roots_monic_d2 (ff2k_mul(t,f[1]),ff2k_mul(t,f[0]));
    } else {
        return ff2k_poly_count_distinct_roots_monic_d2 (f[1],f[0]);
    }
}

// accepts binary strings 1010 or 0b1010 (terminated by whitespace or punctuation), binary arrays [0,1,0,1] (least significant bit first, no whitespace)
// sets z to FF2K_NULL if unable to parse, returns character where parsing stopped
static inline char *ff2k_parse (ff2k_t *z, char *s, int conway)
{
    register ff2k_t x;
    register int i;

    *z = FF2K_NULL;
    if ( ! *s ) return 0;
    if ( ! *(s+1) ) { if ( *s >='0' && *s <= '1' ) { *s = *s-'0'; }; return s+1; }
    if ( s[0] == '[' ) {
        for ( x = 0, i = 0, s++ ; *s >= '0' && *s <= '1' && i < _ff2k_ctx->k ; i++, s++ ) {
            if (*s=='1') x |= ((ff2k_t)1)<<i;
            if (*(++s) != ',') break;
        }
        if ( *s != ']' ) return s;
        s++;
    } else {
        if ( s[0] == 0 && s[1] == 'b' ) s+=2;
        for ( x = 0, i = 0 ; i < _ff2k_ctx->k && *s >= '0' && *s <= '1' ; s++, i++ ) x = (x<<1) | (*s-'0');
        if ( *s && ! isspace(*s) && ! ispunct(*s) ) return s;
    }
    *z = conway ? ff2k_from_conway(x) : x;
    return s;
}

static inline char *ff2k_poly_sprint (char *buf, ff2k_t f[], int d, int conway)
{
    register char *s = buf;
    *s++ = '[';
    for ( register int i = 0 ; i <= d ; i++ ) { if ( i ) *s++ = ',';  s += strlen(ff2k_sprint(s,f[i],conway)); }
    *s++ = ']';  *s++ = '\0';
    return buf;
}

static inline void ff2k_poly_print (ff2k_t f[], int d, int conway)
    { char buf[(d+1)*(2*FF2K_MAXK+3)+32]; printf("%s",ff2k_poly_sprint(buf,f,d,conway)); }

// accepts [f0,f1,...,fd] where each fi is parseable by ff2k_parse
// sets def to -1 for zero poly, -2 if error, returns pointer to end of parse
static inline char *ff2k_poly_parse (ff2k_t f[], int *df, int maxd, char *s, int conway)
{
    register int d;

    while ( isspace (*s) ) s++;
    *df = -2;
    if ( *s != '[' ) return s;
    s++; while ( isspace (*s) ) s++;
    for ( d = -1 ; d <= maxd && *s && *s != ']' ; s++ ) {
        s = ff2k_parse (f + ++d, s, conway);
        if ( f[d] == FF2K_NULL ) return s;
        while ( isspace (*s) ) s++;
        if ( *s != ',' ) break;
    }
    if ( *s != ']' ) return s;
    while ( d >= 0 && !f[d] ) d--;
    *df = d;
    return ++s;
}

// accepts [f,h] defining a curve y^2 + h(x)y = f(x), where f,h are parseable by ff2k_poly_parse
// returns (an upper bound on) the genus (-1 for error)
static inline int ff2k_hyperelliptic_parse (ff2k_t f[], int *df, ff2k_t h[], int *dh, int maxg, char *s, int conway)
{
    register int g;

    while ( isspace (*s) ) s++;
    if ( *s != '[' ) return -1;
    s = ff2k_poly_parse (f, df, 2*maxg+2, ++s, conway);
    if ( *df < -1 ) return -1;
    while ( isspace (*s) ) s++;
    if ( *s != ',' ) return -1;
    s = ff2k_poly_parse (h, dh, maxg+1, ++s, conway);
    if ( *dh < -1 ) return -1;
    while ( isspace (*s) ) s++;
    if ( *s != ']' ) return -1;
    g = (*df-1)/2;
    if ( *dh-1 > g ) g = *dh-1;
    return g;
}

static inline char *ff2k_hyperelliptic_sprint (char *buf, ff2k_t f[], int df, ff2k_t h[], int dh, int conway)
{
    char fbuf[(df+1)*(2*FF2K_MAXK+3)+32], hbuf[(dh+1)*(2*FF2K_MAXK+3)+32];
    sprintf (buf,"[%s,%s]", ff2k_poly_sprint (fbuf, f, df, conway), ff2k_poly_sprint (hbuf, h, dh, conway));
    return buf;
}

static inline void ff2k_hyperelliptic_print (ff2k_t f[], int df, ff2k_t h[], int dh, int conway)
    { char buf[(df+dh+4)*(2*FF2K_MAXK+3)+32]; printf("%s",ff2k_hyperelliptic_sprint(buf,f,df,h,dh,conway)); }

int ff2k_poly_count_distinct_roots_d2 (ff2k_t f[3]);
int ff2k_poly_count_distinct_roots_d3 (ff2k_t f[4]);
int ff2k_poly_count_distinct_roots_d4 (ff2k_t f[5]);

void ff2k_Lmatrix_slow (unsigned M[], ff2k_t a, ff2k_t b);    
void ff2k_Lmatrix_fast (unsigned M[], ff2k_t a, ff2k_t b);    
static inline void ff2k_Lmatrix (unsigned M[], ff2k_t a, ff2k_t b)
    { (_ff2k_ctx->k <= 4) ? ff2k_Lmatrix_slow(M,a,b) : ff2k_Lmatrix_fast(M,a,b); }

void ff2_transpose (unsigned T[], unsigned M[], int rows, int cols);
int ff2_row_reduce (unsigned M[], int rows, int cols);
int ff2_row_reduce_mat (unsigned R[], unsigned M[], int rows, int cols);
static inline int ff2k_Rmatrix (unsigned R[], ff2k_t a, ff2k_t b)
{
    unsigned M[FF2K_MAXK];
    ff2k_Lmatrix(R, a, b);
    ff2_transpose (M, R, _ff2k_ctx->k, _ff2k_ctx->k);
    return ff2_row_reduce_mat (R, M, _ff2k_ctx->k, _ff2k_ctx->k);
}

int ff2k_poly_count_distinct_roots_d3_linearize (ff2k_t L[5], ff2k_t f[4]);
int ff2k_poly_count_distinct_roots_d4_linearize (ff2k_t L[5], ff2k_t f[5]);


#ifdef __cplusplus
}
#endif

#endif
