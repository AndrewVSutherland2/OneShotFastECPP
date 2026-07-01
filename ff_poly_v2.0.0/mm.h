#ifndef _MONTGOMERY64_INCLUDE_
#define _MONTGOMERY64_INCLUDE_

#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

/*
    Copyright 2017-2019 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#ifdef __cplusplus
extern "C" {
#endif

#define SWAP(a,b)           do { typeof(a) _SWAP_TMP_ = a; a = b; b = _SWAP_TMP_; } while (0)

#ifndef MM_BITS
#define MM_BITS 64
#endif

#if MM_BITS == 64
typedef uint64_t mm_t;
typedef int64_t mms_t;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
typedef unsigned __int128 mm2_t; //typedef unsigned int mm2_t __attribute__((mode(TI)));
#pragma GCC diagnostic pop
#elif MM_BITS == 32
typedef uint32_t mm_t;
typedef int32_t mms_t;
typedef uint64_t mm2_t;
#else
#error MM_BITS must be 32 or 64
#endif

// all 64-bit inputs are assumed to be in [0,p-1], where p is an odd integer (p need not be prime) less than 2^(MM_BITS-1)


// computes -1/p mod 2^MM_BITS
static inline mm_t mm_pinv (mm_t p)
{
    mm_t mm_imod256[128] = { // mm_imod256[i] = 1/(2*i+1) mod 256
1,171, 205, 183, 57, 163, 197, 239, 241, 27, 61, 167, 41, 19, 53, 223, 225, 139, 173, 151, 25, 131, 165, 207, 209, 251, 29, 135, 9, 243, 21, 191, 193, 107, 141, 119, 249, 99, 133, 175, 177, 219,
253, 103, 233, 211, 245, 159, 161, 75, 109, 87, 217, 67, 101, 143, 145, 187, 221, 71, 201, 179, 213, 127, 129, 43, 77, 55, 185, 35, 69, 111, 113, 155, 189, 39, 169, 147, 181, 95, 97, 11, 45, 23,
153, 3, 37, 79, 81, 123, 157, 7, 137, 115, 149, 63, 65, 235, 13, 247, 121, 227, 5, 47, 49, 91, 125, 231, 105, 83, 117, 31, 33, 203, 237, 215, 89, 195, 229, 15, 17, 59, 93, 199, 73, 51, 85, 255 };
    register mm_t t;
    
    // Use Jebelean's trick for computing p^{-1} mod 2^64 (see HECHECC p. 190 remark 10.4 (ii))
    t = mm_imod256[((p&0xFFUL)>>1)];
    t = (2*t + (-(p*t*t)))&0xFFFFUL;
    t = (2*t + (-(p*t*t)))&0xFFFFFFFFUL;
#if MM_BITS == 64
    t = (2*t + (-(p*t*t)))&0xFFFFFFFFFFFFFFFFUL;
#endif
    return -t;
}

// returns R := 2^MM_BITS modulo p (R represents 1 in Montgomery rep)
static inline mm_t mm_R (mm_t p)
    { register mm_t x = ((mm_t)1<<(MM_BITS-1)) % p; x += x;  return x > p ? x-p : x; }

// computes x/R mod p.  Provided x < p*2^MM_BITS the loop will execute at most once
// we allow x > p*R, but if x is much larger than this may be very slow
static inline mm_t mm_redc (mm2_t x, mm_t p, mm_t pinv) // pinv = -1/p mod 2^64
    { register mm_t r = ((mm2_t)((mm_t)x*pinv) * p + x) >> MM_BITS;  while ( r >= p ) r -= p;  return r; }

// computes (x*y)/R mod p (if x and y are in Montgomery rep, the output is x*y in Montgomery rep)
static inline mm_t mm_mul (mm_t x, mm_t y, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x*y,p,pinv); }

static inline mm_t mm_sqr (mm_t x, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x*x,p,pinv); }

static inline mm_t mm_cube (mm_t x, mm_t p, mm_t pinv)
    { return mm_mul(mm_sqr(x,p,pinv),x,p,pinv); }

// computes x*y mod p, under the assumption that y*R < 2^MM_BITS (NOT VERIFIED!)
// (it will still work whenever x*y*R < 2^(2*MM_BITS), but may be very slow)
// If x is in Montgomery rep (and y is not) the returned value will be xy in Montgomery rep
static inline mm_t mm_mul_ui (mm_t x, uint64_t y, mm_t R, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x*y*R,p,pinv); }

// computes x*y*z/R mod p under the assumption that x*y*z < p*2^MM_BITS (NOT VERIFIED!)
// (it will still work whenever x*y*R < 2^(2*MM_BITS), but may be very slow)
// If x and y are in Montgomery rep (and z is not) the returned value will be xyz in Montgomery rep
static inline mm_t mm_mul_mul_ui (mm_t x, mm_t y, uint64_t z, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x*y*z,p,pinv); }

// assumes x in [0,p-1]
static inline mm_t mm_neg (mm_t x, mm_t p)
    { x = p-x;  while ( x >= p ) x -= p; return x; } // while is faster than if

// assumes x,y in [0,p-1] (if not this may loop for a *very* long time)
static inline mm_t mm_add (mm_t x, mm_t y, mm_t p)
    { x += y;  while ( x >= p ) x -= p; return x; } // while is faster than if

static inline mm_t mm_div2 (mm_t x, mm_t p)
    { return (x+p*(x&1))>>1; }
    
// computes x+y*z
static inline mm_t mm_addmul (mm_t x, mm_t y, mm_t z, mm_t R, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x*R+(mm2_t)y*z, p,pinv); }

// assumes x,y in [0,p-1]
static inline mm_t mm_sub (mm_t x, mm_t y, mm_t p)
    { mms_t u = x-y; while ( u < 0 ) u += p; return u; }
//    { x -= y; return ((((mms_t)x) >> (MM_BITS-1)) & p) + x; } // this is about the same speed

// Given R=2^MM_BITS mod p (1 in Montgomery rep), computes R2=2^(2*MM_BITS) mod p (2^MM_BITS in Montgomery rep)
static inline mm_t mm_R2 (mm_t R, mm_t p, mm_t pinv)
    { return ((mm2_t)R*R) % p; } // this is faster than powering 2 in Montgomery rep

// given R2=R^2 mod p, computes R3=R^3 mod p (R3 represents R2 in Montgomery rep)
static inline mm_t mm_R3 (mm_t R2, mm_t p, mm_t pinv)
    { return mm_mul(R2,R2,p,pinv); }
 
// computes x/R mod p (the integer represented by x in Montgomery rep)
static inline uint64_t mm_to_ui (mm_t x, mm_t p, mm_t pinv)
    { return mm_redc(x, p, pinv); }

// computes x/R mod p (the integer represented by x in Montgomery rep)
static inline int64_t mm_to_si (mm_t x, mm_t p, mm_t pinv)
    { return (int64_t)mm_redc(x, p, pinv); }

static inline mm_t mm_from_posint (int n, mm_t R, mm_t p) 
{
    mm_t x,y;
    if ( !n ) return 0;
    x = R;
    for (;!(n&1);n>>=1) x = mm_add(x,x,p);
    for (y=x,n>>=1;n;n>>=1) {
        x = mm_add(x,x,p);
        if ( (n&1) ) y = mm_add(x,y,p);
    }
    return y;
}

static inline mm_t mm_from_int (int n, mm_t R, mm_t p) 
    { return n < 0 ? mm_neg(mm_from_posint(-n, R, p),p) : mm_from_posint(n, R, p); }


// computes xR = xR^2/R mod p (the Montgomery representation of the integer x)
static inline mm_t mm_from_ui (uint64_t x, mm_t R2, mm_t p, mm_t pinv)
{
#if MM_BITS == 32
    if ( x > p ) x %= p; // we don't try to be clever here
#endif
    return mm_mul (x,R2,p,pinv);
}

// computes xR = xR^2/R mod p (the Montgomery representation of the integer x)
static inline mm_t mm_from_si (int64_t x, mm_t R2, mm_t p, mm_t pinv)
    { return x < 0 ? mm_neg(mm_from_ui(-x,R2,p,pinv),p) : mm_from_ui(x,R2,p,pinv); }

// given x=aR mod p computes R/a mod p (Montgomery inverse)
static inline mm_t mm_inv_small (mm_t r, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    register mm_t s, t, v;
    register int k;
    
    // using __builtin_ctzl is slower for p < 256, but give almost a 50% speedup when p ~ 2^60

    assert (r); k = 0;
    while ( !(r&1) ) { r >>= 1; k++; }
    t = (p-r)>>1;  v = 1; s = 2; k++;
    while ( !(t&1) ) { t >>= 1; s <<= 1; k++; }
    for (;;) {
        if ( t > r ) {
            t = (t-r)>>1; v += s; s <<= 1; k++;
            while ( !(t&1) ) { t >>= 1; s <<= 1; k++; }
        } else {
            r = (r-t)>>1; s += v; v <<= 1; k++;
            if ( ! r ) break;
            while ( !(r&1) ) { r >>= 1; v <<= 1; k++; }
        }
    }
    while ( v >= p ) v -= p;
    v = p - v; // note that v is not zero
    // we could use a lookup table for this, but the improvment is minimal (typically under 1 nanosecond)
    if ( k <= MM_BITS ) { v = mm_redc((mm2_t)v*R3,p,pinv); k += MM_BITS; } else { v = mm_redc((mm2_t)v*R2,p,pinv); }
    return mm_redc((mm2_t)v*(((mm_t)1)<<(2*MM_BITS-k)),p,pinv);
}

// given x=aR mod p computes R/a mod p (Montgomery inverse)
static inline mm_t mm_inv (mm_t r, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    register mm_t s, t, v;
    register int b,k;
    
    if ( p < 1787 ) return mm_inv_small(r,R2,R3,p,pinv);    // crossover is CPU dependent, this was tested on Intel Xeon E5/E7 v2 and v3
    
    assert (r);
    b = __builtin_ctzl(r); r >>= b;
    t = (p-r)>>1;  v = 1; s = 2; k = b+1;
    b = __builtin_ctzl(t); t >>= b; s <<= b; k+=b;
    do {
        if ( t > r ) {
            t = (t-r)>>1; v += s; s <<= 1; k++;
            b = __builtin_ctzl(t); t >>= b; s <<= b; k+=b;   // updating s and k twice is slightly faster
        } else {
            r = (r-t)>>1; s += v; v <<= 1; k++;
            if ( ! r ) break;
            b = __builtin_ctzl(r); r >>= b; v <<= b; k+=b;
        }
    } while ( r);
    while ( v >= p ) v -= p;
    v = p - v; // note that v is not zero
    // we could use a lookup table for this, but the improvment is minimal (typically under 1 nanosecond)
    if ( k <= MM_BITS ) { v = mm_redc((mm2_t)v*R3,p,pinv); k += MM_BITS; } else { v = mm_redc((mm2_t)v*R2,p,pinv); }
    return mm_redc((mm2_t)v*(((mm_t)1)<<(2*MM_BITS-k)),p,pinv);
}

// simple binary exp is faster than 2-bit sliding window
static inline mm_t mm_exp_ui (mm_t x, uint64_t e, mm_t R, mm_t p, mm_t pinv)
{
    register mm_t y;

    if (!e) return R;
    for (;!(e&1);e>>=1) x = mm_sqr(x,p,pinv);
    for (y=x,e>>=1;e;e>>=1) {
        x = mm_sqr(x,p,pinv);
        if ( (e&1) ) y = mm_mul(x,y,p,pinv);
    }
    return y;
}

// given x[i] = x^(2^i) for i in [1..floor(log_2(e))] and e in [0,..2^64-1], computes x^e
static inline mm_t mm_exp_pow (mm_t x[], uint64_t e, mm_t R, mm_t p, mm_t pinv)
{
    register mm_t y;
    register uint64_t m;
    
    switch (e) {
    case 0: return R;
    case 1: return x[0];
    case 2: return x[1];
    case 3: return mm_mul(x[0],x[1],p,pinv);
    case 4: return x[2];
    case 5: return mm_mul(x[0],x[2],p,pinv);
    case 6: return mm_mul(x[1],x[2],p,pinv);
    case 7: return mm_mul(mm_mul(x[0],x[1],p,pinv),x[2],p,pinv);
    case 8: return x[3];
    case 9: return mm_mul(x[0],x[3],p,pinv);
    case 10: return mm_mul(x[1],x[3],p,pinv);
    }
    for ( m = 1 ; !(e&m) ; m <<= 1, x++ );
    y = *x++;
    for ( m <<= 1 ; m <= e ; m <<= 1, x++ ) if ( (e&m) ) y = mm_mul(y,*x,p,pinv);
    return y;
}

static inline mm_t mm_exp_2k (mm_t x, int k, mm_t p, mm_t pinv)
    {  while ( k-- ) x = mm_sqr (x, p, pinv); return x; }

// computes the Legendre symbol (x/p)
static inline mm_t mm_legendre (mm_t x, mm_t R, mm_t p, mm_t pinv)
    { return x ? (mm_exp_ui(x,p>>1,R,p,pinv) == R ? 1 : -1) : 0; }

// computes the square root of 1/x modulo p=2^e*m+1 (with m odd), given lists of binary powers of generator a for the 2-Sylow and ainv=1/a
// returns 0 if x does not have a square-root in Fp (if ext is set then sqrt(x) = z*sqrt(a[0]) in Fp2)
int mm_invsqrt (mm_t *z, mm_t x, mm_t a[], mm_t ainv[], int e, mm_t R, mm_t p, mm_t pinv, int ext);
    
// o[i] = c*f[i] for i in [0,n-1]
static inline mm_t *mm_smul (mm_t o[], mm_t c, mm_t f[], int n, mm_t p, mm_t pinv)
    { for ( register int i = 0 ; i < n ; i++ ) o[i] = mm_mul(c,f[i],p,pinv); return o; }

// o[i] = c*f[i] for i in [0,n-1]
static inline mm_t *mm_smul_ui (mm_t o[], mm_t c, uint64_t f[], int n, mm_t R, mm_t p, mm_t pinv)
    { for ( register int i = 0 ; i < n ; i++ ) o[i] = mm_mul_ui(c,f[i],R,p,pinv); return o; }

// computes x[0]*y[0] + x[1]*y[1] + ... + x[n-1]*y[n-1] as an mm2_t (no reductions)
static inline mm2_t mm2_dot (mm_t x[], mm_t y[], int n)
    { mm2_t z = 0; for ( register int i = 0 ; i < n ; i++ ) z += (mm2_t)x[i]*y[i]; return z; }

// computes x[0]*y[n-1] + x[1]*y[n-2] + ... + x[n-1]*y[0] as an mm2_t (no reductions)
static inline mm2_t mm2_conv (mm_t x[], mm_t y[], int n)
    { mm2_t z = 0; for ( register int i = 0 ; i < n ; i++ ) z += (mm2_t)x[i]*y[n-1-i]; return z; }

// computes (x[0]*y[0] + x[1]*y[1] + ... + x[n-1]*y[n-1])/R mod p
static inline mm_t mm_dot (mm_t x[], mm_t y[], int n, mm_t p, mm_t pinv)
    { return mm_redc(mm2_dot(x,y,n),p,pinv); }

// computes x[0]*y[0] + x[1]*y[1] + ... + x[n-1]*y[n-1] mod p under the assumption that mm2_dot(x,y,n)*R < p*2^MM_BITS (NOT VERIFIED!)
// (it will still work whenever mm2_dot(x,y,n)*R < 2^(2*MM_BITS), but it may be very slow)
static inline mm_t mm_dot_ui (mm_t x[], mm_t y[], int n, mm_t R, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)R*mm2_dot(x,y,n),p,pinv); }

// computes x[0]*y[n-1] + x[1]*y[n-2] + ... + x[n-1]*y[0] (with redc)
static inline mm_t mm_conv (mm_t x[], mm_t y[], int n, mm_t p, mm_t pinv)
    { return mm_redc(mm2_conv(x,y,n),p,pinv); }

// computes x[0]*x[n-1] + x[1]*x[n-2] + ... + x[n-1]*x[0] as an mm2_t (no redc)
static inline mm2_t mm2_uconv (mm_t x[], int n)
    { mm2_t z = (n&1) ? (mm2_t)x[n/2]*x[n/2] : 0; for ( register int i = 0 ; i < n/2 ; i++ ) z += (mm2_t)x[i]*(2*x[n-1-i]); return z; }

// computes x[0]*x[n-1] + x[1]*x[n-2] + ... + x[n-1]*x[0] (with redc)
static inline mm_t mm_uconv (mm_t x[], int n, mm_t p, mm_t pinv)
    { return mm_redc(mm2_uconv(x,n),p,pinv); }

// computes o(x)=f(x)*g(x) of degree df+dg+1, aliasing is permitted (i.e. o=f or o=g is fine), prefers df <= dg
static inline mm_t *mm_poly_mul (mm_t o[], mm_t f[], int df, mm_t g[], int dg, mm_t p, mm_t pinv)
{
    register int i;
    if ( df > dg ) { SWAP(f,g); SWAP(df,dg); }
    for ( i = 0 ; i <= df ; i++ ) o[df+dg-i] = mm_conv (f+df-i, g+dg-i, i+1, p, pinv);
    for ( ; i <= dg ; i++ ) o[df+dg-i] = mm_conv (f, g+dg-i, df+1, p, pinv);
    for ( ; i <= df+dg ; i++ ) o[df+dg-i] = mm_conv (f, g, df+dg+1-i, p, pinv);
    return o;
}

// computes o(x)=f(x)*g(x) where both inputs and output are implicitly monic (the leading coefficient is not examined nor set)
// this code is not really any faster than mm_poly_mul
static inline mm_t *mm_poly_mul_m (mm_t o[], mm_t f[], int df, mm_t g[], int dg, mm_t R, mm_t p, mm_t pinv)
{
    register int i;
    if ( df > dg ) { SWAP(f,g); SWAP(df,dg); }
    if ( !df ) { memmove (o,g,dg*sizeof(g[0])); return o; }
    o[df+dg-1] = mm_add(f[df-1],g[dg-1],p);
    for ( i = 2 ; i <= df ; i++ ) o[df+dg-i] = mm_redc (mm2_conv(f+df-i+1, g+dg-i+1, i-1) + (mm2_t)R*(f[df-i]+g[dg-i]), p, pinv);
    for ( ; i <= dg ; i++ ) o[df+dg-i] =  mm_redc (mm2_conv(f, g+dg-i+1, df) + (mm2_t)R*g[dg-i], p, pinv);
    for ( ; i <= df+dg ; i++ ) o[df+dg-i] = mm_conv (f, g, df+dg+1-i, p, pinv);
    return o;
}

// computes o[i] = sum_j f[dg+i-j]*g[j] for i in [0,df-dg]
static inline mm_t *mm_poly_midmul (mm_t o[], mm_t f[], int df, mm_t g[], int dg, mm_t p, mm_t pinv)
    { if ( df < dg ) return o; mm_t t[df-dg+1]; for ( register int i = 0 ; i <= df-dg ; i++ ) t[i] = mm_conv (f+i, g, dg+1, p, pinv); memmove (o,t,(df-dg+1)*sizeof(mm_t)); return o; }

// sets o to f mad monic
static inline mm_t *mm_poly_make_monic (mm_t o[], mm_t f[], int d, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    if ( d < 0 ) return o;
    if ( ! d ) { o[0] = R; return o; }
    if ( f[d] == R ) { memmove (o,f,(d+1)*sizeof(f[0])); return o; }
    mm_smul (o,mm_inv(f[d],R2,R3,p,pinv),f,d,p,pinv); f[d] = R;
    return o;
}


// o(x) = f(x+a)
static inline mm_t *mm_poly_translate (mm_t o[], mm_t f[], int d, mm_t a, mm_t R, mm_t p, mm_t pinv)
{   
    if ( o != f ) memmove(o,f,(d+1)*sizeof(f[0]));
    register mm_t w = mm_mul(o[d], a, p, pinv);
    for ( register int i = d-1 ; i >= 0 ; i-- ) { for ( register int j = i ; j < d-1 ; j++ ) { o[j] = mm_addmul (o[j], o[j+1], a, R, p, pinv); } o[d-1] = mm_add(o[d-1], w, p); }
    return o;
}

// depresses monic f by translating, requires deg(f) nonzero mod p, sets a so that o(x+a) = f(x)
static inline mm_t *mm_poly_depress_monic (mm_t o[], mm_t a[1], mm_t f[], int d, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    if ( ! f[d-1] ) { if ( a ) a[0] = 0; memmove (o,f,(d+1)*sizeof(f[0])); return o; }
    mm_t z = mm_mul(mm_inv(mm_from_ui(d,R2,p,pinv),R2,R3,p,pinv),f[d-1],p,pinv);
    mm_poly_translate (o, f, d, mm_neg(z,p), R, p, pinv);
    if ( a ) a[0] = z;
    return 0;
}

static inline mm_t *mm_poly_depress_and_make_monic (mm_t o[], mm_t a[1], mm_t f[], int d, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    register mm_t x,y,z;

    if ( f[d] == R ) { return mm_poly_depress_monic(o,a,f,d,R,R2,R3,p,pinv); }
    if ( !f[d-1] ) { if ( a ) a[0] = 0; return mm_poly_make_monic (o,f,d,R,R2,R3,p,pinv); }
    x = mm_from_ui(d,R2,p,pinv);
    y = mm_inv(mm_mul(x,f[d],p,pinv),R2,R3,p,pinv);
    mm_smul (o,mm_mul(x,y,p,pinv),f,d,p,pinv);   
    z = mm_mul(f[d-1],mm_mul(f[d],y,p,pinv),p,pinv);
    f[d] = R;
    mm_poly_translate (o, f, d, mm_neg(z,p), R, p, pinv);
    if ( a ) a[0] = z;
    return 0;
}

// o(x) = prod_{i=1}^d (x-r[i])
void mm_poly_from_roots (mm_t o[], mm_t r[], int d, mm_t R, mm_t p, mm_t pinv);

// o(x) = f(x)*(x+a)
static inline mm_t *mm_poly_mul_xpa (mm_t o[], mm_t f[], int d, mm_t a, mm_t R, mm_t p, mm_t pinv)
    { if ( d < 0 ) return o;  for ( o[d+1] = f[d] ; d ; d-- ) o[d] = mm_addmul (f[d-1], f[d], a, R, p, pinv);  o[0] = mm_mul(f[0],a,p,pinv); return o; }

// Given f(x) such that f(r)=0, computes o(x) = f(x)/(x-r)
static inline mm_t *mm_poly_remove_root (mm_t o[], mm_t f[], int d, mm_t r, mm_t R, mm_t p, mm_t pinv)
    { register mm_t t0, t1;  for ( t0 = f[d--] ; d ; d-- ) { t1 = mm_addmul (f[d], r, t0, R, p, pinv); o[d] = t0; t0 = t1; }  o[0] = t0; return o; }

// o(x) = f(x)*(x-r)
static inline mm_t *mm_poly_add_root (mm_t o[], mm_t f[], int d, mm_t r, mm_t R, mm_t p, mm_t pinv)
    { return mm_poly_mul_xpa (o, f, d, mm_neg(r, p), R, p, pinv); }

// returns f(a)
static inline mm_t mm_poly_eval (mm_t f[], int d, mm_t a, mm_t R, mm_t p, mm_t pinv)
    { if ( d < 0 ) return 0;  register mm_t y = f[d];  for ( register int i = d-1 ; i >= 0 ; i-- ) y = mm_addmul(f[i], a, y, R, p, pinv);  return y; }

// computes f' (zero filled to degree d-1)
static inline mm_t *mm_poly_derivative (mm_t o[], mm_t f[], int d, mm_t R, mm_t p, mm_t pinv)
{
    if ( d < 1 ) return o;
    o[0] = f[1];
    mm_t x = mm_add(R,R,p);
    for ( register int i = 1 ; i < d ; i++ ) { o[i] = mm_mul (f[i+1],x,p,pinv); x = mm_add(x,R,p); }
    return o;
}    
    
// computes of o(x) = f(x) mod g(x) where deg(f)=n and g(x) = x^d - g[d-1]x^(d-1) - ... - g[1]x - g[0] (note that signs and g is implicitly monic, g[d] is never examined)
// we don't assume g is depressed (doing so really doesn't speed things up in the generic case where n is arbitrary)
// o should have space for d entries and will be zero extended to degree d-1.
// Complexity is approximately 2d(n-d)M +_(2d(n-d)-n)A + nR
static inline void mm_poly_mod_m (mm_t o[], mm_t f[], int n, mm_t g[], int d, mm_t R, mm_t p, mm_t pinv)
{
    if ( n < d ) { memmove(o,f,(n+1)*sizeof(mm_t)); memset (o+n+1,0,(d-n-1)*sizeof(mm_t)); return; }
    
    mm_t t[n-d+1];
    register int m;
    t[n-d] = f[n];
    for ( m = n-1 ; m >= d && m >= n-d ; m-- ) t[m-d] = mm_redc(mm2_conv(g+d+m-n, t+m-d+1, n-m) + (mm2_t)R*f[m],p,pinv);
    for ( ; m >= d ; m-- ) t[m-d] = mm_redc(mm2_conv(g, t+m-d+1, d) + (mm2_t)R*f[m],p,pinv);           // only relevant for n > 2d
    for ( ; m >= n-d ; m-- ) o[m] = mm_redc(mm2_conv(g+d+m-n, t, n-d+1) + (mm2_t)R*f[m],p,pinv);       // only relevant for n < 2d
    for ( ; m >= 0 ; m-- ) o[m] =  mm_redc(mm2_conv(g, t, m+1) + (mm2_t)R*f[m],p,pinv);                // only relevant for d > 0 (we do handle the stupid case d=0 but don't optimized for it)
}

// given monic g of degree d and f of degree n computes q and r such that f = q*g+r with deg(r)<deg(g), r zero-padded to degree d-1
static inline void mm_poly_div_m (mm_t q[], mm_t r[], mm_t f[], int n, mm_t g[], int d, mm_t R, mm_t p, mm_t pinv)
{
    register int m;

    if ( n < d ) { if ( r ) { memmove(r,f,(n+1)*sizeof(mm_t)); memset (r+n+1,0,(d-n-1)*sizeof(mm_t)); } return; }
    
    mm_t t[n-d+1];
    if ( ! q ) q = t;
    q[n-d] = f[n];
    for ( m = n-1 ; m >= d && m >= n-d ; m-- ) q[m-d] = mm_redc(mm2_conv(g+d+m-n, q+m-d+1, n-m) + (mm2_t)R*f[m],p,pinv);
    for ( ; m >= d ; m-- ) q[m-d] = mm_redc(mm2_conv(g, q+m-d+1, d) + (mm2_t)R*f[m],p,pinv);           // only relevant for n > 2d
    if ( ! r ) return;
    for ( ; m >= n-d ; m-- ) r[m] = mm_redc(mm2_conv(g+d+m-n, q, n-d+1) + (mm2_t)R*f[m],p,pinv);       // only relevant for n < 2d
    for ( ; m >= 0 ; m-- ) r[m] =  mm_redc(mm2_conv(g, q, m+1) + (mm2_t)R*f[m],p,pinv);                // only relevant for d > 0 (we do handle the stupid case d=0 but don't optimized for it)
}

// given monic g of degree d and f of degree n computes q and r such that f = q*g+r with deg(r)<deg(g), r zero-padded to degree d-1
static inline void mm_poly_div_exact_m (mm_t q[], mm_t f[], int n, mm_t g[], int d, mm_t R, mm_t p, mm_t pinv)
{
    register int m;

    if ( n < d ) return;
    
    q[n-d] = f[n];
    for ( m = n-1 ; m >= d && m >= n-d ; m-- ) q[m-d] = mm_redc(mm2_conv(g+d+m-n, q+m-d+1, n-m) + (mm2_t)R*f[m],p,pinv);
    for ( ; m >= d ; m-- ) q[m-d] = mm_redc(mm2_conv(g, q+m-d+1, d) + (mm2_t)R*f[m],p,pinv);           // only relevant for n > 2d
}

// computes o=gcd(f,g), returns degree (or a negative value if both f and g are zero), prefers df >= dg, optimized for df = dg+1
int mm_poly_gcd (mm_t o[], mm_t f[], int df, mm_t g[], int dg, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv);

// computes o(x)=f1(x)*f2(x) mod g(x) where g(x)=x^d - g[d-2]x^{d-2} - ... -g[1]x - g_0 (note signs and that g is monic, depressed, and non-constant)
// assumes f1 and f2 are reduced modulo g with coefficients zero-padded to degree d-1.
static inline void mm_poly_mul_mod_md (mm_t o[], mm_t f1[], mm_t f2[], mm_t g[], int d, mm_t p, mm_t pinv)
{
    mm_t t[d-1];

    // calling mm2_conv with n <= 0 is harmless (even with bogus pointers)
    for ( register int i = d-2 ; i >=0 ; i-- ) t[i] = mm_redc (mm2_conv (f1+i+1, f2+i+1, d-i-1) + mm2_conv(t+i+2, g+i+2, d-i-3), p, pinv);
    o[d-1] = mm_redc (mm2_conv(f1,f2,d) + mm2_conv(t+1,g+1,d-2), p, pinv);
    for ( register int i = d-2 ; i >= 0 ; i-- ) o[i] = mm_redc (mm2_conv(f1, f2, i+1) + mm2_conv(t, g, i+1), p, pinv);
}

// computes o(x)=f1(x)*f2(x) mod g(x) where g(x)=x^d - g[d=-1]x^{d-1} - ... -g[1]x - g_0 (note signs and that g is monic, not nescessarily depressed, and non-constant)
// assumes f1 and f2 are reduced modulo g with coefficients zero-padded to degree d-1.
static inline void mm_poly_mul_mod_m (mm_t o[], mm_t f1[], mm_t f2[], mm_t g[], int d, mm_t p, mm_t pinv)
{
    mm_t t[d-1];

    // calling mm2_conv with n <= 0 is harmless (even with bogus pointers)
    for ( register int i = d-2 ; i >=0 ; i-- ) t[i] = mm_redc (mm2_conv (f1+i+1, f2+i+1, d-i-1) + mm2_conv(t+i+1, g+i+2, d-i-2), p, pinv);
    o[d-1] = mm_redc (mm2_conv(f1,f2,d) + mm2_conv(t,g+1,d-1), p, pinv);
    for ( register int i = d-2 ; i >= 0 ; i-- ) o[i] = mm_redc (mm2_conv(f1, f2, i+1) + mm2_conv(t, g, i+1), p, pinv);
}

// computes o(x)=f(x)^2 where deg(f)=deg(g)=d (or zero-padded to deg d, note d = #coeffs + 1)
static inline void mm_poly_sqr (mm_t o[], mm_t f[], int d, mm_t p, mm_t pinv)
{
    for ( register int i = d ; i >=0 ; i-- ) o[d+i] = mm_uconv (f+i, d-i+1, p, pinv);
    for ( register int i = d-1 ; i >= 0 ; i-- ) o[i] = mm_uconv(f, i+1, p, pinv);
}

// computes o(x)=f(x)^2 mod g(x) where g(x)=x^d - g[d-2]x^{d-2} - ... -g[1]x - g_0 (note signs and that g is monic, depressed, and non-constant)
// assumes f is reduced modulo g with coefficients zero-padded to degree d-1.
static inline void mm_poly_sqr_mod_md (mm_t o[], mm_t f[], mm_t g[], int d, mm_t p, mm_t pinv)
{
    mm_t t[d-1];

    // calling mm2_conv with n <= 0 is harmless (even with bogus pointers)
    for ( register int i = d-2 ; i >=0 ; i-- ) t[i] = mm_redc (mm2_uconv (f+i+1, d-i-1) + mm2_conv(t+i+2, g+i+2, d-i-3), p, pinv);
    o[d-1] = mm_redc (mm2_uconv(f,d) + mm2_conv(t+1,g+1,d-2), p, pinv);
    for ( register int i = d-2 ; i >= 0 ; i-- ) o[i] = mm_redc (mm2_uconv(f, i+1) + mm2_conv(t, g, i+1), p, pinv);
}

// computes o(x)=xf(x)^2 mod g(x) where g(x)=x^d - g[d-2]x^{d-2} - ... -g[1]x - g_0 (note signs and that g is monic, depressed, and non-constant)
// assumes f is reduced modulo g with coefficients zero-padded to degree d-1.
static inline void mm_poly_sqr_mul_x_mod_md (mm_t o[], mm_t f[], mm_t g[], int d, mm_t p, mm_t pinv)
{
    mm_t t[d];

    // calling mm2_conv with n <= 0 is harmless (even with bogus pointers)
    for ( register int i = d-1 ; i >=0 ; i-- ) t[i] = mm_redc (mm2_uconv (f+i, d-i) + mm2_conv(t+i+2, g+i+1, d-i-2), p, pinv);
    o[d-1] = mm_redc (mm2_uconv(f,d-1) + mm2_conv(t+1,g,d-1), p, pinv);
    for ( register int i = d-2 ; i >= 0 ; i-- ) o[i] = mm_redc (mm2_uconv(f, i) + mm2_conv(t, g, i+1), p, pinv);
}


// register versions of dot, conv, uconv for small fixed n

static inline mm_t mm_dot_2r (mm_t x0, mm_t x1, mm_t y0, mm_t y1, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y0+(mm2_t)x1*y1,p,pinv); }
static inline mm_t mm_conv_2r (mm_t x0, mm_t x1, mm_t y0, mm_t y1, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y1+(mm2_t)x1*y0,p,pinv); }
static inline mm_t mm_uconv_2r (mm_t x0, mm_t x1, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*(2*x1),  p,  pinv); }

static inline mm_t mm_dot_3r (mm_t x0, mm_t x1, mm_t x2, mm_t y0, mm_t y1, mm_t y2, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y0+(mm2_t)x1*y1+(mm2_t)x2*y2,p,pinv); }
static inline mm_t mm_conv_3r (mm_t x0, mm_t x1, mm_t x2, mm_t y0, mm_t y1, mm_t y2, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y2+(mm2_t)x1*y1+(mm2_t)x2*y0,p,pinv); }
static inline mm_t mm_uconv_3r (mm_t x0, mm_t x1, mm_t x2, mm_t p, mm_t pinv)
    { return mm_redc(2*(mm2_t)x0*x2+(mm2_t)x1*x1,p,pinv); }

static inline mm_t mm_dot_4r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y0+(mm2_t)x1*y1+(mm2_t)x2*y2+(mm2_t)x3*y3,p,pinv); }
static inline mm_t mm_conv_4r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y3+(mm2_t)x1*y2+(mm2_t)x2*y1+(mm2_t)x3*y0,p,pinv); }
static inline mm_t mm_uconv_4r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t p, mm_t pinv)
    { return mm_redc(2*((mm2_t)x0*x3+(mm2_t)x1*x2),p,pinv); }

static inline mm_t mm_dot_5r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t y4, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y0+(mm2_t)x1*y1+(mm2_t)x2*y2+(mm2_t)x3*y3+(mm2_t)x4*y4,p,pinv); }
static inline mm_t mm_conv_5r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t y4, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y4+(mm2_t)x1*y3+(mm2_t)x2*y2+(mm2_t)x3*y1+(mm2_t)x4*y0,p,pinv); }
static inline mm_t mm_uconv_5r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t p, mm_t pinv)
    { return mm_redc(2*((mm2_t)x0*x4+(mm2_t)x1*x3)+(mm2_t)x2*x2,p,pinv); }
    
static inline mm_t mm_dot_6r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t x5, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t y4, mm_t y5, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y0+(mm2_t)x1*y1+(mm2_t)x2*y2+(mm2_t)x3*y3+(mm2_t)x4*y4+(mm2_t)x5*y5,p,pinv); }
static inline mm_t mm_conv_6r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t x5, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t y4, mm_t y5, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y5+(mm2_t)x1*y4+(mm2_t)x2*y3+(mm2_t)x3*y2+(mm2_t)x4*y1+(mm2_t)x5*y0,p,pinv); }
static inline mm_t mm_uconv_6r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t x5, mm_t p, mm_t pinv)
    { return mm_redc(2*((mm2_t)x0*x5+(mm2_t)x1*x4+(mm2_t)x2*x3),p,pinv); }

static inline mm_t mm_dot_7r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t x5, mm_t x6, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t y4, mm_t y5, mm_t y6, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y0+(mm2_t)x1*y1+(mm2_t)x2*y2+(mm2_t)x3*y3+(mm2_t)x4*y4+(mm2_t)x5*y5+(mm2_t)x6*y6,p,pinv); }
static inline mm_t mm_conv_7r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t x5, mm_t x6, mm_t y0, mm_t y1, mm_t y2, mm_t y3, mm_t y4, mm_t y5, mm_t y6, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)x0*y6+(mm2_t)x1*y5+(mm2_t)x2*y4+(mm2_t)x3*y3+(mm2_t)x4*y2+(mm2_t)x5*y1+(mm2_t)x6*y0,p,pinv); }
static inline mm_t mm_uconv_7r (mm_t x0, mm_t x1, mm_t x2,mm_t x3, mm_t x4, mm_t x5, mm_t x6, mm_t p, mm_t pinv)
    { return mm_redc(2*((mm2_t)x0*x6+(mm2_t)x1*x5+(mm2_t)x2*x4)+(mm2_t)x3*x3,p,pinv); }

// disc(f2*x^2+f1*x+f0) = f[1]^2-4*f[0]*f[2]
static inline mm_t mm_poly_disc_2 (mm_t f[3], mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)f[1]*f[1]+(mm2_t)(p-f[2])*(4*f[0]),p,pinv); }

// disc(x^3+f1*x+f0) = -4*f1^3 - 27*f0^2
static inline mm_t mm_poly_disc_md_3 (mm_t f[2], mm_t p, mm_t pinv)
    { register mm_t D = mm_sqr (f[1],p,pinv);  return mm_neg(mm_redc(27*(mm2_t)f[0]*f[0]+(mm2_t)f[1]*(D<<2),p,pinv),p); }

// disc(x^3+f2*x^2+f1*x+f0) = -4*f1^3 - 27*f0^2 + 18*f0*f1*f2 + f1^2*f2^2 -4*f0*f2^3
static inline mm_t mm_poly_disc_m_3 (mm_t f[3], mm_t p, mm_t pinv)
{
    register mm_t t1 = mm_sqr(f[1],p,pinv), t2 = mm_sqr(f[2],p,pinv), t3 = mm_mul(f[0],f[2],p,pinv);
    return mm_redc ((mm2_t)t1*(4*(p-f[1])) + 27*(mm2_t)f[0]*(p-f[0]) + 18*(mm2_t)t3*f[1] + (mm2_t)t2*(t1+4*(p-t3)),p,pinv); // avoids overflow for p < 2^61
}

// disc(f3*x^3+f2*x^2+f1*x+f0) = -4*f1^3*f3 - 27*f0^2*f3^2 + 18*f0*f1*f2*f3 + f1^2*f2^2 -4*f0*f2^3
static inline mm_t mm_poly_disc_3 (mm_t f[4], mm_t p, mm_t pinv)
{
    register mm_t t1 = mm_sqr(f[1],p,pinv), t2 = mm_sqr(f[2],p,pinv), t3 = mm_mul(f[0],f[3],p,pinv), t4 = mm_mul(f[0],f[2],p,pinv), t5 = mm_mul(f[1],f[3],p,pinv);
    return mm_redc ((mm2_t)t1*(4*(p-t5)) + 27*(mm2_t)t3*(p-t3) + 18*(mm2_t)t4*t5 + (mm2_t)t2*(t1+4*(p-t4)),p,pinv); // avoids overflow for p < 2^61
}


// computes bottom two coefficients of f(x-f2/3),given precomputed value 1/3 in Montgomery rep (assumes p != 3)
// without precomputation cost increases by about 4 clock cycles using mul_ui with third = ((3-(p%3))*p+1)/3
static inline mm_t mm_poly_depress_m_3 (mm_t g[2], mm_t f[3], mm_t third, mm_t R, mm_t p, mm_t pinv)
{
    register mm_t c = mm_mul (f[2],third,p,pinv), c2 = mm_sqr(c,p,pinv);
    g[0] = mm_redc ((mm2_t)R*f[0] + (mm2_t)c*(p-f[1]) + (mm2_t)(2*c)*c2,p,pinv);
    g[1] = mm_redc ((mm2_t)R*f[1] + (mm2_t)c*(p-f[2]),p,pinv);
    return c;
}

// computes bottom three coefficients of f(x-f3/4)
static inline mm_t mm_poly_depress_m_4 (mm_t g[4], mm_t f[4], mm_t R, mm_t p, mm_t pinv)
{
    register mm_t c = mm_div2(mm_div2(f[3],p),p), c2 = mm_sqr(c,p,pinv);
    g[0] = mm_redc ((mm2_t)R*f[0] + (mm2_t)(p-c)*f[1] + (mm2_t)c2*(f[2]+3*(p-c2)),p,pinv);
    g[1] = mm_redc ((mm2_t)R*f[1] + (mm2_t)(2*c)*(4*c2+p-f[2]),p,pinv);
    g[2] = mm_redc ((mm2_t)R*(f[2]+(6*(p-c2))),p,pinv);
    return c;
}

// disc(x^4+f2*x^2+f1*x+f0) = f1^2*(4*f2*(36*f0-f2^2)-27*f1^2) + 16*f0*(f2^2*(f2^2-8*f0)+16*f0^2)
static inline mm_t mm_poly_disc_md_4 (mm_t f[3], mm_t R, mm_t p, mm_t pinv)
{
    register mm_t t1 = mm_sqr (f[1],p,pinv), t2 = mm_sqr (f[2],p,pinv);
    register mm_t t3 = mm_redc (36*(mm2_t)f[0]*f[2]+(mm2_t)(p-f[2])*t2,p,pinv);                 // t3 = f2*(36*f0-f2^2)
    register mm_t t4 = mm_redc ((mm2_t)t2*(t2+8*(mm2_t)(p-f[0]))+16*(mm2_t)f[0]*f[0],p,pinv);   // t4 = f2^2*(f2^2-8*f0)+16*f0^2
    return mm_redc ((mm2_t)t1*(4*t3+27*(mm2_t)(p-t1)) + 16*(mm2_t)f[0]*t4,p,pinv);              // t1*(4*t3-27*t1) + 16*f0*t4
}

// disc(x^4+f3*x^3+f2*x^2+f1*x+f0)
static inline mm_t mm_poly_disc_m_4 (mm_t f[4], mm_t R, mm_t p, mm_t pinv)
    { mm_t g[3];  mm_poly_depress_m_4 (g, f, R, p, pinv);  return mm_poly_disc_md_4 (g, R, p, pinv); }

// computes ad-bc
static inline mm_t mm_det_2r (mm_t a, mm_t b, mm_t c, mm_t d, mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)a*d+(mm2_t)b*(p-c),p,pinv); }

static inline mm_t mm_det_2 (mm_t m[4], mm_t p, mm_t pinv)
    { return mm_redc((mm2_t)m[0]*m[3]+(mm2_t)m[1]*(p-m[2]),p,pinv); }

static inline mm_t mm_det_3r (mm_t a, mm_t b, mm_t c, mm_t d, mm_t e, mm_t f, mm_t g, mm_t h, mm_t i, mm_t p, mm_t pinv)
    { return mm_dot_3r (a,d,g,mm_det_2r(e,f,h,i,p,pinv),mm_det_2r(c,b,i,h,p,pinv),mm_det_2r(b,c,e,f,p,pinv),p,pinv); }
    
static inline mm_t mm_det_3 (mm_t m[9], mm_t p, mm_t pinv)
    { return mm_det_3r (m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],p,pinv); }
    
static inline mm_t mm_det_4 (mm_t m[16], mm_t p, mm_t pinv)
{
    return mm_dot_4r (m[0],m[4],m[8],m[12], mm_det_3r(m[5],m[6],m[7],m[9],m[10],m[11],m[13],m[14],m[15],p,pinv),
                                            mm_det_3r(m[1],m[2],m[3],m[13],m[14],m[15],m[9],m[10],m[11],p,pinv),
                                            mm_det_3r(m[1],m[2],m[3],m[5],m[6],m[7],m[13],m[14],m[15],p,pinv),
                                            mm_det_3r(m[1],m[2],m[3],m[9],m[10],m[11],m[5],m[6],m[7],p,pinv), p,pinv);
}
    
// destructively computes determinant of n x n matrix A
mm_t _mm_det (mm_t m[], int n, mm_t R, mm_t p, mm_t pinv);
static inline mm_t mm_det  (mm_t m[], int n, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    switch (n) {
    case 1: return m[0];
    case 2: return mm_det_2 (m,p,pinv);
    case 3: return mm_det_3 (m,p,pinv);
    case 4: return mm_det_4 (m,p,pinv);
    }
    mm_t x = _mm_det (m, n, R, p, pinv);
    return mm_mul (*m,mm_inv(x,R2,R3,p,pinv),p,pinv);
}

// returns degree of poly f given degree bound d
static inline int mm_poly_deg (mm_t f[], int d) { while ( d >= 0 && !f[d] ) d--;  return d; }

//
static inline int mm_poly_squarefree (mm_t f[], int d, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    assert ( d > 0 );
    mm_t g[d];
    mm_poly_derivative (g,f,d,R,p,pinv);
    return mm_poly_gcd (g,f,d,g,mm_poly_deg(g,d-1),R,R2,R3,p,pinv) == 0 ? 1 : 0;
}

// computes Res(f,g) as o[0]/o[1] with o[1] nonzero (caller must invert o[1] to get exact value)
void _mm_poly_res (mm_t o[2], mm_t *f, int df, mm_t *g, int dg, mm_t R, mm_t p, mm_t pinv);

// returns Res(f,g)
static inline mm_t mm_poly_res (mm_t *f, int df, mm_t *g, int dg, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
    { mm_t o[2]; _mm_poly_res (o,f,df,g,dg,R,p,pinv); return ( o[1] == R ? o[0] : mm_mul(o[0],mm_inv(o[1],R2,R3,p,pinv),p,pinv) ); }

static inline mm_t mm_poly_disc (mm_t f[], int d, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    assert (d > 0);
    switch (d) {
    case 1: return R;
    case 2: return mm_poly_disc_2 (f,p,pinv);
    case 3: if ( f[3] == R ) return f[2] ? mm_poly_disc_m_3 (f,p,pinv) : mm_poly_disc_md_3 (f,p,pinv); else return mm_poly_disc_3 (f,p,pinv);
    case 4: if ( f[4] == R ) return f[3] ? mm_poly_disc_m_4 (f,R,p,pinv) : mm_poly_disc_md_4 (f,R,p,pinv);
            { mm_t g[4]; if ( f[4] != R ) mm_smul(g,mm_inv(f[4],R2,R3,p,pinv),f,4,p,pinv);  mm_t x = mm_exp_ui(f[4],6,R,p,pinv);
              return mm_mul(x,f[3] ? mm_poly_disc_m_4 (g,R,p,pinv) : mm_poly_disc_md_4 (g,R,p,pinv),p,pinv); }
    }
    
    mm_t o[2], g[d];
    int e;
    
    mm_poly_derivative (g,f,d,R,p,pinv); e = mm_poly_deg(g,d-1);
    _mm_poly_res (o, f, d, g, e, R, p, pinv);
    if ( ! o[0] ) return 0;
    if ( (d&2) ) o[0] = mm_neg(o[0],p); // (-1)^(d*(d-1)/2) = -1 iff d=2,3 mod 4
    e = d-e-2;
    if ( e > 0 && f[d] != R ) o[0] = mm_mul(o[0],mm_exp_ui(f[d],e,R,p,pinv),p,pinv);
    if ( e < 0 && f[d] != R ) o[1] = mm_mul(o[1],mm_exp_ui(f[d],-e,R,p,pinv),p,pinv);
    return o[1] == R ? o[0] : mm_mul(o[0],mm_inv(o[1],R2,R3,p,pinv),p,pinv);
}

int mm_poly_res_legendre (mm_t f[], int df, mm_t g[], int dg, mm_t R, mm_t p, mm_t pinv);

static inline mm_t mm_poly_disc_legendre (mm_t f[], int d, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    assert (d > 0);
    if ( d <= 4 ) return mm_legendre (mm_poly_disc(f,d,R,R2,R3,p,pinv),R,p,pinv);
    
    mm_t g[d],o[2];
    int k,e;
    
    mm_poly_derivative (g,f,d,R,p,pinv); e = mm_poly_deg(g,d-1);
    _mm_poly_res(o,f,d,g,e,R,p,pinv);
    k = mm_legendre (mm_mul(o[0],o[1], p, pinv), R , p, pinv);
    if ( ! k ) return 0;
    if ( (d&2) && (p&3)==3 ) k = -k;                    // account for factor of (-1)^(d*(d-1)/2), only relevant if p = 3 mod 4
    if ( ((d^e)&1) ) k *= mm_legendre (f[d],R,p,pinv);    // account for factor of f[d]^(d-e-2), only relevant if f[d] is not square
    return k;
}

// computes the square root of 1/x modulo p=2^e*m+1 (with m odd), given lists of binary powers of generator a for the 2-Sylow and ainv=1/a
// returns 0 if x does not have a square-root in Fp (if ext is set then sqrt(x) = z*sqrt(a[0]) in Fp2)
int mm_invsqrt (mm_t *z, mm_t x, mm_t a[], mm_t ainv[], int e, mm_t R, mm_t p, mm_t pinv, int ext);

static inline int mm_sqrt (mm_t *z, mm_t x, mm_t a[], mm_t ainv[], int e, mm_t R, mm_t p, mm_t pinv)
    { if ( ! mm_invsqrt(z,x,a,ainv,e,R,p,pinv,0) ) return 0;  *z = mm_mul (*z,x,p,pinv); return 1; }


#ifdef __cplusplus
}
#endif

#endif
