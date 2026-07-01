#ifndef _NTUTIL_INCLUDE_
#define _NTUTIL_INCLUDE_

/*
    Copyright (c) 2007-2019 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <assert.h>
#include <math.h>
#include "cstd.h"

/*
    Elementary number-theoretic functions implemented as inlines
*/

// Algorithm 1.4.10, p. 29 of [CANT] trimmed down to handle a >= 0, b > 0 odd (for the sake of space and speed)
static inline int ui_legendre (unsigned long a, unsigned long b)
{
    register unsigned long r;
    register int k,v;
    
    if ( a > 2*b ) a %= b;
    while ( a >= b ) a -= b;
    k = 1;
    while ( a ) {
        v = __builtin_ctzl(a); a >>= v;
        if ( (v&1) && (b&1) && ((b&6) == 2 || (b&6) == 4) ) k = -k;
        if ( (a&b&3) == 3 ) k = -k;
        r = a;   a = (b < 2*r ? b-r : b%r);  b = r;
    }
    return ( b == 1 ? k : 0 );
}

// Given an odd prime p (not checked!) returns a non-residue mod p
// This takes less than 10ns on average
static inline unsigned long ui_nonresidue (unsigned long p)
{
    if ((p&3)==3) return p-1;
    if ((p&7)==5) return 2;
    if ((p%3)==2) return p-3;
    if ((p%5)&2) return 5;
    if (__builtin_popcount(p%7)>1 ) return p-7;
    for ( int n = 11 ; ; n += 2 ) if ( ui_legendre(n,p) < 0 ) return n;
}

// computes the legendre symbol for any integer a and odd prime b (in fact b can be any positive odd number)
static inline int legendre (long a, long b)
{
    register int s;
    
    s = 1;
    if ( a < 0 ) { s = ((b&3)==1 ? 1 : -1); a = -a; }
    return s*ui_legendre (a,b);
}


// Algorithm 1.4.10, p. 29 of [CANT] for a,b >= 0
static inline int ui_kronecker (unsigned long a, unsigned long b)
{
    register unsigned long r;
    register int k, v;
    
    if ( ! b ) return ( a==1 ? 1 : 0 );
    k = 1;
    if ( !(b&1) ) {
        if ( !(a&1) ) return 0;
        v = __builtin_ctzl(b); b >>= v;
        if ( (v&1) && ( (a&6) == 2 || (a&6) == 4 ) ) k = -1;
    }
    if ( a > b ) a %= b;
    while ( a ) {
        v = __builtin_ctzl(a); a >>= v;
        if ( (v&1) && (b&1) && ((b&6) == 2 || (b&6) == 4) ) k = -k;
        if ( (a&b&3) == 3 ) k = -k;
        r = a;   a = (b < (r+r) ? b-r : b%r);  b = r;
    }
    return ( b == 1 ? k : 0 );
}

// computes the kronecker symbol for arbitrary integers a and b
static inline int kronecker (long a, long b)
{
    int s;
    
    if ( b < 0 ) { s = (a<0?-1:1); b = -b; } else s = 1;
    if ( a < 0 ) { long c; if ( ! b ) return (a==-1?1:0); for ( c = b ; !(c&1) ; c >>= 1 ); s *= ((c&3)==3 ? -1 : 1); a = -a; }
    return s*ui_kronecker (a,b);
}

long i_sqrt_modp (long a, long p);  // returns -1 if  a is not a square mod p, otherwise return value is a sqrt in [0,p-1]

int i_norm_equation (long *x, long *y, long d, long p, long pe); // finds solution to 4pe=x^2+y^2d where pe is a power of p and d is coprime to p

#define _mod64res2  0x202021202030213
#define _mod63res2  0x402483012450293
#define _mod25res2  0x1294a53
#define _mod13res2  0x161b
#define _mod27res3  0x40A0503
#define _mod7res3   0x43
#define _mod63res3  0x4080001818000103
#define _mod13res3  0x1123
#define _mod19res3  0x41983
#define _mod31res3  0x68818117
#define _mod37res3  0x10ac804d43
#define _mod11res5  0x403
#define _mod31res5  0x46000063
#define _mod41res5  0x1410800420b
#define _mod61res5  0x1005810120206803
#define _mod29res7  0x10021003
#define _mod43res7  0x430000000c3
#define _mod49res7  0x10000c00c0003
#define _mod19res9  0x40003
#define _mod27res9  0x4000003
#define _mod37res9  0x1080000043
#define _mod73res9  0x8008400008400403  // also need to check 72
#define _mod23res11 0x400003
#define _mod67res11 0x6060000003        // also need to check 66
#define _mod89res11 0x90002400001003    // also need to check 77 and 88
#define _mod53res13 0x10000040800003
#define _mod79res13 0x180000001800003   // also need to check 78
#define _mod103res17 0x300C00000000003  // also need to check 102
#define _mod137res17 0x22000000403      // also need to check -n

// return sqrt(n) if n is a perfect nonzero square, 0 ow.
static inline long i_sqrt (long n)
{
    if ( n <= 0 ) return 0;
    if ( ! (_mod64res2&(1UL<<(n&0x3F))) ) return 0;
    long x = sqrt(n)+0.1;        // avoid potential rounding problems
    return (x*x == n ? x : 0);
}

static inline unsigned long ui_sqrt (unsigned long n)
{    
    if ( ! (_mod64res2&(1UL<<(n&0x3F))) ) return 0;
    unsigned long x = sqrt(n)+0.1;        // avoid potential rounding problems
    return (x*x == n ? x : 0);
}

static inline uint64_t ui128_sqrt (uint128_t n)
{
    if ( ! (n>>64) ) return ui_sqrt(n);
    if ( ! (_mod64res2&(1UL<<(n&0x3F))) ) return 0;
    if ( ! (_mod63res2&(1UL<<(n%63))) ) return 0;
    if ( ! (_mod25res2&(1UL<<(n%25))) ) return 0;
    uint128_t x = sqrtl(n)+0.1;        // avoid potential rounding problems
    return (x*x == n ? x : 0);
}

static inline unsigned long ui_cbrt (unsigned long n)
{
    if ( ! (_mod63res3&(1UL<<(n%63))) ) return 0;
    if ( ! (_mod13res3&(1UL<<(n%13))) ) return 0;
    if ( ! (_mod19res3&(1UL<<(n%19))) ) return 0;
    if ( ! (_mod31res3&(1UL<<(n%31))) ) return 0;
    if ( ! (_mod37res3&(1UL<<(n%37))) ) return 0;
    // using reduction mod 43 does not help
    unsigned long x = cbrt(n)+0.1;
    return (x*x*x == n ? x : 0);
}

static inline uint64_t ui128_cbrt (uint128_t n)
{
    if ( !(n>>64) ) return ui_cbrt(n);
    if ( ! (_mod63res3&(1UL<<(n%63))) ) return 0;
    if ( ! (_mod13res3&(1UL<<(n%13))) ) return 0;
    if ( ! (_mod19res3&(1UL<<(n%19))) ) return 0;
    if ( ! (_mod31res3&(1UL<<(n%31))) ) return 0;
    if ( ! (_mod37res3&(1UL<<(n%37))) ) return 0;
    // using reduction mod 43 does not help
    uint64_t x = cbrtl(n)+0.1;
    return ((uint128_t)x*x*x == n ? x : 0);
}

static inline unsigned long ui_root4 (unsigned long n)
    { n = ui_sqrt(n); return n ? ui_sqrt(n) : 0; }
static inline uint64_t ui64_root4 (uint64_t n)
    { n = ui64_sqrt(n); return n ? ui64_sqrt(n) : 0; }

static inline uint64_t ui128_root4 (uint128_t n)
    { uint64_t m = ui128_sqrt(n); return m ? ui_sqrt(m) : 0; }

static inline unsigned long ui_root5 (long n)
{
    if ( ! (_mod11res5&(1UL<<(n%11))) ) return 0;
    if ( ! (_mod31res5&(1UL<<(n%31))) ) return 0;
    if ( ! (_mod41res5&(1UL<<(n%41))) ) return 0;
    if ( ! (_mod61res5&(1UL<<(n%61))) ) return 0;
    unsigned long x = pow(n,0.2)+0.1, y = x*x;
    return (y*y*x == n ? x : 0);
}

static inline uint64_t ui128_root5 (uint128_t n)
{
    if ( ! (n>>64) ) return ui_root5(n);
    if ( ! (_mod11res5&(1UL<<(n%11))) ) return 0;
    if ( ! (_mod31res5&(1UL<<(n%31))) ) return 0;
    if ( ! (_mod41res5&(1UL<<(n%41))) ) return 0;
    if ( ! (_mod61res5&(1UL<<(n%61))) ) return 0;
    uint64_t x = pow(n,0.2)+0.1, y = x*x;
    return (y*(uint128_t)y*x == n ? x : 0);
}

static inline unsigned long ui_root6 (unsigned long n)
    { n = ui_sqrt(n); return n ? ui_cbrt(n) : 0; }

static inline uint64_t ui128_root6 (uint128_t n)
    { uint64_t m = ui128_sqrt(n); return m ? ui_cbrt(m) : 0; }

static inline unsigned long ui_root7 (unsigned long n)
{
    if ( ! (_mod29res7&(1UL<<(n%29))) ) return 0;
    if ( ! (_mod43res7&(1UL<<(n%43))) ) return 0;
    if ( ! (_mod49res7&(1UL<<(n%49))) ) return 0;
    unsigned long x = pow(n,1.0/7)+0.1, y = x*x*x;
    return (y*y*x == n ? x : 0);
}

static inline uint128_t ui128_root7 (uint128_t n)
{
    if ( ! (n>>64) ) return ui_root7(n);
    if ( ! (_mod29res7&(1UL<<(n%29))) ) return 0;
    if ( ! (_mod43res7&(1UL<<(n%43))) ) return 0;
    if ( ! (_mod49res7&(1UL<<(n%49))) ) return 0;
    uint64_t x = pow(n,1.0/7)+0.1, y = x*x*x;
    return (y*(uint128_t)y*x == n ? x : 0);
}

static inline unsigned long ui_root8 (unsigned long n)
    { n = ui_sqrt(n); return n ? ui_root4(n) : 0; }

static inline uint128_t ui128_root8 (uint128_t n)
    { uint64_t m = ui128_root4(n); return m ? ui_sqrt(m) : 0; }

static inline unsigned long ui_root9 (unsigned long n)
{
    if ( ! (_mod19res9&(1UL<<(n%19))) ) return 0;
    if ( ! (_mod27res9&(1UL<<(n%27))) ) return 0;
    if ( ! (_mod37res9&(1UL<<(n%37))) ) return 0;
    int i = n%73; if ( !(_mod73res9&1UL<<i) && i != 72 ) return 0;
    unsigned long x = pow(n,1.0/9)+0.1, y = x*x, z = y*y;
    return (z*z*x == n ? x : 0);
}

static inline unsigned long ui128_root9 (uint128_t n)
{
    if ( ! (n>>64) ) return ui_root9(n);
    if ( ! (_mod19res9&(1UL<<(n%19))) ) return 0;
    if ( ! (_mod27res9&(1UL<<(n%27))) ) return 0;
    if ( ! (_mod37res9&(1UL<<(n%37))) ) return 0;
    int i = n%73; if ( !(_mod73res9&1UL<<i) && i != 72 ) return 0;
    uint64_t x = pow(n,1.0/9)+0.1, y = x*x, z = y*y;
    return (z*(uint128_t)z*x == n ? x : 0);
}

static inline unsigned long ui_root10 (unsigned long n)
    { n = ui_sqrt(n); return n ? ui_root5(n) : 0; }

static inline uint128_t ui128_root10 (uint128_t n)
    { uint64_t m = ui128_sqrt(n); return m ? ui_root5(m) : 0; }

static inline unsigned long ui_root11 (unsigned long n)
{
    if ( ! (_mod23res11&(1UL<<(n%23))) ) return 0;
    int i = n%67; if ( !(_mod67res11&(1UL<<i)) && i != 66 ) return 0;
    i = n%89; if ( !(_mod89res11&(1UL<<i)) && i != 77 && i != 88 ) return 0;
    unsigned long x = pow(n,1.0/11)+0.1, y = x*x, z=x*y*y;
    return (z*z*x == n ? x : 0);
}

static inline uint128_t ui128_root11 (uint128_t n)
{
    if ( ! (n>>64) ) return ui_root11(n);
    if ( ! (_mod23res11&(1UL<<(n%23))) ) return 0;
    int i = n%67; if ( !(_mod67res11&(1UL<<i)) && i != 66 ) return 0;
    i = n%89; if ( !(_mod89res11&(1UL<<i)) && i != 77 && i != 88 ) return 0;
    uint64_t x = pow(n,1.0/11)+0.1, y = x*x, z = x*y*y;
    return (z*(uint128_t)z*x == n ? x : 0);
}

static inline unsigned long ui_root12 (unsigned long n)
    { n = ui_sqrt(n); return n ? ui_root6(n) : 0; }

static inline uint64_t ui128_root12 (uint128_t n)
    { uint64_t m = ui128_sqrt(n); return m ? ui_root6(m) : 0; }

static inline unsigned long ui_root13 (unsigned long n)
{
    if ( ! (_mod53res13&(1UL<<(n%53))) ) return 0;
    int i = n%79; if ( !(_mod79res13&(1UL<<i)) && i != 78 ) return 0;
    unsigned long x = pow(n,1.0/13)+0.1, y = x*x*x, z = y*y;
    return (z*z*x == n ? x : 0);
}

static inline uint128_t ui128_root13 (uint128_t n)
{
    if ( ! (n>>64) ) return ui_root13(n);
    if ( ! (_mod53res13&(1UL<<(n%53))) ) return 0;
    int i = n%79; if ( !(_mod79res13&1UL<<i) && i != 78 ) return 0;
    uint64_t x = pow(n,1.0/13)+0.1, y = x*x*x, z = y*y;
    return (z*(uint128_t)z*x == n ? x : 0);
}

static inline unsigned long ui_root17 (unsigned long n)
{
    int i = n%103; if ( !(_mod103res17&(1UL<<i)) && i != 102 ) return 0;
    i = n%137;  if ( i > 64 ) i = 137-i; if ( !(_mod137res17&(1UL<<i)) ) return 0;
    unsigned long x = pow(n,1.0/17)+0.1, w = x*x, y = w*w, z = y*y;
    return (z*z*x == n ? x : 0);
}

static inline uint128_t ui128_root17 (uint128_t n)
{
    if ( ! (n>>64) ) return ui_root17(n);
    int i = n%103; if ( !(_mod103res17&(1UL<<i)) && i != 102 ) return 0;
    i = n%137;  if ( i > 64 ) i = 137-i; if ( !(_mod137res17&(1UL<<i)) ) return 0;
    uint64_t x = pow(n,1.0/17)+0.1, w = x*x, y = w*w, z = y*y;
    return (z*(uint128_t)z*x == n ? x : 0);
}

static inline unsigned long ui_root19 (unsigned long n)
{
    if ( n == 524288 ) return 2;
    if ( n == 1162261467 ) return 3;
    if ( n == 19073486328125 ) return 5;
    if ( n == 11398895185373143 ) return 7;
    return 0;
}

static char mod191res19[191] = {1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1};
static char mod229res19[229] = {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};

static inline uint128_t ui128_root19 (uint128_t n)
{
    if ( !mod191res19[n%191] ) return 0;
    if ( !mod229res19[n%229] ) return 0;
    uint64_t x = pow(n,1.0/19)+0.1, w = x*x, y = w*w, z = y*y*x ;
    return (z*(uint128_t)z*x == n ? x : 0);
}

// The ui_ppbase_12 functions below are more than 10 times as fast as calling mpz_perfect_power_p

static inline unsigned long ui_ppbase_12 (unsigned long n) // returns least m for which n = m^e with e > 1, assuming e < 12, equivalently, m > 40.
{
    unsigned long j,k,m;
    if ( (m=ui_sqrt(n)) ) {
        if ( (k=ui_sqrt(m)) ) {
            if ( (j=ui_sqrt(k)) ) return j;
            if ( (j=ui_cbrt(k)) ) return j;
            return k;
        }
        if ( (k=ui_cbrt(m)) ) return k;
        if ( (k=ui_root5(m)) ) return k;
        return m;
    }
    if ( (m=ui_cbrt(n)) ) return ( (k=ui_cbrt(m)) ? k : m );
    if ( (m = ui_root5(n)) ) return m;
    if ( (m = ui_root7(n)) ) return m;
    if ( (m = ui_root11(n)) ) return m;
    return 0;
}

static inline uint64_t ui128_ppbase_24 (uint128_t n) // returns least m for which n = m^e with e > 1, assuming e < 24, equivalently, m > 40.
{
    if (! (n>>64) ) return ui_ppbase_12(n);
    unsigned long i,j,k,m;
    if ( (m=ui128_sqrt(n)) ) {
        if ( (k=ui_sqrt(m)) ) {
            if ( (j=ui_sqrt(k)) ) return ( (i=ui_sqrt(j)) ? i : j);
            if ( (j=ui_cbrt(k)) ) return j;
            if ( (j=ui_root5(k)) ) return j;
            return k;
        }
        if ( (k=ui_cbrt(m)) )  return ( (j=ui_cbrt(k)) ? j : k);
        if ( (k=ui_root5(m)) ) return k;
        if ( (k=ui_root7(m)) ) return k;
        if ( (k=ui_root11(m)) ) return k;
        return m;
    }
    if ( (m=ui128_cbrt(n)) ) {
        if ( (k=ui_cbrt(m)) ) return k;
        if ( (k=ui_root5(m)) ) return k;
        if ( (k=ui_root7(m)) ) return k;
    }
    if ( (m = ui128_root5(n)) ) return m;
    if ( (m = ui128_root7(n)) ) return m;
    if ( (m = ui128_root11(n)) ) return m;
    int b = ui128_len(n);
    if ( b > 69 && (m = ui128_root13(n)) ) return m;
    if ( b > 91 && (m = ui128_root17(n)) ) return m;
    if ( b > 101 && (m = ui128_root19(n)) ) return m;
    if ( (n&0xFFFFFFFFFFFFFFFFUL) == 6681134104860204761UL && (n>>64) == 673145554550018685UL ) return 41;
    return 0;
}

#endif
