/*

zp_poly_main.h: public interface to zp_poly library

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

#ifndef ZP_POLY_MAIN_H
#define ZP_POLY_MAIN_H

// TODO: cleanup and organize this header file to improve readability

#include <gmp.h>

#define ZP_POLY_EXP_TIME_TRACE		0

/*
  The following applies to all of the functions defined below.
    1) p is an odd prime
    2) all mpz_t's hold integers in [0,p-1]
    3) degree d polys in Z/pZ[x] are represented as arrays
       of d+1 coefficients in [0,p-1], indexed by degree.
    4) in zp_array functions, array lengths (degree+1) are specified
    5) in zp_poly functions, poly degrees (length-1) are specified
    6) in zp_poly_mod functions, coefficient arrays all have length
       deg(g), where g is the modulus, and are zero-filled as required.
*/

// g = f, optimizes aliased case
static inline void zp_poly_copy (zp_t g[], zp_t f[], int d)
{
	register int i;
	
	if ( f != g ) for ( i = 0 ; i <= d ; i++ ) zp_set(g[i],f[i]);
}

// returns exact degree of f of degree <= d, zero poly has degree -1
static inline int zp_poly_degree (zp_t f[], int d)
{
	register int i;
	
	for ( i = d ; i >= 0 && zp_is_zero(f[i]) ; i-- );
	return i;
}

// ( f1 == f2 ? 1 : 0 )
static inline int zp_poly_equal (zp_t f1[], int d1, zp_t f2[], int d2)
{
	register int i;
	
	if ( d1 != d2 ) return 0;
	for ( i = 0 ; i <= d1 ; i++ ) if ( !zp_equal(f1[i],f2[i]) ) return 0;
	return 1;
}

// ( f == 0 ? 1 : 0 )
static inline int zp_poly_is_zero (zp_t f[], int d)
	{ return ( zp_poly_degree(f,d)==-1 ? 1 : 0 ); }


// h=0 padded to degree d
static inline void zp_poly_zero_fill (zp_t h[], int d)
{
	register int i;
	for ( i = 0 ; i <=d ; i++ ) zp_set_zero(h[i]);
}

// zero extend h from degree d1 to d2
static inline void zp_poly_zero_pad (zp_t h[], int d1, int d2)
{
	register int i;
	for ( i = d1 ; i <=d2 ; i++ ) zp_set_zero(h[i]);
}
	
// ( f1 == 1 ? 1 : 0 )
static inline int zp_poly_is_one (zp_t f[], int d)
	{ return ( zp_poly_degree(f,d)==0 && zp_is_one(f[0]) ); }

// returns zero if lc(f)!=1, nonzero o.w.  zero poly is (vacuously) monic
static inline int zp_poly_is_monic (zp_t f[], int d)
	{ return ( d >= 0 && zp_is_one(f[d]) ); }

// y = f(x), y may alias x
static inline void zp_poly_eval (zp_t y, zp_t f[], int d, zp_t x, zp_p p)
{
	zp_t w;
	register int i;

	if ( d < 0 ) { zp_set_zero(y); return; }
	zp_init_w (w, p);
	zp_set (w,f[d]);
	for ( i = d-1 ; i >= 0 ; i-- ) { mpz_mul(w,w,x); mpz_add(w,w,f[i]); mpz_mod(w,w,p); }
	zp_set (y, w);
	zp_clear (w);
}

// h = c*f
static inline void zp_poly_smul (zp_t h[], zp_t f[], int d, zp_t c, zp_p p)
{
	mpz_t w;
	register int i;

	if ( zp_is_zero(c) ) { zp_poly_zero_fill (h, d); return; }						// check for multiplication by 0
	if ( zp_is_one(c) ) { zp_poly_copy (h, f, d); return; }						// check for multiplication by 1
	zp_init_w (w, p);
	for ( i = d ; i >= 0 ; i-- ) { mpz_mul(w,c,f[i]); mpz_mod(h[i],w,(zp_p_to_mpz)p); }	// work top down (in case c is f[0])
	zp_clear (w);
}

// sets c = 1/lc(f) and replaces f with c*f
static inline int zp_poly_make_monic (zp_t f[], int d, zp_t c, zp_p p)
{
	if ( d < 0 ) return 1;								// zero poly is already monic (vacuously)
	if ( zp_is_one(f[d]) ) { zp_set_one(c); return 1; }			// handle 1 efficiently
	if ( ! zp_invert(c,f[d],p) ) return 0;						// with p prime, this fails only if lc(f)=0
	zp_set_one(f[d]);  zp_poly_smul (f,f,d-1,c,p);
	return 1;
}

// g = f'
static inline void zp_poly_derivative (zp_t g[], zp_t f[], int d, zp_p p)
{
	register int i;
	
	if ( d <= 0 ) return;
	// work bottom up in case g and f coincide
	zp_set (g[0], f[1]);
	for ( i = 1 ; i < d ; i++ ) { zp_mul_ui (g[i], f[i+1], i+1, p); }
}

// replace poly f of degree d by x^d*f(1/x)
static inline void zp_poly_rev_inplace (zp_t f[], int d)
{
	mpz_t x;
	register int i;
	
	// note that this doesn't actually move any limbs, it just swaps __mpz_t structs
	// one might consider this cheating, but it is fast...
	for ( i = 0 ; i < (d+1)/2 ; i++ ) { *x = *f[i]; *f[i] = *f[d-i]; *f[d-i] = *x; }
}

// h = f1 + f2, if dh is not null, *dh is set to deg(h)
static inline int zp_poly_add (zp_t h[], zp_t f1[], int d1, zp_t f2[], int d2, zp_p p)
{
	register int i,d;
	
	d = ( d1 < d2 ? d1 : d2 );
	for (  i = 0 ; i <= d ; i++ ) zp_add (h[i],f1[i],f2[i],p);
	for ( ; i <= d1 ; i++ ) zp_set(h[i],f1[i]);
	for ( ; i <= d2 ; i++ ) zp_set(h[i],f2[i]);
	return zp_poly_degree (h,i-1);
}

// h = 3 * f
static inline void zp_poly_triple (zp_t *h, zp_t *f, int d, zp_p p)
{
	zp_t w;
	register int i;
	
	if ( h == f ) {
		zp_init_w (w, p); 
		for ( i = 0 ; i <= d ; i++ ) { zp_add (w, f[i], f[i], p);  zp_add (h[i], f[i], w, p); }
		zp_clear (w);
	} else {
		for ( i = 0 ; i <= d ; i++ ) { zp_add (h[i], f[i], f[i], p);  zp_add (h[i], h[i], f[i], p); }
	}
}

// h = f1 - f2, if dh is not null, *dh is set to deg(h)
static inline int zp_poly_sub (zp_t h[], zp_t f1[], int d1, zp_t f2[], int d2, zp_p p)
{
	register int i,d;
	
	d = ( d1 < d2 ? d1 : d2 );
	for (  i = 0 ; i <= d ; i++ ) zp_sub (h[i],f1[i],f2[i],p);
	for ( ; i <= d1 ; i++ ) zp_set (h[i],f1[i]);
	for ( ; i <= d2 ; i++ ) zp_neg (h[i],f2[i],p);
	return zp_poly_degree (h,i-1);
}

// h = -h
static inline void zp_poly_neg (zp_t h[], zp_t f[], int d, zp_p p)
	{ register int i;  for ( i = 0 ; i <= d ; i++ ) zp_neg(h[i],f[i],p); }

// h = f1*f2.  Assumes either d1 or d2 (or both) is small.  Same functionality as zp_poly_mul
// w is a temporary provide by the caller.  up to 2*log2(p) + log2(min(d1,d2)) bits may be written to w,
// but the entries of h will be smaller than p at all times
void zp_poly_mul_mod_xn_small (zp_t h[], zp_t  f1[], int d1, zp_t  f2[], int d2, int n, zp_p p);

// h = f1*f2 mod x^n
// the output h will have n coefficients written to it (possibly padded with zeros)
static inline void zp_poly_mul_mod_xn (zp_t h[], zp_t f1[], int d1, zp_t f2[], int d2, int n, zp_p p)
{
	if ( d1 <= ZP_POLY_VERY_SMALL_DEGREE || d2 <= ZP_POLY_VERY_SMALL_DEGREE ) {
		zp_poly_mul_mod_xn_small (h, f1, d1, f2, d2, n, p);
	} else {
		zp_array_mul_mod_xn (h, f1, d1+1, f2, d2+1, 0, n, p);
	}
}

// h = f1*f2.   (any or all of f1,f2,h may coincide)  dh = d1 + d2 (or -1 if d1 or d2=-1)
// the case f1=f2 will be detected and optimized.  returns degree of product
static inline int zp_poly_mul (zp_t h[], zp_t f1[], int d1, zp_t f2[], int d2, zp_p p)
	{ if ( d1 < 0 || d2 < 0 ) return -1; zp_poly_mul_mod_xn (h, f1, d1, f2, d2, d1+d2+1, p); return d1+d2; }

//  h = f1*f2, where f1 and f2 are monic polynomials.  The leading monic coefficient of h is not written, so only d1+d2 coefficients are needed
// h cannot overlap f1 or f2, but f1 and f2 may coincide (and this case is handled efficiently)
static inline void zp_poly_mul_monic (zp_t h[], zp_t f1[], int d1, zp_t f2[], int d2, zp_p p)
{
	register int i;
	
	if ( d1 < 0 || d2 < 0 ) return;
	if ( ! d1 ) { zp_poly_copy (h, f2, d2-1); return; }
	if ( ! d2 ) { zp_poly_copy (h, f1, d1-1); return; }
	zp_poly_mul (h, f1, d1-1, f2, d2-1, p);
	zp_add (h[d1+d2-1], f1[d1-1], f2[d2-1], p);
	for ( i = d1-2 ; i >= 0 ; i-- ) zp_add (h[d2+i], h[d2+i], f1[i], p);
	for ( i = d2-2 ; i >= 0 ; i-- ) zp_add (h[d1+i], h[d1+i], f2[i], p);
	return;
}

// h += f1*f2 with at least one of d1 or d2 very small
void zp_poly_addmul_small (zp_t h[], int *dh, zp_t f1[], int d1, zp_t f2[], int d2, zp_p p);

// h -= f1*f2 with at least one of d1 or d2 very small
void zp_poly_submul_small (zp_t h[], int *dh, zp_t f1[], int d1, zp_t f2[], int d2, zp_p p);

// h = h + f1 * f2
static inline void zp_poly_addmul (zp_t h[], int *dh, zp_t f1[], int d1, zp_t f2[], int d2, zp_p p)
{
	int n;
	
	if ( d1 <= ZP_POLY_VERY_SMALL_DEGREE || d2 <= ZP_POLY_VERY_SMALL_DEGREE ) {
		zp_poly_addmul_small (h, dh, f1, d1, f2, d2, p);
	} else {
		n = zp_imax(*dh,d1+d2);
		zp_array_mul_mod_xn (h, f1, d1+1, f2, d2+1, *dh+1, n+1, p);
		*dh = zp_poly_degree(h, n);
	}
}

// h = h - f1 * f2
static inline void zp_poly_submul (zp_t h[], int *dh, zp_t f1[], int d1, zp_t f2[], int d2, zp_p p)
{
	int n; 

	if ( d1 <= ZP_POLY_VERY_SMALL_DEGREE || d2 <= ZP_POLY_VERY_SMALL_DEGREE ) {
		zp_poly_submul_small (h, dh, f1, d1, f2, d2, p);
	} else {
		if ( *dh < 0 ) { if ( d1 < 0 || d2 < 0 ) return; mpz_set_ui(h[0],0); *dh = 0; }	  // hack to get a zero poly with non-zero length
		n = zp_imax(*dh,d1+d2);
		zp_array_mul_mod_xn (h, f1, d1+1, f2, d2+1, -((*dh)+1), n+1, p); 
		*dh = zp_poly_degree(h, n);
	}
}

// h = f1*f2 + f3*f4, h must be large enough to hold either product, even if the sum has lower degree than this
// use zp_poly_mul2_mod_xn when you know the result will be smaller.
static inline void zp_poly_mul2 (zp_t h[], int *dh, zp_t f1[], int d1, zp_t f2[], int d2, zp_t f3[], int d3, zp_t f4[], int d4, zp_p p)
 { int n = zp_imax(d1+d2,d3+d4)+1;  zp_array_mul2_mod_xn (h, f1, d1+1, f2, d2+1, f3, d3+1, f4, d4+1, n, p);  *dh = zp_poly_degree (h, n-1); }
 
// h = f1*f2 + f3*f4 mod x^n
static inline void zp_poly_mul2_mod_xn (zp_t h[], zp_t f1[], int d1, zp_t f2[], int d2,
								 zp_t f3[], int d3, zp_t f4[], int d4, int n, zp_p p)
 { zp_array_mul2_mod_xn (h, f1, d1+1, f2, d2+1, f3, d3+1, f4, d4+1, n, p); }

// h = f1*f2 + f3*f4 where n is known to be an upper bound on the degree of h (results undefined otherwise)
static inline void zp_poly_mul42n (zp_t c0[], zp_t c1[], zp_t a00[], int d00, zp_t a01[], int d01, zp_t a10[], int d10, zp_t a11[], int d11,
                                                                   zp_t b0[], int e0, zp_t b1[], int e1, int n, zp_p p)
 { zp_array_mul42n (c0, c1, a00, d00+1, a01, d01+1, a10, d10+1, a11, d11+1, b0, e0+1, b1, e1+1, n, p); }
	 
// c0 = a00*b0 + a01*b1 mod x^n and c1 = a10*b0 + a11*b1 mod x^n (full aliasing supported)
static inline void zp_poly_mul42_mod_xn (zp_t c0[], zp_t c1[], zp_t a00[], int d00, zp_t a01[], int d01, zp_t a10[], int d10, zp_t a11[], int d11,
                                                                   zp_t b0[], int e0, zp_t b1[], int e1, int n, zp_p p)
 { zp_array_mul42_mod_xn (c0, c1, a00, d00+1, a01, d01+1, a10, d10+1, a11, d11+1, b0, e0+1, b1, e1+1, n, p); }

// c = a * b where a, b, and c are 2x2 matrices (c may alias a or b or both)
static inline void zp_poly_mul4_mod_xn (zp_t c00[], zp_t c01[], zp_t c10[], zp_t c11[],
	                                                         zp_t a00[], int d00, zp_t a01[], int d01, zp_t a10[], int d10, zp_t a11[], int d11,
	                                                         zp_t b00[], int e00, zp_t b01[], int e01, zp_t b10[], int e10, zp_t b11[], int e11, int n, zp_p p)
  { zp_array_mul4_mod_xn (c00, c01, c10, c11, a00, d00+1, a01, d01+1, a10, d10+1, a11, d11+1, b00, e00+1, b01, e01+1, b10, e10+1, b11, e11+1, n, p); }
 
// q = a div b with b not nesc. monic (faster if it is).  Assumes da-db is very small.  q may alias a but not b
void zp_poly_div_small  (zp_t q[], zp_t a[], int da, zp_t b[], int db, zp_p p);

// r = a - qb where q is the Euclidean quotient of a/b and da >= db > dr, r is zero filled to db-1, degree of r returned
int zp_poly_rem_small (zp_t r[], zp_t a[], int da, zp_t b[], int db, zp_t q[], zp_p p);

// compute h = 1/f mod x^n, given f(0)=1
void zp_poly_inverse_mod_xn (zp_t h[], zp_t f[], int d, int n, zp_p p);

// compute the Euclidean quotient q = x^n / g for g monic
// if n < d then the quotient is zero and this function does nothing
// q and g may not overlap
void zp_poly_div_xn (zp_t q[], zp_t g[], int d, int n, zp_p p);

// compute the Euclidea quotient q = a / b (b need not be monic)
// q may equal a but not b.  Use zp_poly_rem to get the remainder
int zp_poly_div (zp_t q[], zp_t a[], int da, zp_t b[], int db, zp_p p);

// given Euclidean qoutient q = a / b, computes r = a - qb with deg r < deg b, r is zero filled to db-1, degree of r returned
// this function does not require b monic
int zp_poly_rem (zp_t r[], zp_t a[], int da, zp_t b[], int db, zp_t q[], zp_p p);

// computes r = a mod b by calling zp_poly_div and zp_poly_rem
int zp_poly_mod (zp_t r[], zp_t a[], int da, zp_t b[], int db, zp_p p);

// r = radix rep of a in base b, that is a = sum r_i b^i
// where each r_i has deg < db and there are ceil(da/db) terms
// the output is the concatenation of da/db+1 polys of deg db-1 (each zero padded)
// with total size db * (da/db+1) < da + db
void zp_poly_radix (zp_t r[], zp_t a[], int da, zp_t b[], int db, zp_p p);

// compute r, s, t so that r = s*f + t*g, r is a monic divisor of f and g,
// ds < max(dg,1), dt < max(df,1).  (f and g need not be monic)
// both s and t are optional and may be null if not required
void zp_poly_xgcd (zp_t r[], int *dr, zp_t s[], int *ds, zp_t t[], int *dt,
			      zp_t f[], int df, zp_t g[], int dg, zp_p p);

// r = gcd(f,g) normalized so that r is monic (f and g need not be monic)
void zp_poly_gcd (zp_t r[], int *dr, zp_t f[], int df,zp_t g[], int dg, zp_p p);

// h = 1/f mod g, where f and h are zero-filled to degree d-1
// returns 0 if no inverse of f exists, 1 otherwise (g need not be monic)
// see below for zp_poly_mod interface to the same function
int zp_poly_inv_mod (zp_t h[], zp_t f[], zp_t g[], int d, zp_p p);

// r = f(g)
void zp_poly_compose_horner (zp_t r[], zp_t f[], int df, zp_t g[], int dg, zp_p p);
static inline void zp_poly_compose (zp_t r[], zp_t f[], int df, zp_t g[], int dg, zp_p p)
	{ zp_poly_compose_horner (r, f, df, g, dg, p); }

// possible values of algo parameter to zp_poly_mod_init
#define ZP_POLY_MOD_ALGO_AUTO		       	0  // select automatically
#define ZP_POLY_MOD_ALGO_STUPID 		1  // quadratic-time
#define ZP_POLY_MOD_ALGO_BARRETT    		2  // variant of Barrett's algorithm
#define ZP_POLY_MOD_ALGO_MONTGOMERY 	3  // variant of Montgomery's algorithm -- not used, requires import/export
#define ZP_POLY_MOD_ALGO_KRONECKER	4  // naive Kronecker substitution

typedef zp_poly_mod_struct zp_poly_mod_t[1];	// defined in zp_poly_mul.h

// initilizes mod structure, input poly g must be monic
void zp_poly_mod_init (zp_poly_mod_t mod, zp_t g[], int d, zp_p p, int algo);
void zp_poly_mod_clear (zp_poly_mod_t mod);

/*
	All arrays passed in and out of zp_poly_mod functions have length deg(g)
	(with leading zeros as required)
*/

static inline zp_t *zp_poly_mod_alloc (zp_poly_mod_t g)
	{ return zp_poly_alloc (g->d-1, g->p); }

static inline void zp_poly_mod_free (zp_t *f, zp_poly_mod_t g)
	{ zp_poly_free (f, g->d-1); }
	
// h = f
static inline void zp_poly_mod_copy (zp_t h[], zp_t f[], zp_poly_mod_t g)
{
	register int i;
	if ( f != h ) for ( i = 0 ; i < g->d ; i++ ) mpz_set(h[i],f[i]);
}

// ( f1 == f2 ? 1 : 0 )
static inline int zp_poly_mod_equal (zp_t f1[], zp_t f2[], zp_poly_mod_t g)
{
	
	register int i;
	for ( i = 0 ; i < g->d ; i++ ) if ( ! zp_equal (f1[i],f2[i]) ) return 0;
	return 1;
}
	
// h = x^e * f1 * f2 mod g, where e = 0 or e = 1
void zp_poly_mod_mulx (zp_t h[], zp_t f1[], zp_t f2[], int e, zp_poly_mod_t g);

// h = f1 * f2 mod g
static inline void zp_poly_mod_mul (zp_t h[], zp_t f1[], zp_t f2[], zp_poly_mod_t g)
	{ zp_poly_mod_mulx (h, f1, f2, 0, g); }

// h = 1 / f mod g
static inline int zp_poly_mod_inv (zp_t h[], zp_t f[], zp_poly_mod_t g)
	{ return zp_poly_inv_mod (h, f, g->g, g->d, g->p); }

// h[i] = 1/ f[i] mod g for i in [0,n-1]
int zp_poly_mod_array_inv (zp_t *h[], zp_t *f[], int n, zp_poly_mod_t g);

// ( f == 0 ? 1 : 0 )
static inline int zp_poly_mod_is_zero (zp_t f[], zp_poly_mod_t g)
	{ return zp_poly_is_zero (f, g->d-1); }
	
// h = 0
static inline void zp_poly_mod_set_zero (zp_t h[], zp_poly_mod_t g)
	{ zp_poly_zero_fill (h, g->d-1); }
	
// ( f == 1 ? 1 : 0 )
static inline int zp_poly_mod_is_one (zp_t f[], zp_poly_mod_t g)
	{ return zp_poly_is_one (f, g->d-1); }

// h = 1
static inline void zp_poly_mod_set_one (zp_t h[], zp_poly_mod_t g)
	{ zp_set_one (h[0]);  zp_poly_zero_pad (h, 1, g->d-1); }

// h = random poly mod g
static inline void zp_poly_mod_randomize (zp_t h[], zp_poly_mod_t g)
	{ zp_poly_randomize (h, g->d-1, g->p); }	
	
// h = f1 + f2
static inline void zp_poly_mod_add (zp_t h[], zp_t f1[], zp_t f2[], zp_poly_mod_t g)
{
	register int i;
	for ( i = 0 ; i < g->d ; i++ ) zp_add (h[i], f1[i], f2[i], g->p);
}

// h = 3 * f
static inline void zp_poly_mod_triple (zp_t h[], zp_t f[], zp_poly_mod_t g)
	{ zp_poly_triple (h, f, g->d-1, g->p); }

// h = f1 - f2
static inline void zp_poly_mod_sub (zp_t h[], zp_t f1[], zp_t f2[], zp_poly_mod_t g)
{
	register int i;
	for ( i = 0 ; i < g->d ; i++ ) zp_sub (h[i], f1[i], f2[i], g->p);
}

// h = -f
static inline void zp_poly_mod_neg (zp_t h[], zp_t f[], zp_poly_mod_t g)
{
	register int i;
	for ( i = 0 ; i < g->d ; i++ ) zp_neg (h[i], f[i], g->p);
}

// h = c * h
static inline void zp_poly_mod_smul (zp_t h[], zp_t f[], zp_t c, zp_poly_mod_t g)
	{ zp_poly_smul (h, f, g->d-1, c, g->p); }
	
// h = f^n mod g (where f and h are both zero padded to degree deg(g)-1)
void zp_poly_mod_pow (zp_t h[], zp_t f[], zp_t n, zp_poly_mod_t g);

// h = x^n mod g (zero padded to degree deg(g) -1)
void zp_poly_mod_pow_xn (zp_t h[], zp_t n, zp_poly_mod_t g);

// h = (x+a)^n mod g (zero padded to degree deg(g) -1)
void zp_poly_mod_pow_xan (zp_t h[], zp_t a, zp_t n, zp_poly_mod_t g);

// h = (x^3 + a*x + b) * f mod g (zero padded to deg(g)-1)
// h may coincide with f
void zp_poly_mod_mul_by_x3axb (zp_t h[], zp_t f[], zp_t a, zp_t b, zp_poly_mod_t g);

// computes h = (x^3+ax+b)^n mod g (zero padded to degree deg(g) -1)
void zp_poly_mod_pow_x3axbn (zp_t h[], zp_t a, zp_t b, zp_t n, zp_poly_mod_t mod);

// finds a root of a monic poly f if one exists.  returns total number of distinct roots
int zp_poly_find_root (zp_t r, zp_t f[], int d, zp_mod_t p);

// finds up to k distinct roots of a monic poly f 
// returns the total number of distinct roots (possibly greater than k, but only k roots will be placed in r[])
int zp_poly_find_roots (zp_t r[], int k, zp_t f[], int d, zp_mod_t p);
	
// finds a root of a polynomial that is assumed to be a product of linear factors, destroying it in the process
void zp_poly_find_split_root (zp_t r, zp_t h[], int d, zp_mod_t p);
void zp_poly_find_split_roots (zp_t r[], int k, zp_t h[], int d, zp_mod_t p);

// factors monic poly f into distinct irreducible monic polys of degree n[i] and multiplicity e[i], returns the number of factors
// the distinct factors of f are stored in the array r as a concatenated list of implicitly monic polys with leading coefficient omitted 
// so r should have space for d coefficients (if f is square-free, this is exactly how many coefficients will be present)
int zp_poly_factor (zp_t r[], int n[], int e[], zp_t f[], int d, zp_mod_t p);

// Finds a factor a monic poly f of degree d that has minimal positive degree.  Returns the degree
int zp_poly_min_degree_factor (zp_t r[], zp_t f[], int d, zp_mod_t p);

// g = f * (x - r) where f is a poly in Z/pZ[x] and r is an element of Z/pZ
void zp_poly_add_root (zp_t g[], zp_t f[], int d, zp_t r, zp_p p);

// computes g=f/(x-r) for f in Z/pZ[x] assuming f(r) = 0.  f and g may point to the same place
void zp_poly_remove_root (zp_t g[], zp_t f[], int d, zp_t r, zp_p p);

// replace monic f(x) by f(x-f[d-1]/d), killing x^(d-1) coefficient
// sets s to f[d-1]/d, if r is a root of the new f then r-s is a root of the old f
void zp_poly_depress_monic_inplace (zp_t f[], int d, zp_t s, zp_p p);

// f = \prod_i  (x-r[i]), f and r may point to the same location
static inline void zp_poly_from_roots_naive (zp_t f[], zp_t r[], int d, zp_p p)
{
	zp_t t,ri;
	register int i,j;

	mpz_init_modp (t, p);  mpz_init_modp (ri, p);
	mpz_set_one (f[d]);
	mpz_neg_mod (f[d-1], r[d-1], p);
	for ( i = d-2 ; i >= 0 ; i-- ) {
		mpz_set (ri, r[i]);
		mpz_mul_mod (t,ri,f[i+1], p);
		mpz_neg_mod (f[i],t,p);
		for ( j = i+1 ; j < d-1 ; j++ ) {
			mpz_mul_mod (t,ri,f[j+1],p);
			mpz_sub_mod (f[j],f[j],t,p);
		}
		mpz_sub_mod (f[j],f[j],ri,p);
	}
	mpz_clear (t);  mpz_clear (ri);
}

// f = \prod_i  (x-r[i]), f and r may point to the same location
void zp_poly_from_roots (zp_t f[], zp_t r[], int d, zp_p p);

#define ZP_POLY_MOD_COMP_WINDOW_SIZE		16

// r = f(g) mod h
void zp_poly_mod_compose_window (zp_t r[], zp_t f[], int d, zp_t g[], zp_poly_mod_t h, int k);
static inline void zp_poly_mod_compose (zp_t r[], zp_t f[], int d, zp_t g[], zp_poly_mod_t h)
{
	int k;

	// window size should never be bigger than sqrt(d)
	k = ZP_POLY_MOD_COMP_WINDOW_SIZE;
	while ( k*k > d ) k--;
	zp_poly_mod_compose_window (r, f, d, g, h, k);
}

int zp_poly_jacobi_symbol (zp_t u[], int du, zp_t v[], int dv, zp_p p);

 #endif
