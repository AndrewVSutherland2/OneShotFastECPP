/*

zp_poly_util.h: miscellaneous utility/support functions for zp_poly

Copyright (C) 2010, David Harvey and Andrew V. Sutherland

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

#ifndef ZP_POLY_UTIL_H
#define ZP_POLY_UTIL_H

#include <assert.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>

#define MPZ_EXP2_MINPBITS		3000

#define delta_secs(s,e)	((double)(e-s)/(CLOCKS_PER_SEC))

// TODO: make crossover points depend on the size of p
// For log2(p) > 7500, polys of deg > 4 are already better handled
// by zp_array_mul than schoolbook.
#define ZP_POLY_VERY_SMALL_DEGREE		1

// we try to avoid gmp reallocating mpz_t's by initializing mod p values
// to pbits plus a fudge factor (this is bigger than it should need to be,
// but GMP is not always so reliable on this point)
#define ZP_PBITS_FUDGE		48

// prime modulus structure -- not required by all functions but needed for zp2 functions that work Fp^2, and also to compute square and cube roots in Fp.
typedef struct zp_mod_struct {
	mpz_t p;			// a prime > 3 (FOR BACKWARD COMPATIBILITY THIS MUST BE THE FIRST ELEMENT)
	mpz_t pm1;		// p-1
	mpz_t o;			// odd part of p-1
	int k;			// 2-adic valuation of p-1, so p-1=2^k*o
	short p1mod3;	// (p % 3) == 1 ? 1 : 0
	short p1mod4;	// (p % 4) == 1 ? 1 : 0 
	unsigned long s;	// least positive non-residue mod p
	mpz_t si;			// inverse of s mod p
	mpz_t t1;			// work variable for zp_poly_util
	mpz_t w[2], t[2];	// handy work-variables for zp2 (THESE ARE RESERVED FOR THE SOLE USE OF INLINES IN ZP2.H - ADD MORE WORK VARIABLES IF NEEDED)
	// items below are used in zp_roots.h
	mpz_t g2;		// generator of the 2-Sylow subgoup of Fp
	mpz_t g3[2];		// generator of the 3-Sylow subgroup of Fp^2
	mpz_t z3[2];		// cube-root of unity in Fp^2
} zp_mod_t[1];

typedef mpz_t zp_t;				// zp_t type indicates an mpz_t in the range [0,p-1]
typedef void *zp_p;				// Functions that require the modulus p as an mpz_t by not necessarily a zp_mod_t are declared using zp_p so that either an mpz_t or zp_mod_t may be specified
typedef __mpz_struct *zp_p_to_mpz;	// used to cast a zp_p to an mpz_t


void zp_mod_init (zp_mod_t p, mpz_t P);
void zp_mod_clear (zp_mod_t p);
void zp_randomb (mpz_t a, int b);
int zp_random_bit ();
void zp_randomm (mpz_t a, mpz_t m);
static inline void zp_random (zp_t a, zp_p p) { zp_randomm (a, (zp_p_to_mpz)p); }
void zp_cleanup (void) ;		// frees rand state

// try to ensure our integer min/max functions don't collide with anything
static inline int zp_imax (int a, int b) { return ( a > b ? a : b ); }
static inline int zp_imin (int a, int b) { return ( a < b ? a : b ); }

// ui_bits returns the least i for which 2^i > n (# bits in binary rep of n)
// this function could be made faster within an inline assembly directive
static inline int ui_bits (unsigned long n)
	{ register int i;  for ( i = 0 ; (1UL<<i) <= n ; i++ );  return i; }
static inline long mpz_bits (zp_t n) { return mpz_sizeinbase (n, 2); }
static inline long zp_bits (zp_p p) { return mpz_sizeinbase ((zp_p_to_mpz)p, 2); }

#ifndef ZP_POLY_UTIL_C
extern long const zp_mem_used, zp_mem_peak, zp_malloc_count, zp_realloc_count, zp_free_count;
#endif

// malloc/realloc/free memory wrappers to track memory usage
void *zp_malloc (size_t n);
void *zp_realloc (void *p, size_t n1, size_t n2);
void zp_free (void *p, size_t n);

// outputs memory stats to stdout
// by default (verbosity=0) no output when mem_used = 0
void zp_mem_report (int verbosity);

// allocates and inits an array of n mpz_t's with specified initial allocation size
static inline mpz_t *zp_array_alloc (int n, int bits)
{
	mpz_t *v;
	int i;
	
	v = (mpz_t *) zp_malloc (n * sizeof(mpz_t));
	for ( i = 0 ; i < n ; i++ ) mpz_init2 (v[i], bits);
	return v;
}

static inline void zp_array_free (mpz_t v[], int n)
{
	int i;
	
	for ( i = 0 ; i < n ; i++ ) mpz_clear (v[i]);
	zp_free (v, n * sizeof (mpz_t));
}

// allocate and init d+1 coeffs for degree d poly (coeffs will be zero)
static inline mpz_t *zp_poly_alloc (int d, zp_p p)
	{ return zp_array_alloc (d+1, zp_bits(p)+ZP_PBITS_FUDGE); }
	
// free and clear the d+1 coeffs of a degree d poly
static inline void zp_poly_free (zp_t f[], int d)
	{ zp_array_free (f, d+1); }

// stack for more efficient memory allocation
typedef struct _zp_poly_stack_struct {
	mpz_t *beg, *end, *next, *init, *high;
	int bits, size;
} zp_poly_stack_t[1];

void zp_poly_stack_init (zp_poly_stack_t stack, int size, mpz_t p);
void zp_poly_stack_clear (zp_poly_stack_t stack);

// allocates a poly of degree d (with d+1 coeffs) off the stack (contents arbitrary!)
static inline mpz_t *zp_poly_stack_alloc (int d, zp_poly_stack_t stack)
{
	mpz_t *f, *z;

	f = stack->next;
//printf ("Request to allocate %d coefficients at %lx\n", d+1, f);
	stack->next += d+1;
	if ( stack->high < stack->next ) stack->high = stack->next;
	assert (stack->next <= stack->end);
	for ( z = stack->init ; z < stack->next ; z++ ) mpz_init2(*z,stack->bits+ZP_PBITS_FUDGE);
	stack->init = z;
//printf ("stack alloced %d coefficients at %lx\n", d+1, f);
	return f;
}

static inline mpz_t *zp_poly_stack_mark (zp_poly_stack_t stack)
	{ return stack->next; }
	
static inline void zp_poly_stack_pop (zp_poly_stack_t stack, mpz_t *mark)
{
	assert ( mark >= stack->beg && mark <= stack->next );
	stack->next = mark;
//printf ("stack popped %d entries to %lx\n", stack->next - mark, mark);
}

// reduces b mod p assuming b lies in [1-p,2p-1]
static inline void mpz_mod_small (mpz_t a, mpz_t b, mpz_t p)
{
	if ( mpz_sgn(b) < 0 ) mpz_add(a,b,p);
	else if ( mpz_cmp(b,p) >= 0 ) mpz_sub(a,b,p);
	else mpz_set(a,b);
}


// for backward compatability we duplicate many of the functions below
// note that mpz functions take p as an mpz_t whereas zp functions take p as a zp_mod_t

// init for standard values mod p, not used to store products
static inline void mpz_init_modp (mpz_t c, mpz_t p)
	{ mpz_init2(c,zp_bits(p)+ZP_PBITS_FUDGE); }
static inline void zp_init (zp_t c, zp_p p) { mpz_init_modp (c, (zp_p_to_mpz)p); }

// init for working variables that may need to store products
static inline void mpz_init_modp_w (mpz_t w, mpz_t p)
	{ mpz_init2(w,2*zp_bits(p)+ZP_PBITS_FUDGE); }
static inline void zp_init_w (mpz_t w, zp_p p) { mpz_init_modp (w, (zp_p_to_mpz)p); }

static inline void zp_clear (zp_t c) { mpz_clear (c); }
	
// sets a = -b mod p, assuming b lies in [0,p-1]
static inline void mpz_neg_mod(mpz_t a, mpz_t b, mpz_t p)
	{ mpz_neg(a,b); if ( mpz_sgn(a) < 0 ) mpz_add(a,a,p); }
static inline void zp_neg (zp_t a, zp_t b, zp_p p) { mpz_neg_mod (a, b, (zp_p_to_mpz)p); }
	
// c = a+b mod p, assuming a,b in [0,p-1]
static inline void mpz_add_mod (mpz_t c, mpz_t a, mpz_t b, mpz_t p)
	{ mpz_add(c,a,b); if ( mpz_cmp(c,p) >= 0 ) mpz_sub(c,c,p); }
static inline void zp_add (zp_t c, zp_t a, zp_t b, zp_p p) { mpz_add_mod (c, a, b, (zp_p_to_mpz)p); }
static inline void zp_add_ui (zp_t c, zp_t a, unsigned long b, zp_p p)
	{ mpz_add_ui(c,a,b); mpz_mod (c, c, (zp_p_to_mpz)p); }

// c = a-b mod p, assuming a,b in [0,p-1]
static inline void mpz_sub_mod (mpz_t c, mpz_t a, mpz_t b, mpz_t p)
	{ mpz_sub(c,a,b); if ( mpz_sgn(c) < 0 ) mpz_add(c,c,p); }
static inline void zp_sub  (zp_t c, zp_t a, zp_t b, zp_p p) { mpz_sub_mod (c, a, b, (zp_p_to_mpz)p); }
static inline void zp_sub_ui   (zp_t c, zp_t a, unsigned long b, zp_p p)
	{ mpz_sub_ui(c,a,b); mpz_mod (c, c, (zp_p_to_mpz)p); }

// a++ mod p, assuming a,b in [0,p-1]
static inline void mpz_inc_mod (mpz_t a, mpz_t p)
	{ mpz_add_ui(a,a,1); if ( mpz_cmp(a,p) >= 0 ) mpz_sub (a,a,p); }
static inline void zp_inc  (zp_t a, zp_p p) { mpz_inc_mod (a, (zp_p_to_mpz)p); }

// a-- mod p, assuming a,b in [0,p-1]
static inline void mpz_dec_mod (mpz_t a, mpz_t p)
	{ mpz_sub_ui(a,a,1); if ( mpz_sgn(a) < 0 ) mpz_add (a,a,p); }
static inline void zp_dec (zp_t a, zp_p p) { mpz_dec_mod (a, (zp_p_to_mpz)p); }
	
// c = a*b mod p (doesn't require a or b reduced, may expand c)
static inline void mpz_mul_mod (mpz_t c, mpz_t a, mpz_t b, mpz_t p)
	{ mpz_mul (c, a, b); mpz_mod (c, c, p); }
static inline void zp_mul  (zp_t c, zp_t a, zp_t b, zp_p p) { mpz_mul_mod (c, a, b, (zp_p_to_mpz)p); }
static inline void zp_mul_ui  (zp_t c, zp_t a, unsigned long b, zp_p p)
	{ mpz_mul_ui (c, a, b); mpz_mod (c, c, (zp_p_to_mpz)p); }
// b = a^2 mod p (doesn't require a reduced, may expand b)
static inline void zp_sqr (zp_t b, zp_t a, zp_p p) { zp_mul (b, a, a, p); }
static inline void zp_cube (zp_t b, zp_t a, zp_p p)
{
	if ( b == a ) {zp_t w;  zp_init_w (w,p); zp_sqr (w,a,p); zp_mul (b,w,a,p); zp_clear (w); }
	else { zp_sqr (b, a, p); zp_mul (b, b, a, p); }
}

// b = a^e mod p
static inline void  zp_exp (zp_t b, zp_t a, mpz_t e, zp_p p)
	{ if ( mpz_cmp (e, (zp_p_to_mpz)p) >= 0 ) { zp_t w; mpz_sub_ui(w,(zp_p_to_mpz)p,1); mpz_mod (w, e, w); mpz_powm (b, a, w, (zp_p_to_mpz)p); zp_clear(w); } else mpz_powm (b, a, e, (zp_p_to_mpz)p); }
static inline void  zp_exp_ui (zp_t b, zp_t a, unsigned long e, zp_p p)
	{ mpz_powm_ui (b, a, e, (zp_p_to_mpz)p); }

// returns the Jacobi symbo (a/p)
static inline int zp_jacobi_symbol (zp_t a, zp_p p)
	{ return mpz_jacobi (a, (zp_p_to_mpz)p); }
	
// b1 = a^e1, b2 = a^e2 (for modest p it is faster to just call powm twice, but for large p it isn't)
void _zp_exp2 (mpz_t b1, mpz_t b2, mpz_t a, mpz_t e1, mpz_t e2, zp_mod_t p);
	
// a may coincide with b1 or b2 (but not both!)
static inline void zp_exp2 (mpz_t b1, mpz_t b2, mpz_t a, mpz_t e1, mpz_t e2, zp_mod_t p)
{
	if ( zp_bits(p) < MPZ_EXP2_MINPBITS ) {
		if ( a == b1 ) { zp_exp (b2, a, e2, p); zp_exp (b1, a, e1, p); }
		else { zp_exp (b1, a, e1, p); zp_exp (b2, a, e2, p); }
	} else {
		_zp_exp2 (b1, b2, a, e1, e2, p);
	}
}

// b = a/2 mod p
static inline void mpz_div2_mod (mpz_t b, mpz_t a, mpz_t p)
{
	assert (mpz_tstbit(p,0));	// p must be odd
	if ( mpz_tstbit(a,0) ) {
		mpz_add (b,a,p);
		mpz_div_2exp (b,b,1);
	} else {
		mpz_div_2exp (b,a,1);
	}
}static inline void zp_div2 (zp_t b, zp_t a, zp_p p) { mpz_div2_mod (b, a, (zp_p_to_mpz)p); }

// b = a/3 mod p (output is reduced if  a is but otherwise not necessarily)
static inline void mpz_div3_mod (mpz_t b, mpz_t a, mpz_t p)
{
	assert (mpz_cmp_ui(p,3) > 0 );	// p must be > 3
	mpz_set (b, a);
	while ( ! mpz_divisible_ui_p (b, 3) ) mpz_add (b, b, p);
	mpz_divexact_ui (b, b, 3);
}
static inline void zp_div3 (zp_t b, zp_t a, zp_p p) { mpz_div3_mod (b, a,(zp_p_to_mpz)p); }

// b = 1/a mod p
static inline int zp_invert (zp_t b, zp_t a, zp_p p) { return mpz_invert (b, a, (zp_p_to_mpz)p); }
	
// o[i] = 1 / a[i] mod p for i in [0,n-1].  returns 0 if any a[i] is 0 mod p, 1 ow.
int mpz_array_invert (mpz_t o[], mpz_t a[], int n, mpz_t p);
static inline int zp_array_invert (mpz_t o[], mpz_t a[], int n, zp_p p) { return mpz_array_invert (o, a, n,(zp_p_to_mpz)p); }

// handy functions for readability, could be made slightly faster by going underneath
static inline void zp_set (zp_t b, zp_t a) { mpz_set (b, a); }
static inline void zp_set_ui (zp_t b, unsigned long a, zp_p p) { mpz_set_ui (b, a); mpz_mod (b, b, (zp_p_to_mpz)p); }
static inline void zp_set_mpz (zp_t b, mpz_t a, zp_p p) { mpz_mod (b, a, (zp_p_to_mpz)p); }
static inline void zp_set_str (zp_t b, char *str, zp_p p) { mpz_set_str (b, str, 0); mpz_mod (b, b, (zp_p_to_mpz)p); }
static inline void zp_set_zero (zp_t a) { mpz_set_ui(a,0); }
static inline void mpz_set_zero (mpz_t a) { mpz_set_ui(a,0); }
static inline void zp_set_one (zp_t a) { mpz_set_ui(a,1); }
static inline void mpz_set_one (mpz_t a) { mpz_set_ui(a,1); }
static inline int zp_is_one (zp_t a) { return !mpz_cmp_ui(a,1); }
static inline int mpz_is_one(mpz_t a) { return !mpz_cmp_ui(a,1); }
static inline int zp_is_zero (zp_t a) { return !mpz_sgn(a); }
static inline int mpz_is_zero (mpz_t a) { return !mpz_sgn(a); }
static inline int zp_equal (zp_t a, zp_t b) { return mpz_cmp(a,b) == 0; }

void zp_random_prime (mpz_t p, int b);

static inline void zp_poly_randomize (zp_t f[], int d, zp_p p)
{
	register int i;
	
	if ( d < 0 ) return;
	for ( i = 0 ; i <= d ; i++ ) zp_randomm (f[i], p);
	// make lc nonzero
	while ( zp_is_zero(f[d]) ) zp_randomm (f[d], p);
}

// Updates h to the next poly of degree d in the lexicographic ordering, returns 0 if h is maximal poly of degree d
static inline int zp_poly_next (zp_t h[], int d, zp_p p)
{
	register int i;
	
	for ( i = 0 ; i <= d ; i++ ) {
		mpz_add_ui(h[i],h[i],1);
		if ( mpz_cmp(h[i],(zp_p_to_mpz)p) ) break;
		else mpz_set_ui(h[i],0);
	}
	return ( i <= d ? 1 : 0 );
}

static inline void zp_poly_lift_to_Z (mpz_t f[], zp_t g[], int d, zp_p p)
{
	register int i;
	mpz_t x;
	
	mpz_init_set (x, (zp_p_to_mpz)p);
	mpz_div_2exp (x, x, 1);
	for ( i = 0 ; i <= d ; i++ ) if ( mpz_cmp(g[i],x) > 0 ) { mpz_sub (f[i],p,g[i]); mpz_neg(f[i],f[i]); } else mpz_set (f[i],g[i]);
}


void zp_poly_print (zp_t f[], int d);

#endif
