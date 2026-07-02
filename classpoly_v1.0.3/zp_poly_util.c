/*

zp_poly_util.c: miscellaneous utility/support functions for zp_poly

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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>

// the following globals are exported in zp_poly_util.h
long zp_mem_used;
long zp_mem_peak;
long zp_malloc_count;
long zp_realloc_count;
long zp_free_count;

#define ZP_POLY_UTIL_C
#include "zp_poly_util.h"

gmp_randstate_t zp_randstate;
int zp_randstate_init;

static inline void zp_randinit (void) { if ( ! zp_randstate_init ) { gmp_randinit_default (zp_randstate); zp_randstate_init = 1; } }
void zp_randomb (mpz_t a, int b) { zp_randinit();  mpz_urandomb (a, zp_randstate, b); }
int zp_random_bit () { return gmp_urandomb_ui (zp_randstate, 1); }
void zp_randomm (mpz_t a, mpz_t p) { zp_randinit(); mpz_urandomm (a, zp_randstate, p); }
void zp_cleanup (void) { if ( zp_randstate_init) gmp_randclear (zp_randstate); }

void *zp_malloc (size_t n)
{
//printf ("zp_malloc %d\n", n);
	zp_mem_used += n;
	if ( zp_mem_used > zp_mem_peak ) zp_mem_peak = zp_mem_used;
	zp_malloc_count++;
	return malloc (n);
}

void *zp_realloc (void *p, size_t n1, size_t n2)
{
//printf ("realloc %d to %d\n", n1, n2);
	zp_mem_used += n2 - n1;
	if ( zp_mem_used > zp_mem_peak ) zp_mem_peak = zp_mem_used;
	zp_realloc_count++;
	return realloc (p, n2);
}

void zp_free (void *p, size_t n)
{
//printf ("zp_free %d\n", n);
	zp_mem_used -= n;
	zp_free_count++;
	free (p);
}

void zp_mem_report (int verbosity)
{
	if ( !verbosity && !zp_mem_used ) return;
	printf ("zp_mem report:\n");
	printf ("    %ld allocs, %ld reallocs, %ld frees\n",
		zp_malloc_count, zp_realloc_count, zp_free_count);
	printf ("    %ld bytes outstanding in %ld blocks\n", zp_mem_used,
		zp_malloc_count - zp_free_count);
	printf ("    %ld bytes peak usage\n", zp_mem_peak);
}

void zp_poly_stack_init (zp_poly_stack_t stack, int size, mpz_t p)
{
	stack->beg = (mpz_t *) zp_malloc (size*sizeof(mpz_t));
	stack->bits = zp_bits(p)+2;
	stack->size = size;
	stack->end = stack->beg+stack->size;
	stack->next = stack->high = stack->init = stack->beg;
//printf ("Created stack of size %d\n", stack->size);
}

void zp_poly_stack_clear (zp_poly_stack_t stack)
{
	register mpz_t *z;

//printf ("Used  %d of %d stack entries, freeing stack\n", stack->high-stack->beg, stack->size);
	for ( z = stack->beg ; z < stack->init ; z++ ) mpz_clear(*z);
	zp_free (stack->beg, stack->size*sizeof(mpz_t));
	// for safety's sake, clear these values
	stack->beg = stack->end = stack->next = 0;
	stack->size = 0;
}

// ptab[i] = least c s.t. 2^{1000*i}+c is a probable prime
int zp_poly_ptab[21] =
{
   1, 297, 841, 3993, 63, 2157, 9901, 5443, 3081, 9465,
   177, 10147, 13275, 16905, 34615, 14635, 19545, 24043, 12331, 1575,
   13071
};

void zp_random_prime (mpz_t p, int b)
{
	// hardwire some big known primes to avoid the slow call mpz_nextprime
	if ( !(b%1000)&& b <= 20000) {
		mpz_set_ui(p,1); mpz_mul_2exp(p,p,b); mpz_add_ui(p,p,zp_poly_ptab[b/1000]);
//		printf ("using p=2^%d+%d\n", 1000*(b/1000), zp_poly_ptab[b/1000]);
	} else if ( b == 1024 ) {
		mpz_ui_pow_ui (p, 2, 1023); mpz_add_ui(p,p,1155);
//		printf ("using p=2^1023+1155\n");
	} else if ( b == 2048 ) {
		mpz_ui_pow_ui (p, 2, 2047); mpz_add_ui(p,p,1919);
//		printf ("using p=2^2047+1919\n");
	} else if ( b == 4096 ) {
		mpz_ui_pow_ui (p, 2, 4095); mpz_add_ui(p,p,579);
//		printf ("using p=2^4095+579\n");
	} else if ( b == 8192 ) {
		mpz_ui_pow_ui (p, 2, 8191); mpz_add_ui(p,p,1911);
//		printf ("using p=2^8191+1911\n");
	} else if ( b == 9966 ) {
		mpz_ui_pow_ui (p,10,3000); mpz_add_ui(p,p,1027);
//		printf ("using p=10^3000+1027\n");
	} else if ( b == 16384 ) {
		mpz_ui_pow_ui (p, 2, 16383); mpz_add_ui(p,p,20253);
//		printf ("using p=2^16383+20253\n");
	} else if ( b == 16646 ) {
		mpz_ui_pow_ui (p,2,16612); mpz_mul_ui(p,p,16219299585); mpz_sub_ui(p,p,1);
//		printf ("using p=16219299585 * 2^16612 - 1\n");
	} else if ( b == 32768 ) {
		mpz_ui_pow_ui (p, 2, 32767); mpz_add_ui(p,p,34769);
//		printf ("using p=2^32767+34769\n");
	} else if ( b == 33549 ) {
		mpz_ui_pow_ui(p,2,33548); mpz_add_ui(p,p,4471);
//		printf ("using p=2^33548+4471\n");
	} else if ( b == 65536 ) {
		mpz_ui_pow_ui (p, 2, 65535); mpz_add_ui(p,p,37355);
//		printf ("using p=2^65535+37355\n");		
	} else {
		if ( b > 512 ) printf ("Generating %d-bit prime modulus...\n", b);
		do { zp_randomb (p,b); } while ( zp_bits(p) < b );
		mpz_nextprime (p,p);
		if ( mpz_cmp_ui(p,2) == 0 ) mpz_set_ui(p,3);	// don't use 2
//		if ( b > 512 ) gmp_printf ("using p=%Zd\n", p);
	}
}

void zp_poly_print (mpz_t f[], int d)
{
	register int i;

	if ( d < 0 ) { printf ("[0]"); return; }
	if ( d == 0 ) { gmp_printf ("[%Zd]", f[0]); return; }
	printf ("[");
	if ( ! mpz_is_one(f[d]) ) gmp_printf ("%Zd*", f[d]);
	if ( d > 0 ) printf ("x");
	if ( d >1 ) printf ("^%d", d);
	for ( i = d-1 ; i > 1 ; i-- ) gmp_printf (" + %Zd*x^%d", f[i], i);
	if ( i == 1 ) gmp_printf (" + %Zd*x", f[i--]);
	if ( i == 0 ) gmp_printf (" + %Zd", f[0]);
	printf ("]");
}

void zp_mod_init (zp_mod_t p, mpz_t P)
{
//	assert ( mpz_sizeinbase(P,2) >= 32 );

	mpz_init_set (p->p, P);  mpz_init_modp (p->pm1, P); mpz_init_modp (p->o, P);  mpz_init_modp (p->g2, P); mpz_init_modp (p->si, P); mpz_init_modp (p->t1, P);
	mpz_init_modp_w (p->w[0], P), mpz_init_modp_w (p->w[1], P);  mpz_init_modp_w (p->t[0], P), mpz_init_modp_w (p->t[1], P);
	mpz_init_modp (p->g3[0], P); mpz_init_modp (p->g3[1], P); mpz_init_modp (p->z3[0], P); mpz_init_modp (p->z3[1], P);

	mpz_sub_ui (p->pm1, P, 1);  mpz_set (p->o, p->pm1);
	for ( p->k = 1 ; !mpz_tstbit (p->o, p->k) ; p->k++ );
	mpz_div_2exp (p->o, p->o, p->k);
	p->p1mod3 = ( mpz_congruent_ui_p (P, 1, 3) ? 1 : 0 );
	p->p1mod4 = ( mpz_congruent_ui_p (P, 1, 4) ? 1 : 0 );
	// find the least positive non-residue mod P, which we assume fits in an unsigned long (assuming the GRH,  this will be true for any P that is likely to fit in memory)
	if ( mpz_ui_kronecker (2,P) < 0 ) p->s = 2;
	else for ( p->s = 3 ; mpz_ui_kronecker (p->s,P) >= 0 ; p->s+=2 );
	mpz_set_ui (p->t[0], p->s);
	mpz_invert (p->si, p->t[0], P);

	mpz_set_ui (p->g2, p->s);
	zp_exp (p->g2, p->g2, p->o, p);
	mpz_sub_ui (p->o, p->o, 1);
	mpz_div_2exp (p->o, p->o, 1);		// P->o = (m-1)/2, where p-1=2^k*m
	// g3 and z3 are set if and when needed	
}

void zp_mod_clear (zp_mod_t p)
{
	mpz_clear (p->p); mpz_clear (p->pm1); mpz_clear (p->o); mpz_clear (p->g2); mpz_clear (p->si); mpz_clear (p->t1);
	mpz_clear (p->t[0]); mpz_clear (p->w[0]); mpz_clear (p->g3[0]); mpz_clear (p->z3[0]);
	mpz_clear (p->t[1]); mpz_clear (p->w[1]); mpz_clear (p->g3[1]); mpz_clear (p->z3[1]);
}


// o[i] = 1/a[i] mod p  for i = 0 to n-1.  Returns 0 if any a[i] is 0 mod p, 1 otherwise
int mpz_array_invert (mpz_t o[], mpz_t a[], int n, mpz_t p)
{
	mpz_t *w;
	int i, sts;
	
	assert ( n >= 0 );
	if ( !n ) return 1;
	if ( n==1 ) return mpz_invert (o[0], a[0], p);
	if ( o == a ) { w = (mpz_t *) zp_malloc (n * sizeof(mpz_t)); for ( i = 0 ; i < n ; i++ ) mpz_init_modp_w (w[i], p); } else w = o;
	mpz_mul_mod (w[1], a[0], a[1], p);
	for ( i = 2 ; i < n ; i++ ) mpz_mul_mod (w[i], w[i-1], a[i], p);
	sts = mpz_invert (w[0], w[n-1], p);
	if ( sts ) {
		for ( i = n-1 ; i > 1 ; i-- ) {
			mpz_mul_mod (w[i], w[i-1], w[0], p);
			mpz_mul_mod (w[0], w[0], a[i], p);
			mpz_set (o[i], w[i]);
		}
		mpz_mul_mod (w[1], a[0], w[0], p);
		mpz_mul_mod (w[0], a[1], w[0], p);
		mpz_set (o[1], w[1]);
		mpz_set (o[0], w[0]);
	}
	if ( w != o ) { for ( i = 0 ; i < n ; i++ ) mpz_clear (w[i]); zp_free (w, n* sizeof(mpz_t)); }
	return sts;
}


// simple dual binary exponentiation, faster than 2 calls to mpz_powm when p is large (say > 3000 bits)
// a may coincide with b1 or b2 (but not both!)
void _zp_exp2 (zp_t b1, zp_t b2, zp_t a, mpz_t e1, mpz_t e2, zp_mod_t p)
{
	mpz_t x, b3;
	int i, j, n;

	if ( zp_is_zero(a) ) { mpz_set_ui(b1,e1?0:1); mpz_set_ui(b2,e2?0:1); return; }
	zp_init (x, p);  zp_init (b3, p);
	i = mpz_sizeinbase (e1, 2);
	n = mpz_sizeinbase (e2,2);
	if ( i > n ) n = i;
	zp_set (x, a);
	zp_set_one (b1);  zp_set_one (b2); zp_set_one (b3);
	for ( i = 0 ; i < n ; i++ ) {
		j = mpz_tstbit(e1,i) + (mpz_tstbit(e2,i)<<1);
		switch (j) {
		case 1: zp_mul (b1, b1, x, p); break;
		case 2: zp_mul (b2, b2, x, p); break;
		case 3: zp_mul (b3, b3, x, p); break;
		}
		if ( i < n-1 ) zp_mul (x, x, x, p);
	}
	zp_mul (b1, b1, b3, p);
	zp_mul (b2, b2, b3, p);
	zp_clear (b3); zp_clear (x);
}

	