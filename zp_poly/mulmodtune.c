/*

modmultune.c: performance testing program for arithmetic in the ring Fp[X]/(g(X))

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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"

int MINB = 1000;
#define MAXB			2000
#define DELTAB		MINB
#define MAXN			8

static inline double secs (clock_t start, clock_t end)	{ return (double)(end-start)/CLOCKS_PER_SEC; }

int main (int argc, char *argv[])
{
	clock_t start, end;
	double cost;
	zp_poly_mod_t mod;
	mpz_t p, *g, *f1, *f2, *h, *f[MAXN];
	int d, pbits;
	register int i, n;
	
	if ( argc > 1 ) MINB = atoi(argv[1]);
	
	mpz_init2(p,MAXB);
	puts ("Timing loop starting...\n");
	for ( pbits = MINB ; pbits <= MAXB ; pbits+=DELTAB ) {
		zp_random_prime(p,pbits);
		f1 = zp_poly_alloc (pbits, p);  f2 = zp_poly_alloc (pbits, p); h = zp_poly_alloc (pbits, p); g = zp_poly_alloc (pbits+1, p); 
		for ( i = 0 ; i < MAXN ; i++ ) f[i] = zp_poly_alloc (pbits, p);
		for ( d = pbits/2 ; d <= pbits ; d += pbits/4 ) {
			zp_poly_randomize (f1, d, p); zp_poly_randomize (f2, d, p);  zp_poly_randomize (g, d, p);  zp_poly_randomize (h, d, p);
			for ( i = 0 ; i < MAXN ; i++ ) zp_poly_randomize (f[i], d, p);
			mpz_set_one (g[d]);
			start = clock();
			zp_poly_mod_init (mod, g, d, p, 0);
			end = clock();
			printf ("zp_poly_mod_init    deg %5d at %d bits: %.3f secs\n", d, pbits, secs(start,end));
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_mod_add (h,f1,f2,mod);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			cost = secs(start,end)/i;
			printf ("zp_poly_mod_add     deg %5d at %d bits: %.3f usecs\n", d, pbits, cost);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_mod_mul (h,f1,f1,mod);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			cost = secs(start,end)/i;
			printf ("zp_poly_mod_sqr     deg %5d at %d bits: %.3f usecs\n", d, pbits, cost);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_mod_mul (h,f1,f2,mod);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			cost = secs(start,end)/i;
			printf ("zp_poly_mod_mul     deg %5d at %d bits: %.3f usecs\n", d, pbits, cost);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_mod_inv (h,f1,mod);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			cost = secs(start,end)/i;
			printf ("zp_poly_mod_inv     deg %5d at %d bits: %.3f usecs\n", d, pbits, cost);
			for ( n = 2 ; n <= MAXN ; n *= 2 ) {
				start = clock();
				for ( i = 1 ; ; i++ ) {
					zp_poly_mod_array_inv (f,f,n,mod);
					end = clock();
				       if ( end-start >= 500000L ) break;
				}
				cost = secs(start,end)/(i*n);
				printf ("zp_poly_mod_inv(%d)  deg %5d at %d bits: %.3f usecs\n", n, d, pbits, cost);
			}
			zp_poly_mod_clear (mod);
		}
		zp_poly_free (f1, pbits);  zp_poly_free (f2, pbits);  zp_poly_free (h, pbits);  zp_poly_free (g, pbits+1);  
		for ( i = 0 ; i < MAXN ; i++ ) zp_poly_free (f[i], pbits);
	}
	zp_mem_report(1);
}
