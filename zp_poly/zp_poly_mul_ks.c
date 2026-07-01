/*

zp_poly_mul_ks.c: auxiliary multiplication module for zp_poly using naive KS

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
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include "zp_poly_main.h"
#include "zp_poly_util.h"

void zp_array_mul_mod_xn_ks (mpz_t h[], mpz_t f1[], int d1, mpz_t f2[], int d2, int m, int n, mpz_t p);


int main (int argc, char *argv[])
{
	clock_t start, end;
	zp_poly_mod_t mod;
	int m, n, s, d, k1, k2;
	mpz_t *f, *g, *h, *h1, *h2, p;
	unsigned long t1, t2;
	int i;
	
	if ( argc < 2 ) { puts ("Please specify log_2(deg*pbits)"); return 0; }
	s = atoi (argv[1]);
	if ( s < 6 || s> 30 ) { puts ("s must be at least 6 and no more than 30"); return 0; }
	
	mpz_init(p);
	for ( m = 6 ; m <= 16 ; m++ ) {
		zp_random_prime (p, 1<<m);
		for ( n = 6 ; m+n <= s && n <= m+8 ; n++ ) {
			d = 1<<n;
			f = zp_poly_alloc (d-1, p);  zp_poly_randomize (f, d-1, p);
			g = zp_poly_alloc (d-1, p); zp_poly_randomize (g, d-1, p);
			h = zp_poly_alloc (d, p); zp_poly_randomize (h, d-1, p);  mpz_set_ui (h[d], 1);
			h1 = zp_poly_alloc (d-1, p);  h2 = zp_poly_alloc (d-1, p);
			start = clock();  k1 = 0;
			do {
				zp_array_mul_mod_xn (h1, f, d, g, d, 0, d, p);  k1++;
				end = clock();
			} while ( (end-start) < CLOCKS_PER_SEC/10);
			t1 = (unsigned long) (end-start);
			start = clock();  k2 = 0;
			do {
				zp_array_mul_mod_xn_ks (h2, f, d, g, d, 0, d, p);  k2++;
				end = clock();
			} while ( (end-start) < CLOCKS_PER_SEC/10 );
			t2 = (unsigned long) (end-start);
			printf ("pbits=%d, deg=%d, array mul time (std) =%f  array mul time (ks) = %f\n", 1<<m, d, (double)t1/(k1*CLOCKS_PER_SEC), (double)t2/(k2*CLOCKS_PER_SEC));
			if ( ! zp_poly_equal (h1, d-1, h2, d-1) ) {
				for ( i = 0 ; i < d ; i++ ) if ( mpz_cmp(h1[i], h2[i]) ) break;
				printf ("ks array mul failed at coefficient %d!\n", i);
			}
			start = clock();
			zp_poly_mod_init (mod, h, d, p, ZP_POLY_MOD_ALGO_BARRETT);
			end = clock();
			printf ("pbits = %d, deg=%d, mod init time (std) %f\n", 1<<m, d, (double)(end-start)/CLOCKS_PER_SEC);
			start = clock();  k1 = 0;
			do {
				zp_poly_mod_mul (h1, f, g, mod);  k1++;
				end = clock();
			} while ( (end-start) < CLOCKS_PER_SEC/10 );
			t1 = (unsigned long) (end-start);
			zp_poly_mod_clear (mod);
			start = clock();
			zp_poly_mod_init (mod, h, d, p, ZP_POLY_MOD_ALGO_KRONECKER);
			end = clock();
			printf ("pbits = %d, deg=%d, mod init time (ks) %f\n", 1<<m, d, (double)(end-start)/CLOCKS_PER_SEC);
			start = clock();  k2 = 0;
			do {
				zp_poly_mod_mul (h2, f, g, mod);  k2++;
				end = clock();
			} while ( (end-start) < CLOCKS_PER_SEC/10 );
			t2 = (unsigned long) (end-start);
			zp_poly_mod_clear (mod);
			printf ("pbits=%d, deg=%d, mod mul time (std) =%f  mod mul time (ks) = %f\n", 1<<m, d, (double)t1/(k1*CLOCKS_PER_SEC), (double)t2/(k2*CLOCKS_PER_SEC));
			if ( ! zp_poly_equal (h1, d-1, h2, d-1) ) printf ("ks mod mul failed!\n");

			zp_poly_free (f, d-1); zp_poly_free (g,d-1); zp_poly_free (h, d); zp_poly_free (h1,d-1); zp_poly_free (h2,d-1);
		}
	}
}
