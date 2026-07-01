/*

multune.c: very basic performance testing program for zp_poly mul/div/gcd operations

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
#include <math.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"

#define LOG2			0.69314718055994530941723212145818
#define log2(x)		(log(x)/LOG2)

int MINB = 1000;		// default, may be overridden on the command line
#define MAXB		20000
#define DELTAB	MINB

int main (int argc, char *argv[])
{
	clock_t start, end;
	long mcost;
	mpz_t p, w, *f1, *f2, *h1, *h2, *h3, *h4, *a00, *a01, *a10, *a11, *b00, *b01, *b10, *b11;
	int d, d1, d2, d3, pbits;
	register int i;
	
	if ( argc > 1 ) MINB = atoi(argv[1]);
	
	if ( ZP_POLY_STRASSEN_MIN_PBITS >= 0 ) {
		printf ("Warning, ZP_POLY_STRASSEN_MIN_PBITS = %d\n", ZP_POLY_STRASSEN_MIN_PBITS);
		printf ("#define ZP_POLY_STRASSEN_MIN_PBITS  -1 in zp_poly_mul.h to disable automatic use of Strassen\n");
	}
	
	mpz_init2(p,MAXB); mpz_init2(w,2*MAXB+64);
	puts ("Timing loop starting...\n");
	for ( pbits = MINB ; pbits <= MAXB ; pbits+=DELTAB ) {
		zp_random_prime(p,pbits);
		a00 =  zp_poly_alloc(pbits,p); a01 =  zp_poly_alloc(pbits,p); a10 =  zp_poly_alloc(pbits,p); a11 =  zp_poly_alloc(pbits,p); 
		b00 =  zp_poly_alloc(2*pbits,p); b01 =  zp_poly_alloc(2*pbits,p); b10 =  zp_poly_alloc(2*pbits,p); b11 =  zp_poly_alloc(2*pbits,p); 
		f1 = zp_poly_alloc(pbits,p); f2 = zp_poly_alloc(pbits,p); h1 = zp_poly_alloc(2*pbits,p);  h2 = zp_poly_alloc(2*pbits,p);  h3 = zp_poly_alloc(2*pbits,p);  h4 = zp_poly_alloc(2*pbits,p);
		for ( d = pbits/2 ; d <= pbits ; d += pbits/4 ) {
			zp_poly_randomize (a00,d,p);  zp_poly_randomize (a01,d,p);  zp_poly_randomize (a10,d,p);  zp_poly_randomize (a11,d,p);
			zp_poly_randomize (b00,d,p);  zp_poly_randomize (b01,d,p);  zp_poly_randomize (b10,d,p);  zp_poly_randomize (b11,d,p);
			zp_poly_randomize (f1,d,p); zp_poly_randomize (f2,d,p);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_mul (h1,f1,d,f2,d,p);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			mcost = (end-start)/i;
			printf ("zp_poly_mul   deg %5d at %d pbits: %.2f msecs\n", d, pbits, (double)mcost/1000);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_mul2 (h1,&d1,a00,d,f1,d,a01,d,f2,d,p);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			printf ("zp_poly_mul2  deg %5d at %d pbits: %.3f M\n", d, pbits, (double)((end-start)/i)/mcost);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_array_mul42_mod_xn (h3,h4,a00,d+1,a01,d+1,a10,d+1,a11,d+1,f1,d+1,f2,d+1,d+d+1,p); 
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			printf ("zp_poly_mul42 deg %5d at %d pbits: %.3f M\n", d, pbits, (double)((end-start)/i)/mcost);
			if ( ZP_POLY_STRASSEN_MIN_PBITS < 0 || pbits < ZP_POLY_STRASSEN_MIN_PBITS ) {
				start = clock();
				for ( i = 1 ; ; i++ ) {
					zp_array_mul4_mod_xn (b00,b01,b10,b11,a00,d+1,a01,d+1,a10,d+1,a11,d+1,b00,d+1,b01,d+1,b10,d+1,b11,d+1,d+d+1,p); 
					end = clock();
				       if ( end-start >= 500000L ) break;
				}
				printf ("zp_poly_mul4 deg %5d at %d pbits: %.3f M\n", d, pbits, (double)((end-start)/i)/mcost);
			}
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_array_mul4_mod_xn_strassen (b00,b01,b10,b11,a00,d+1,a01,d+1,a10,d+1,a11,d+1,b00,d+1,b01,d+1,b10,d+1,b11,d+1,d+d+1,p); 
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			printf ("zp_poly_mul4s deg %5d at %d pbits: %.3f M\n", d, pbits, (double)((end-start)/i)/mcost);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_div (h2,h1,2*d,f1,d,p);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			printf ("zp_poly_div   deg %5d at %d pbits: %.3f M\n", d, pbits, (double)((end-start)/i)/mcost);
			start = clock();
			for ( i = 1 ; ; i++ ) {
				zp_poly_xgcd (h1,&d1,h2,&d2,h3,&d3,f1,d,f2,d,p);
				end = clock();
			       if ( end-start >= 500000L ) break;
			}
			printf ("zp_poly_xgcd  deg %5d at %d pbits: %.3f M log2(d)\n", d, pbits, (double)((end-start)/i) / (mcost * log2(d)));
			puts ("");
		}
		zp_poly_free (a00, pbits); zp_poly_free (a01, pbits); zp_poly_free (a10, pbits); zp_poly_free (a11, pbits); 
		zp_poly_free (b00, pbits); zp_poly_free (b01, pbits); zp_poly_free (b10, pbits); zp_poly_free (b11, pbits); 
		zp_poly_free (f1,pbits);   zp_poly_free (f2,pbits);   zp_poly_free (h1,2*pbits); ;   zp_poly_free (h2,2*pbits);  zp_poly_free (h3,2*pbits);  zp_poly_free (h4,2*pbits);
	}
}
