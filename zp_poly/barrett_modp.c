#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include "zp_poly_util.h"

// Algorithm 9.2.8 in Crandall+Pomerance
void barrett_z_recip (mpz_t r, mpz_t n)
{
	mpz_t s,y;
	int b;

	b = mpz_sizeinbase(n,2);						// we want sizeinbase(n-1,2), but we'll assume n is not a power of 2
	mpz_init2(s,2*b);  mpz_init2 (y,2*b);
	mpz_set_ui(r,1);  mpz_mul_2exp (r, r, b);			// r = 2^b
	// Newton loop
	do {
		mpz_set (s,r) ;							// s = r
		mpz_mul (y, r, r);  mpz_div_2exp (y, y, b);
		mpz_mul (y, y, n); mpz_div_2exp (y, y, b);		// y = floor((n*floor(r^2/2^b))/2^b)
		mpz_mul_2exp (r,r, 1); mpz_sub (r, r, y);		// r = 2r - y
	} while ( mpz_cmp(r,s) > 0 );
	// final adjustment
	mpz_set_ui(y,1); mpz_mul_2exp (y, y, 2*b);
	mpz_mul (s, n, r);  mpz_sub(y, y, s);
	while ( mpz_cmp_ui (y,0) < 0 ) {
		mpz_sub_ui (r, r, 1);
		mpz_add (y, y, n);
	}
	mpz_clear (s); mpz_clear (y);
}

void barrett_mpz_mod (mpz_t x, mpz_t n, mpz_t R, mpz_t d)
{
	int s;
	
	s = 2*(mpz_sizeinbase(R,2)-1);
	for (;;) {
		mpz_mul (d, x, R);  mpz_div_2exp (d,d,s);
		mpz_submul (x, n, d);
		if ( mpz_cmp (x,n) < 0 ) return;
		mpz_sub (x, x, n);
		mpz_add_ui (d, d, 1);
	}
}


int main (int argc, char *argv[])
{
	clock_t start, end;
	mpz_t p, x, y, z, R, d;
	int pbits;
	register int i, cnt;
	
	if ( argc < 2 ) { puts ("barrett_modp pbits iterations"); return 0; }
	pbits = atoi(argv[1]);  cnt = atoi(argv[2]);
	mpz_init2 (p, pbits);  zp_random_prime (p, pbits);
	mpz_init2 (x, 2*pbits); mpz_init2 (y, 2*pbits); mpz_init2 (z, 2*pbits); mpz_init2(R, 2*pbits); mpz_init2 (d, 2*pbits);
	zp_random (x, p); zp_random (y, p);
	start = clock();
	for ( i = 0 ; i < cnt ; i++ ) mpz_add (z, x, y);
	end = clock();
	printf ("%d adds took %.3f secs\n", cnt, (end-start)/1000000.0);
	start = clock();
	for ( i = 0 ; i < cnt ; i++ ) mpz_add_mod (z, x, y, p);
	end = clock();
	printf ("%d add mods took %.3f secs\n", cnt, (end-start)/1000000.0);
	start = clock();
	for ( i = 0 ; i < cnt ; i++ ) mpz_mul (z, x, y);
	end = clock();
	printf ("%d muls took %.3f secs\n", cnt, (end-start)/1000000.0);
	mpz_mul (x, x, y);
	start = clock();
	for ( i = 0 ; i < cnt ; i++ ) { mpz_set (y, x); mpz_mod (z, y, p); }
	end = clock();
	printf ("%d mods took %.3f secs\n", cnt, (end-start)/1000000.0);
	barrett_z_recip (R, p);
	start = clock();
	for ( i = 0 ; i < cnt ; i++ ) { mpz_set (y, x); barrett_mpz_mod (y, p, R, d); }
	end = clock();
	printf ("%d Barrett mods took %.3f secs\n", cnt, (end-start)/1000000.0);
	
}
