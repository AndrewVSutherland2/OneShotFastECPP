/*

zp_poly_div.c: division module for zp_poly

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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"

/*
    Compute the inverse of f mod x^n using Alg 9.3 in [GG],
    with modification to more efficiently handle n not a power of 2.
    Does not require f monic
*/
void zp_poly_inverse_mod_xn (zp_t h[], zp_t f[], int d, int n, zp_p p)
{
    zp_t *w, *g;
    register int e, i, j, k, m, r;

    assert ( n > 0 && d >= 0 );
    assert ( mpz_sgn(f[0]) != 0 );
    // handle h = f, but otherwise assume h and f don't overlap
    if ( h == f ) g = zp_poly_alloc(n,p); else g=h; 
    mpz_invert(g[0],f[0],p);    // we assume mpz_invert is smart about 1
    w = zp_poly_alloc (2*n-1,p);
    for ( r = 1 ; (1<<r) < n ; r++ );
    for ( m=1,r-- ; m < n; m=e, r-- ) {             // g  is the inverse of f mod x^m
        e = n>>r;  if ( n & ((1<<r)-1) ) e++;       // e = ceil (n/2^r) will be the next m
        if ( e > n ) e = n;                     // but e never needs to exceed n
        zp_poly_mul(w,g,m-1,g,m-1,p);           // w = g^2
        j = ( d < e ? d : e-1 );                    // j = deg (f mod x^e)
        k = ( 2*m-2 < e ? 2*m-2 : e-1 );            // k = deg (w mod x^e)
        zp_poly_mul (w,f,j,w,k,p);              // w = fg^2 mod x^e
        for ( i = 0 ; i < m ; i++ ) {               // set g = 2g - fg^2 mod x^e
            mpz_mul_2exp (g[i],g[i],1);
            mpz_sub (g[i],g[i],w[i]);
            mpz_mod_small (g[i],g[i],(zp_p_to_mpz)p);
        }
        for ( ; i < e ; i++ ) zp_neg(g[i],w[i],p);
    }
    zp_poly_free (w,2*n-1);
    if ( h==f ) { zp_poly_copy(h,g,n-1); zp_poly_free(g,n); }
}


/*
    Given a poly f of deg d, computes the quotient in the Euclidean
    division of x^n by f, i.e.  x^n = h*f + r with deg r < deg f.
    Implemented as a special case of Alg 9.5 in [GG], doesn't require f monic.
*/
void zp_poly_div_xn (zp_t q[], zp_t f[], int d, int n, zp_p p)
{
    zp_poly_rev_inplace (f,d);                          // modify f
    zp_poly_inverse_mod_xn (q, f, d, n-d+1, p);
    if ( f != q ) zp_poly_rev_inplace (f,d);                // restore f if needed
    zp_poly_rev_inplace (q, n-d);
}

/*
    Computes the quotient q in the Euclidean division a/b,
    i.e. the poly q such that a = qb+r with deg r < deg b. 
    Use zp_poly_rem to get r, if needed.  Based on Alg 9.5 in [GG],
           but does not assume b is monic.
*/
int zp_poly_div (zp_t q[], zp_t a[], int da, zp_t b[], int db, zp_p p)
{
    zp_t *w;
    int m;
    
    assert (db >= 0);
    if ( da < db ) return -1;
    m = da - db;
    if ( m <= ZP_POLY_VERY_SMALL_DEGREE ) {
        if ( q ==b ) w = zp_poly_alloc(m,p); else w = q;        // handle aliasing case not handled in div_small
        zp_poly_div_small (w, a, da, b, db, p);
        if ( q == b ) { zp_poly_copy(q,w,m); zp_poly_free(w,m); }
        return m;
    }
    if ( q == a || q==b ) w = zp_poly_alloc(m,p); else w = q;   // handle aliasing
    zp_poly_rev_inplace (b,db);                         // modify b
    zp_poly_inverse_mod_xn (w, b, db, m+1, p);
    if ( b != q ) zp_poly_rev_inplace (b,db);                   // restore b if needed
    zp_poly_rev_inplace (a,da);                         // modify a
    zp_poly_mul_mod_xn (w, a, m, w, m, m+1, p);
    if ( a != q ) zp_poly_rev_inplace (a,da);                   // restore a if needed
    zp_poly_rev_inplace (w,m);
    if ( q == a || q == b ) { zp_poly_copy(q,w,m); zp_poly_free(w,m); }
    return m;
}

// r = a - qb, where q is the Euclidean quotient of a/b, so q=ab+r with deg r < deg b, r may coincide with a
// r is zero-padded to degree db-1, returns deg r
int zp_poly_rem (zp_t r[], zp_t a[], int da, zp_t b[], int db, zp_t q[], zp_p p)
{
    zp_t *w;
    int d;
    register int i;
    
    assert (db >= 0);
    if ( da < db ) { zp_poly_copy (r, a, da);  for ( i = da+1 ; i < db ; i++ ) mpz_set_zero (r[i]);  return da; }   // handles a=0 case correctly
    if ( da-db <= ZP_POLY_VERY_SMALL_DEGREE ) return zp_poly_rem_small (r, a, da, b, db, q, p);
    if ( r == a ) w = zp_poly_alloc(db-1,p); else w = r;            // handle r==a
    d = da-db+1;  if ( d > db ) d = db;                     // d = min(dq+1,db)
    zp_array_mul_mod_xn (w, b, db, q, d, 0, db, p);
    for ( i = 0 ; i < db ; i++ ) mpz_sub_mod(w[i],a[i],w[i],p);
    if ( r == a ) { zp_poly_copy(r,w,db-1); zp_poly_free(w,db-1); }
    return zp_poly_degree (r, db-1);
}

// r = a mod b, r is zero-padded to degree db-1, returns deg r, r may coincide with a
int zp_poly_mod (zp_t r[], zp_t a[], int da, zp_t b[], int db, zp_p p)
{
    zp_t *q;
    int dr;
    
    if ( da < db ) { zp_poly_copy (r, a, da); return da; }
    q = zp_poly_alloc (da-db, p);
    zp_poly_div (q, a, da, b, db, p);
    dr = zp_poly_rem (r, a, da, b, db, q, p);
    zp_poly_free (q, da-db);
    return dr;
}

// input: r = poly a of deg < d*2^n, zero padded to deg d*2^n-1
//           b[i] = b^{2^i} for i < n, where b is a poly of deg d
// output: r = the radix b rep of a
void _zp_poly_radix (zp_t *r, zp_t *b[], int d, int n, zp_p p)
{
    zp_t *w;
    int t;

    // base case
    if ( ! n ) return;
    
    // replace r by its 2-digit radix b^{2^{n-1}} representation
    t = 1<<(n-1);
    w = zp_poly_alloc (d*t-1, p);
    zp_poly_div (w, r, d*2*t-1, b[n-1], d*t, p);
    zp_poly_rem (r, r, d*2*t-1, b[n-1], d*t, w, p);
    zp_poly_copy (r+d*t, w, d*t-1);
    zp_poly_free (w, d*t-1);
    
    // recursively represent both digits in radix b
    _zp_poly_radix (r, b, d, n-1, p);
    _zp_poly_radix (r+d*t, b, d, n-1, p);
}

// r = radix rep of a in base b, that is a = sum r_i b^i
// where each r_i has deg < db and there are ceil(da/db) terms
// the output is the concatenation of da/db+1 polys of deg db-1 (each zero padded)
// with total size db * (da/db+1) < da + db
void zp_poly_radix (zp_t r[], zp_t a[], int da, zp_t b[], int db, zp_p p)
{
    zp_t **e, *c, *ep, *w;
    int i, d, m, n, k;
    
    assert (db > 0);
    if ( da < db ) { zp_poly_copy (r, a, da); return; }
    
    // compute k=2^n > da/db
    k = da/db +1;
    n = ui_bits (k);  k = 1<<n;
    
    // compute b, b^2, b^4, ..., b^{k/2}, storing coefficients in the array c indexed by e (so e[i] points to b^{2^i} of deg 2^i*db
    e = zp_malloc (n*sizeof(mpz_t *));
    c = zp_array_alloc (m=db*(k-2)+n-1, zp_bits(p)+ZP_PBITS_FUDGE);
    e[0] = b;  d = db;  ep = c;
    for ( i = 1 ; i < n ; i++ ) { e[i] = ep;  zp_poly_mul (e[i], e[i-1], d, e[i-1], d, p);  d += d; ep += d+1; }
    assert ( ep-c == m );
    w = zp_array_alloc (k*db, zp_bits(p)+ZP_PBITS_FUDGE);
    zp_poly_copy (w, a, da);  zp_poly_zero_pad (w, da+1, k*db-1);
    _zp_poly_radix (w, e, db, n, p);
    zp_poly_copy (r, w, (da/db+1)*db-1);
    zp_array_free (w, k*db);
    zp_array_free (c, m);
    zp_free (e, n*sizeof(mpz_t *));
}

#ifdef ZP_POLY_RADIX_TEST
#include <time.h>
 
// TODO: switch mpz_t to zp_t in the code below

int main (int argc, char *argv[])
{
    clock_t start, end;
    zp_mod_t pmod;
    mpz_t p, *r, *a, *b, *c, *g, *h, x, y, z, w, s;
    int i, j, k, d, mind, maxd, db, delta, pbits;
    
    if ( argc < 5 ) { puts ("zp_poly_radix_test mind maxd radixd pbits"); return 0; }
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);
    mind = atoi (argv[1]); maxd = atoi (argv[2]);  db = atoi(argv[3]);  pbits = atoi(argv[4]);
    mpz_init2(p,pbits);  zp_random_prime (p,pbits);
    zp_mod_init (pmod, p);
    mpz_init_modp_w (x, p);  mpz_init_modp_w (y, p);  mpz_init_modp_w (z, p);  mpz_init_modp_w (w, p);  mpz_init_modp_w (s, p);
    
    r = zp_poly_alloc (maxd + db, p); a = zp_poly_alloc (maxd, p);  b = zp_poly_alloc (db, p); c = zp_poly_alloc (maxd, p);  g = zp_poly_alloc (maxd + db, p);  h = zp_poly_alloc (maxd, p);
    zp_poly_randomize (b, db, p);
    
    for ( i = delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i; 
    if ( ! delta ) delta = 1;  // a really stupid GCD   
    for ( d = mind ; d <= maxd ; d += delta ) {
        zp_poly_randomize (a, d, p);
        k = d/db + 1;   // # of digits in radix b rep of a
        start = clock();
        zp_poly_radix (r, a, d, b, db, p);
        end = clock();
        printf ("Writing degree %d poly in radix of degree %d took %.3f secs\n", d, db, (double)(end-start)/CLOCKS_PER_SEC);
        for ( i = 0 ; i < 10 ; i++ ) {
            zp_randomm (x, p);
            zp_poly_eval (y, b, db, x, pmod);
            mpz_set_one (w);  mpz_set_zero (s);
            for ( j = 0 ; j < k ; j++ )
                { zp_poly_eval (z, r + j*db, db-1, x, pmod);  mpz_mul_mod (z, z, w, p);  mpz_add_mod (s, s, z, p);  mpz_mul_mod (w, w, y, p);  }
            zp_poly_eval (y, a, d, x, pmod);
            if ( mpz_cmp (y, s) ) { printf ("Random verification failed!\n"); abort(); }
        }
        start = clock();
        zp_poly_zero_fill (g, k*db);
        for ( i = 0 ; i < db ; i++ ) {
            for ( j = 0 ; j < k ; j++ ) mpz_set (c[j],r[j*db+i]);
            zp_poly_compose_horner (h, c, k-1, b, db, pmod);
            zp_poly_add (g+i, g+i, (k-1)*db, h, (k-1)*db, p);
        }
        end = clock();
        printf ("Reading degree %d poly in radix of degree %d took %.3f secs\n", d, db, (double)(end-start)/CLOCKS_PER_SEC);
        if ( ! zp_poly_equal (g, d, a, d) ) { printf ("composition of radix rep failed!\n"); abort(); }
    }
    zp_poly_free (r, maxd+db);  zp_poly_free (a, maxd);  zp_poly_free (b, db); zp_poly_free (c, maxd); zp_poly_free (g, maxd+db); zp_poly_free (h, maxd);
    mpz_clear (x); mpz_clear (y); mpz_clear (z); mpz_clear (w); mpz_clear (s); mpz_clear (p);
    zp_mod_clear (pmod);
    zp_cleanup();
    zp_mem_report (1);
}

#endif

#ifdef ZP_POLY_INV_XN_TEST
#include <time.h>
#define cpusecs(beg,end)    (((double)(end-beg))/CLOCKS_PER_SEC)

int main (int argc, char *argv[])
{
    clock_t start,stop;
    mpz_t p, *f, *g, *h;
    int d, b, n, dh;
    
    if ( argc < 4 ) { puts ("zp_poly_inv_xn_test d n pbits"); return 0; }
    
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);
    
    d = atoi(argv[1]);  n = atoi(argv[2]);  b = atoi(argv[3]);
    mpz_init2(p,b);  zp_random_prime (p,b);
    
    f = zp_poly_alloc(d,p);  zp_poly_randomize (f,d, p);
    if ( ! mpz_sgn(f[0]) ) mpz_set_ui(f[0],1);
    g = zp_poly_alloc(n,p);
    h = zp_poly_alloc(d+n,p);
    
    mpz_set_ui(g[0],1); mpz_set_ui(g[1],1);
    zp_poly_mul(h,g,1,g,1,p);
    gmp_printf ("(x+1)^2 = %Zd*x^2 +%Zd*x + %Zd\n", h[2], h[1], h[0]);
    
    printf ("calling zp_poly_inverse_mod_xn(g,f,%d,%d,p)\n", d, n);
    
    start = clock();
    zp_poly_inverse_mod_xn (g, f, d, n, p);
    stop = clock();
    printf ("zp_poly_inverse_mod_xn(g,f,%d,%d,%d-bit prime) took %.3f secs\n", d,n,b,cpusecs(start,stop));
    zp_poly_mul (h, f, d, g, n-1,p);  dh = zp_poly_degree (h,n-1);
    if ( ! zp_poly_is_one(h,dh) ) {
        printf ("zp_poly_inverse_mod_xn(g,f,%d,%d,p) failed!\n", d, n);
        printf ("f =  ");  zp_poly_print(f,d); puts ("");
        printf ("g =  ");  zp_poly_print(g,n-1); puts ("");
        printf ("f*g mod x^n =  ");  zp_poly_print(h,dh); puts ("");
        gmp_printf ("Modulus p = %Zd\n", p);
    }
    zp_cleanup();
    zp_mem_report(1);
}

#endif

#ifdef ZP_POLY_DIV_TEST
#include <time.h>
#define cpusecs(beg,end)    (((double)(end-beg))/CLOCKS_PER_SEC)

void zp_poly_div_verify (zp_t w[], zp_t a[], int da, zp_t b[], int db, zp_t q[], zp_p p)
{
    int i;
    
    zp_poly_mul (w, b, db, q, da-db, p);
    for ( i = da ; i >= db ; i-- ) if ( mpz_cmp(w[i],a[i]) ) break;
    if ( i >= db ) {
        printf ("zp_poly_div(q,a,%d,b,%d,p) failed!\n", da, db);
        printf ("a=  ");  zp_poly_print(a,da); puts ("");
        printf ("b=  ");  zp_poly_print(b,db); puts ("");
        printf ("q =  ");  zp_poly_print(q,da-db); puts ("");
        printf ("q*b =  ");  zp_poly_print(w,da); puts ("");
        gmp_printf ("Modulus p = %Zd\n", p);
        abort();
    }
}

void random_div_test (zp_t q[], zp_t a[], int da, zp_t b[], int db, zp_p p, zp_t w[])
{
    zp_poly_randomize (a,da,p);  zp_poly_randomize (b,db,p);
    zp_poly_div (q,a,da,b,db,p);
    zp_poly_div_verify(w,a,da,b,db,q,p);
}

// test division a/b using every poly of a of deg <=d and monic b with deg b <= deg a
void full_div_test (int d, zp_p p)
{
    zp_t *q, *a, *b, *w;
    int da, db;
    
    q = zp_poly_alloc(d,p);  a = zp_poly_alloc(d,p);  b = zp_poly_alloc(d,p);  w = zp_poly_alloc(d,p);
    for ( da = 1 ; da <= d ; da++ ) {
        zp_poly_zero_fill(a,da-1); mpz_set_ui(a[da],1);
        do {
            for ( db = 1 ; db <= da ; db++ ) {
                zp_poly_zero_fill(b,db-1);  mpz_set_ui(b[db],1);
                do {
                    zp_poly_div (q,a,da,b,db,p);
                    zp_poly_div_verify (w,a,da,b,db,q,p);
                } while ( zp_poly_next (b,db,p) );
            }
        } while ( zp_poly_next (a,da,p) );
        printf ("\r%d\r",da); fflush (stdout);
    }
    zp_poly_free(q,d);  zp_poly_free(a,d);  zp_poly_free(b,d);  zp_poly_free(w,d);
}

int main (int argc, char *argv[])
{
    clock_t start,stop;
    mpz_t p, *a, *b, *q, *w;
    int da, db, pbits, n, i, cnt, monic;
    
    if ( argc >1 && argc < 4 ) { puts ("zp_poly_div_test da db pbits [monic iterations]"); return 0; }
    
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);
    
    if ( argc == 1 ) {
        mpz_init2(p,2000); mpz_set_one(p);;  mpz_mul_2exp(p,p,2000);
        q = zp_poly_alloc(2000,p);  a = zp_poly_alloc(2000,p);  b= zp_poly_alloc(2000,p);  w = zp_poly_alloc(2000,p);
        mpz_set_ui(p,3);
        puts ("Running exhaustive division tests of small degrees up to 7 using p = 3...");
        full_div_test (7, p);
        puts ("Running one random division test with all degree combinations up to 300 using p=3...");
        for ( da = 0 ; da <= 300 ; da++ ) { for ( db = 0 ; db <= da ; db++ ) random_div_test (q,a,da,b,db,p,w); printf ("\r%d\r", db); fflush(stdout); }
        puts ("Running exhaustive division tests of small degrees up to 5 using p = 5...");
        mpz_set_ui(p,5);
        full_div_test (4, p);
        puts ("Running one random division test with all degree combinations up to 300 using p=5...");
        for ( da = 0 ; da <= 300 ; da++ ) { for ( db = 0 ; db <= da ; db++ ) random_div_test (q,a,da,b,db,p,w); printf ("\r%d\r", db); fflush(stdout); }
        puts ("Running one random division test 20n / 10n for n up to 100 and primes of size 100, 200, ..., 2000 bits...");
        for ( pbits = 100; pbits <= 2000 ; pbits+=100 ) {
            zp_random_prime (p,pbits); printf ("%lu-bit prime tests...\n", zp_bits(p));
            for ( n = 1 ; n <= 100 ; n++ ) { random_div_test (q,a,20*n,b,10*n,p,w); printf ("\r%d/%d\r", 20*n,10*n); fflush(stdout); }
        }
        puts ("All tests passed!");
        mpz_clear (p); zp_poly_free (q,2000);  zp_poly_free (a,2000);  zp_poly_free (b,2000);  zp_poly_free (w,2000);  
        zp_mem_report (1);
        return 0;
    }
    
    da = atoi(argv[1]);  db = atoi(argv[2]);  pbits = atoi(argv[3]);
    if ( da < db ) { puts ("Please make da >= db"); return 0; }
    if ( argc > 4 ) monic = atoi(argv[4]); else monic = 1;
    if ( argc > 5 ) cnt = atoi(argv[5]); else cnt = 1;
    
    mpz_init2(p,pbits);
    zp_random_prime (p,pbits);  pbits = zp_bits(p);
    
    a = zp_poly_alloc(da,p);
    b = zp_poly_alloc(db,p);
    q = zp_poly_alloc(da-db,p);
    w = zp_poly_alloc(da,p);

    printf ("calling zp_poly_div(q,a,%d,b,%d,%d-bit prime)", da, db, pbits);
    if ( cnt > 1 ) printf (" %d times\n", cnt); else puts ("");
    
    if ( cnt == 1 ) {
        mpz_set_ui(b[db],1);
        zp_poly_randomize (b, (monic?db-1:db), p);
        zp_poly_randomize (a, da, p);
        start = clock();
        zp_poly_div (q,a,da,b,db,p);
        stop = clock();
        printf ("zp_poly_div(q,a,%d,b,%d,%d-bit prime) took %.3f secs\n", da,db,pbits,cpusecs(start,stop));
        zp_poly_div_verify (w, a, da, b, db, q, p);
    } else {
        start = clock();
        mpz_set_ui(b[db],1);
        for ( i = 0 ; i < cnt ; i++ ) {
            zp_poly_randomize (b, (monic?db-1:db), p);
            zp_poly_randomize (a, da, p);
            zp_poly_div (q, a, da, b, db, p);
            zp_poly_div_verify (w, a, da, b, db, q, p);
        }
        stop = clock();
        printf ("%d verified calls to zp_poly_div(q,a,%d,b,%d,%d-bit prime) took %.3f secs\n", cnt,da,db,pbits,cpusecs(start,stop));
    }
    puts ("All tests passed");
    mpz_clear (p); zp_poly_free (a,da);  zp_poly_free (b,db);  zp_poly_free (q,da-db);  zp_poly_free (w,da);
    zp_cleanup();
    zp_mem_report (1);
}

#endif
