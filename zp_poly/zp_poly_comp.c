#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"

// r(X) = f(g(X)) with degree df*dg (negative degrees => zero poly)
void zp_poly_compose_horner (zp_t r[], zp_t f[], int df, zp_t g[], int dg, zp_p p)
{
    zp_t *wf, *wg;
    int i, dr;
    
    // deal with trivial cases with df, dg <= 0
    if ( df < 0 || dg < 0 ) return;
    if ( df == 0 ) { zp_set (r[0], f[0]); return; }
    if ( dg == 0 ) { zp_poly_eval (r[0], f, df, g[0], p); return; }
    
    // copy inputs if necessary to deal with aliasing
    if ( r == f ) { wf = zp_poly_alloc (df, p);  zp_poly_copy (wf, f, df); f = wf; } else wf = 0;
    if ( r == g ) { wg = zp_poly_alloc (dg, p);  zp_poly_copy (wg, g, dg); g = wg; } else wg = 0;
    
    // compute r=f(g) using Horner's rule, with df-1 calls to zp_poly_mul
    zp_poly_copy (r, g, dg);
    zp_poly_smul (r, r, dg, f[df], p);  dr = dg;
    for ( i = df-1 ; i > 0 ; i-- ) {
        zp_add (r[0], r[0], f[i], p);
        zp_poly_mul (r, r, dr, g, dg, p);
        dr += dg;
    }
    zp_add (r[0], r[0], f[0], p);
    if ( wf == f ) zp_poly_free (wf, df);
    if ( wg == g ) zp_poly_free (wg, dg);
}

// r(X) = f(g(X)) mod h(X) with f of degree d, g reduced mod h (zero filled to deg(h)-1)
void zp_poly_mod_compose_horner (zp_t r[], zp_t f[], int d, zp_t g[], zp_poly_mod_t h)
{
    mpz_t *wf, *wg;
    int i;
    
    // deal with trivial cases with d <= 0
    if ( d < 0 ) { zp_poly_mod_set_zero (f,h); return; }
    if ( d == 0 ) { mpz_set (r[0], f[0]); zp_poly_zero_pad (r, 1, h->d-1); return; }
    
    // copy inputs if necessary to deal with aliasing
    if ( r == f ) { wf = zp_poly_alloc (d, h->p);  zp_poly_copy (wf, f, d); f = wf; } else wf = 0;
    if ( r == g ) { wg = zp_poly_mod_alloc (h);  zp_poly_mod_copy (wg, g, h); g = wg; } else wg = 0;
    
    // compute r=f(g) using Horner's rule, with df-1 calls to zp_poly_mul
    zp_poly_mod_copy (r, g, h);
    zp_poly_mod_smul (r, r, f[d], h);
    for ( i = d-1 ; i > 0 ; i-- ) {
        mpz_add_mod (r[0], r[0], f[i], h->p);
        zp_poly_mod_mul (r, r, g, h);
    }
    mpz_add_mod (r[0], r[0], f[0], h->p);
    if ( wf == f ) zp_poly_free (wf, d);
    if ( wg == g ) zp_poly_mod_free (wg, h);
}



// h = c[0]*f[0] + ... + c[n-1]*f[m-1] where c[i] are scalars[i] and f[i] are polys mod x^n (filled to deg n-1)
static inline void _zp_poly_linear_combo_mod_xn (mpz_t h[], mpz_t *f[], int m, mpz_t c[], int n, mpz_t p, mpz_t w)
{
    register int i, k;
    if ( m <= 0 || n <= 0 ) { zp_poly_zero_fill (h, n-1); return; }
    for ( k = 0; k < n ; k++) {
        mpz_mul (w, f[0][k], c[0]);
        for ( i = 1 ; i < m ; i++) mpz_addmul (w, f[i][k], c[i]);
        mpz_mod (h[k], w, p);
    }
}


// r(X) = f(g(X)) with degree df*dg (negative degrees => zero poly)
void zp_poly_mod_compose_window (zp_t r[], zp_t f[], int d, zp_t g[], zp_poly_mod_t h, int k)
{
    mpz_t *wf, *wg, **pow, *w, x;
    int i, j;

    // deal with trivial cases with d <= 0
    if ( d < 0 ) { zp_poly_mod_set_zero (f,h); return; }
    if ( d == 0 ) { mpz_set (r[0], f[0]); zp_poly_zero_pad (r, 1, h->d-1); return; }
    
    // copy inputs if necessary to deal with aliasing
    if ( r == f ) { wf = zp_poly_alloc (d, h->p);  zp_poly_copy (wf, f, d); f = wf; } else wf = 0;
    if ( r == g ) { wg = zp_poly_mod_alloc (h);  zp_poly_mod_copy (wg, g, h); g = wg; } else wg = 0;
    
    // pow[i] = g^{i+1} mod h
    if ( k > d ) k = d;
    pow = zp_malloc (k*sizeof(mpz_t *));
    pow[0] = g;
    for ( i = 1 ; i < k ; i++ ) {
        pow[i] = zp_poly_mod_alloc (h);
        zp_poly_mod_mul (pow[i], pow[i-1], pow[0], h);
    }

    // apply Horner's method using a window size of k
    mpz_init_modp_w (x, h->p);
    w = zp_poly_mod_alloc (h);
    j = k * (d/k);
//printf ("k=%d, d=%d, d/k=%d, j=%d\n", k, d, d/k, j);
    _zp_poly_linear_combo_mod_xn (r, pow, d-j, f+j+1, h->d, h->p, x);  mpz_add_mod (r[0], r[0], f[j], h->p);
//printf ("j=%d r = ", j); zp_poly_print (r, h->d-1); puts ("");
    for ( j -= k ; j >= 0 ; j -= k ) {
        zp_poly_mod_mul (r, r, pow[k-1], h);
//printf ("j=%d r = ", j); zp_poly_print (r, h->d-1); puts ("");
        _zp_poly_linear_combo_mod_xn (w, pow, k-1, f+j+1, h->d, h->p, x);  mpz_add_mod (w[0], w[0], f[j], h->p);
//printf ("j=%d w = ", j); zp_poly_print (w, h->d-1); puts ("");
        zp_poly_mod_add (r, r, w, h);
//printf ("j=%d r = ", j); zp_poly_print (r, h->d-1); puts ("");
    }
    zp_poly_mod_free (w, h);
    mpz_clear (x);
    for ( i = 1 ; i < k ; i++ ) zp_poly_mod_free (pow[i], h);
    zp_free (pow, k*sizeof(mpz_t *));
    if ( wf == f ) zp_poly_free (wf, d);
    if ( wg == g ) zp_poly_mod_free (wg, h);
}

#ifdef ZP_POLY_COMP_TEST
#include <time.h>
 
int main (int argc, char *argv[])
{
    clock_t start, end;
    zp_mod_t pmod;
    zp_poly_mod_t mod;
    mpz_t p, *r, *f, *g, *h, *w, *q, x, y, z;
    int i, k, d, mind, maxd, df, delta, pbits, modonly;
    
    if ( argc < 5 ) { puts ("zp_poly_comp_test mind maxd df pbits [modonly]\n    computes f(g) and f(g) mod h with f of degree d and deg(h) in [mind,maxd], deg(g) = deg(h)-1"); return 0; }
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);
    mind = atoi (argv[1]); maxd = atoi (argv[2]);  df = atoi(argv[3]);  pbits = atoi(argv[4]);
    if ( argc > 5 ) modonly = atoi (argv[5]); else modonly = 0;
    mpz_init2(p,pbits);  zp_random_prime (p,pbits);
    zp_mod_init (pmod, p);
    mpz_init_modp_w (x, p);  mpz_init_modp_w (y, p);  mpz_init_modp_w (z, p);
    
    r = zp_poly_alloc (df*maxd, p);  q = zp_poly_alloc (df*maxd, p);
    f = zp_poly_alloc (df, p);  g = zp_poly_alloc (maxd, p); h = zp_poly_alloc (maxd, p);  w = zp_poly_alloc (maxd, p);
    zp_poly_randomize (f, df, p);
//printf ("f = "); zp_poly_print (f, df); puts ("");
    
    for ( i =  delta = 1 ; i <= mind && i <= maxd-mind ; i++ ) if ( !(mind%i) && !((maxd-mind)%i) ) delta = i;
    if ( ! delta ) delta = 1;  // a really stupid GCD   
    for ( d = mind ; d <= maxd ; d += delta ) {
        zp_poly_randomize (g, d-1, p);
        zp_poly_randomize (h, d-1, p);  mpz_set_one (h[d]);
        if ( modonly ) goto skip;
//printf ("g = "); zp_poly_print (g, d-1); puts ("");
        start = clock();
        zp_poly_compose_horner (r, f, df, g, d-1, pmod);
//printf ("r = "); zp_poly_print (r, df*(d-1)); puts ("");
        end = clock();
        printf ("Computing f(g) with deg(f)=%d and deg(g)=%d took %.3f secs\n", df, d-1, (double)(end-start)/CLOCKS_PER_SEC);
        for ( i = 0 ; i < 10 ; i++ ) {
            zp_randomm (x, p);
            zp_poly_eval (z, g, d-1, x, pmod);
            zp_poly_eval (z, f, df, z, pmod);
            zp_poly_eval (y, r, df*(d-1), x, pmod);
            if ( mpz_cmp (y, z) ) { printf ("Random verification failed!\n"); abort(); }
        }
//printf ("h = "); zp_poly_print (h, d); puts ("");
        zp_poly_div (q, r, df*(d-1), h, d, p);
        zp_poly_rem (r, r, df*(d-1), h, d, q, p);
//printf ("r mod h = "); zp_poly_print (r, d-1); puts ("");
skip:
        zp_poly_mod_init (mod, h, d, p, 0);
        start = clock();
        zp_poly_mod_compose_horner (h, f, df, g, mod);
        end = clock();
        printf ("Computing f(g) mod h with deg(f)=%d and deg(h)=%d took %.3f secs\n", df, d, (double)(end-start)/CLOCKS_PER_SEC);
//printf ("f(g) mod h = "); zp_poly_print(h, d-1); puts ("");
        if ( ! modonly && ! zp_poly_mod_equal (r, h, mod) ) { printf ("zp_poly_mod_compose_horner failed!\n"); abort(); }
        for ( k = 1 ; k < 20 ; k++ ) {
            start = clock();
            zp_poly_mod_compose_window (w, f, df, g, mod, k);
            end = clock();
            printf ("Computing f(g) mod h with window size %d took %.3f secs\n", k, (double)(end-start)/CLOCKS_PER_SEC);
//printf ("k=%d: f(g) mod h = ", k); zp_poly_print(w, d-1); puts ("");
            if ( ! zp_poly_mod_equal (w, h, mod) ) { printf ("zp_poly_mod_compose_window failed!\n"); abort(); }
        }
        zp_poly_mod_clear (mod);
    }
    zp_poly_free (r, df*maxd); zp_poly_free (q, df*maxd);
    zp_poly_free (f, df);  zp_poly_free (g, maxd); zp_poly_free (h, maxd); zp_poly_free (w, maxd);
    mpz_clear (x); mpz_clear (y); mpz_clear (z); mpz_clear (p);
    zp_mod_clear (pmod);
    zp_cleanup();
    zp_mem_report (1);
}

#endif
