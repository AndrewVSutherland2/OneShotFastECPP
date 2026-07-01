#include <assert.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "gmp.h"
#include "ff.h"
#include "ffpoly.h"
#include "ffpolysmall.h"
#include "ffpolyalloc.h"
#include "ffpolybig.h"
#include "polyparse.h"
#include "cstd.h"

/*
    Copyright 2008-2019 Andrew V. Sutherland
    See LICENSE file for license details.
*/


// find a rational point (x0:y0:1) on conic c = c0*x^2 + c1*x*y + c2*x*z + c3*y^2 + c4*y*z + c5*z^2 (monomial lex order)
void ff_poly_conic_rational_point (ff_t *x0, ff_t *y0, ff_t c[6])
{
    register ff_t t;
    ff_t f[3];

    assert ( ! _ff_zero(c[1]) || !_ff_zero(c[3]) || !_ff_zero(c[4]) );
    _ff_set_zero(*x0);  _ff_set(f[0],c[3]);  _ff_set(f[1],c[4]);  _ff_set(f[2],c[5]);
    while ( !ff_poly_find_root (y0, f, ff_poly_degree(f,2)) ) {
        _ff_inc(*x0);
        _ff_mult (t, c[0], *x0); _ff_addto (t, c[2]); _ff_multby (t, *x0); _ff_add (f[0], t, c[5]);
        _ff_mult (t, c[1], *x0); _ff_add (f[1], c[4], t);
        _ff_set (f[2], c[3]);
    }
}

// given conic c, finds rational point (x0:y0:1) and computes cc = c(x+x0,y+y0,z) with (0:0:1) as a rational point (so z^2 coeff cc[5] is zero)
void ff_poly_conic_normalize (ff_t cc[6], ff_t *x0, ff_t *y0, ff_t c[6])
{
    register ff_t s, t;
    
    ff_poly_conic_rational_point (x0, y0, c);
    memmove (cc, c, 6*sizeof(cc[0]));
    _ff_double(t,*x0); _ff_sum_2_mults(s,t,*y0,c[1],c[0]); _ff_addto(cc[2],s); // cc2 = c2 + 2*c0*x0 + c1*y0
    _ff_double(t,*y0); _ff_sum_2_mults(s,*x0,t,c[3],c[1]); _ff_addto(cc[4],s); // cc4 = c4 + 2*c3*y0 + c1*x0
    _ff_set_zero (cc[5]);
}

// c0*x^4 + c1*x^3*y + c2*x^3*z + c3*x^2*y^2 + c4*x^2*y*z + c5*x^2*z^2 + c6*x*y^3 + c7*x*y^2*z + c8*x*y*z^2 + c9*x*z^3
// + c10*y^4 + c11*y^3*z + c12*y^2*z^2 + c13*y*z^3 + c14*z^4 (monomial lex order)

// computes translated plane quartic cc(x,y,z) = c(x+x0*z,y+y0*z,z)
void ff_poly_plane_quartic_translate (ff_t cc[15], ff_t c[15], ff_t x0, ff_t y0)
{
    register ff_t t, w;

    memmove (cc, c, 15*sizeof(cc[0]));
    if ( cc[0] || cc[1] || cc[2] || cc[3] || cc[4] || cc[5] ) {
        // general case

        // cc14 = c14 + c13*y0 + c12*y0^2 + c11*y0^3 + c10*y0^4 + c9*x0 + c8*x0*y0 + c7*x0*y0^2 + c6*x0*y0^3 + c5*x0^2 + c4*x0^2*y0 + c3*x0^2*y0^2 + c2*x0^3 + c1*x0^3*y0 + c0*x0^4
        _ff_mult (t, c[10], y0); _ff_addto (t, c[11]); _ff_multby (t, y0);  _ff_addto (t, c[12]);  _ff_multby (t, y0); _ff_addto (t, c[13]);  _ff_multby (t, y0); _ff_addto (cc[14], t);
        _ff_mult (t, c[6], y0); _ff_addto (t, c[7]); _ff_multby (t, y0);  _ff_addto (t, c[8]);  _ff_multby (t, y0); _ff_addto (t, c[9]);  _ff_multby (t, x0); _ff_addto(cc[14], t);
        _ff_mult (t, c[3], y0); _ff_addto (t, c[4]); _ff_multby (t, y0);  _ff_addto (t, c[5]); _ff_multby(t, x0); _ff_multby(t, x0); _ff_addto(c[14],t);
        _ff_mult (t, c[1], y0); _ff_addto (t, c[2]); _ff_multby (t, x0); _ff_multby (t, x0); _ff_multby (t, x0); _ff_addto(cc[14], t);
        _ff_mult (t, c[0], x0); _ff_square (w, t); _ff_mult (t, w, c[0]);  _ff_addto(cc[14], t);
        // cc13 = c13 + 2*c12*y0 + 3*c11*y0^2 + 4*c10*y0^3 + c8*x0 + 2*c7*x0*y0 + 3*c6*x0*y0^2 + c4*x0^2 + 2*c3*x0^2*y0 + c1*x0^3
        _ff_quadruple (w, c[10]); _ff_mult (t, w, y0); _ff_triple (w, c[11]); _ff_addto (t, w);  _ff_multby (t, y0); _ff_double (w, c[12]); _ff_addto(t, w)  _ff_multby (t, y0);  _ff_addto (cc[13], t);
        _ff_triple (w, c[6]);  _ff_mult (t, w, y0);  _ff_double (w, c[7]); _ff_addto(t, w); _ff_multby (t, y0);  _ff_addto (t, c[8]); _ff_multby (t, x0); _ff_addto (cc[13], t);
        _ff_mult (t, c[1], x0);  _ff_mult (w, c[3], y0); _ff_x2 (w); _ff_addto (t, w); _ff_addto (t, c[4]); _ff_multby (t, x0); _ff_multby(t, x0); _ff_addto (cc[13], t);
        // cc12 = c12 + 3*c11*y0 + 6*c10*y0^2 + c7*x0 + 3*c6*x0*y0 + c3*x0^2
        _ff_triple (w, c[10]); _ff_x2(w);  _ff_mult (t, w, y0);  _ff_triple (w, c[11]);  _ff_addto (t, w);  _ff_multby (t, y0);  _ff_addto (cc[12], t);
        _ff_mult (t, c[6], y0); _ff_addto (t, c[7]); _ff_multby (t, x0); _ff_addto (cc[12], t);
        _ff_square (t, x0); _ff_multby (t, c[3]); _ff_addto (cc[12], t);
        // cc11 = c11 + 4*c10*y0 + c6*x0
        _ff_quadruple (w, c[10]); _ff_mult(t,w,y0); _ff_mult (w, c[6], x0);  _ff_addto (t, w); _ff_addto (cc[11], t);
        // q10 = c10
        // cc9 = c9 + c8*y0 + c7*y0^2 + c6*y0^3 + 2*c5*x0 + 2*c4*x0*y0 + 2*c3*x0*y0^2 + 3*c2*x0^2 + 3*c1*x0^2*y0 + 4*c0*x0^3
        _ff_mult (t, c[6], y0); _ff_addto (t, c[7]); _ff_multby (t, y0); _ff_addto (t, c[8]); _ff_multby (t, y0); _ff_addto (cc[9], t);
        _ff_mult (t, c[3], y0); _ff_addto (t, c[4]); _ff_multby (t, y0); _ff_addto(t, c[5]); _ff_x2(t); _ff_multby (t, x0); _ff_addto (cc[9], t);
        _ff_quadruple (t, c[0]); _ff_multby (t, x0); _ff_triple (w, c[1]); _ff_multby(w, y0); _ff_addto (t, w); _ff_multby (t, x0);  _ff_triple (w, c[2]);  _ff_addto (t, w); _ff_multby (t, x0); _ff_addto (cc[9], t);
        // cc8 = c8 + 2*c7*y0 + 3*c6*y0^2 + 2*c4*x0 + 4*c3*x0*y0 + 3*c1*x0^2
        _ff_triple (t, c[6]); _ff_multby (t, y0); _ff_double (w, c[7]); _ff_addto (t, w); _ff_multby (t, y0); _ff_addto (cc[8], t);
        _ff_double (t, c[3]); _ff_multby (t, y0);  _ff_addto(t, c[4]); _ff_x2(t); _ff_multby (t, x0); _ff_addto (cc[8], t);
        _ff_triple (t, c[1]); _ff_square (w, x0); _ff_multby (t, w); _ff_addto (cc[8], t);
        // cc7 = c7 + 3*c6*y0 + 2*c3*x0,
        _ff_triple (t, c[6]); _ff_multby (t, y0); _ff_addto(cc[7], t);
        _ff_double (t, c[3]); _ff_multby (t, x0); _ff_addto (cc[7], t);
        // cc6 = c6
        // cc5 = c5 + c4*y0 + c3*y0^2 + 3*c2*x0 + 3*c1*x0*y0 + 6*c0*x0^2,
        _ff_mult (t, c[3], y0); _ff_addto (t, c[4]); _ff_multby (t, y0); _ff_addto (cc[5], t);
        _ff_triple (t, c[1]); _ff_multby (t, x0); _ff_multby (t, y0); _ff_addto (cc[5], t);
        _ff_triple (t, c[0]); _ff_double(w,x0); _ff_multby (t, w); _ff_triple (w, c[2]); _ff_addto (t, w); _ff_multby (t, x0); _ff_addto (cc[5], t);
        // cc4 = c4 + 2*c3*y0 + 3*c1*x0
        _ff_mult (t, c[3], y0); _ff_x2(t); _ff_addto (cc[4], t);
        _ff_triple (t, c[1]); _ff_multby (t, x0); _ff_addto (cc[4], t);
        // cc3 = c3
        // cc2 = c2 + c1*y0 + 4*c0*x0,
        _ff_mult (t, c[1], y0); _ff_mult (w, c[0], x0); _ff_x2(w); _ff_x2(w); _ff_addto (t, w); _ff_addto (cc[2], t);
        // cc1 = c1
        // cc0 = c0
    } else {
        // expected special case f(x,y,z) = x*(c6*y^3 + c7*y^2*z + c8*y*z^2 + c9*z^3) + c10*y^4 + c11*y^3*z + c12*y^2*z^2 + c13*y*z^3 + c14*z^4
        // cc14 = c14 + c13*y0 + c12*y0^2 + c11*y0^3 + c10*y0^4 + c9*x0 + c8*x0*y0 + c7*x0*y0^2 + c6*x0*y0^3 + c5*x0^2 + c4*x0^2*y0 + c3*x0^2*y0^2 + c2*x0^3 + c1*x0^3*y0 + c0*x0^4
        _ff_mult (t, c[10], y0); _ff_addto (t, c[11]); _ff_multby (t, y0);  _ff_addto (t, c[12]);  _ff_multby (t, y0); _ff_addto (t, c[13]);  _ff_multby (t, y0); _ff_addto (cc[14], t);
        _ff_mult (t, c[6], y0); _ff_addto (t, c[7]); _ff_multby (t, y0);  _ff_addto (t, c[8]);  _ff_multby (t, y0); _ff_addto (t, c[9]);  _ff_multby (t, x0); _ff_addto(cc[14], t);
        // cc13 = c13 + 2*c12*y0 + 3*c11*y0^2 + 4*c10*y0^3 + c8*x0 + 2*c7*x0*y0 + 3*c6*x0*y0^2 + c4*x0^2 + 2*c3*x0^2*y0 + c1*x0^3
        _ff_quadruple (w, c[10]); _ff_mult (t, w, y0); _ff_triple (w, c[11]); _ff_addto (t, w);  _ff_multby (t, y0); _ff_double (w, c[12]); _ff_addto(t, w)  _ff_multby (t, y0);  _ff_addto (cc[13], t);
        _ff_triple (w, c[6]);  _ff_mult (t, w, y0);  _ff_double (w, c[7]); _ff_addto(t, w); _ff_multby (t, y0);  _ff_addto (t, c[8]); _ff_multby (t, x0); _ff_addto (cc[13], t);
        // cc12 = c12 + 3*c11*y0 + 6*c10*y0^2 + c7*x0 + 3*c6*x0*y0 + c3*x0^2
        _ff_triple (w, c[10]); _ff_x2(w);  _ff_mult (t, w, y0);  _ff_triple (w, c[11]);  _ff_addto (t, w);  _ff_multby (t, y0);  _ff_addto (cc[12], t);
        _ff_mult (t, c[6], y0); _ff_addto (t, c[7]); _ff_multby (t, x0); _ff_addto (cc[12], t);
        // cc11 = c11 + 4*c10*y0 + c6*x0
        _ff_quadruple (w, c[10]); _ff_mult(t,w,y0); _ff_mult (w, c[6], x0);  _ff_addto (t, w); _ff_addto (cc[11], t);
        // q10 = c10
        // cc9 = c9 + c8*y0 + c7*y0^2 + c6*y0^3 + 2*c5*x0 + 2*c4*x0*y0 + 2*c3*x0*y0^2 + 3*c2*x0^2 + 3*c1*x0^2*y0 + 4*c0*x0^3
        _ff_mult (t, c[6], y0); _ff_addto (t, c[7]); _ff_multby (t, y0); _ff_addto (t, c[8]); _ff_multby (t, y0); _ff_addto (cc[9], t);
        // cc8 = c8 + 2*c7*y0 + 3*c6*y0^2 + 2*c4*x0 + 4*c3*x0*y0 + 3*c1*x0^2
        _ff_triple (t, c[6]); _ff_multby (t, y0); _ff_double (w, c[7]); _ff_addto (t, w); _ff_multby (t, y0); _ff_addto (cc[8], t);
        // cc7 = c7 + 3*c6*y0 + 2*c3*x0,
        _ff_triple (t, c[6]); _ff_multby (t, y0); _ff_addto(cc[7], t);
        // cc6 = c6
    }
}


// compute f(t) = q(x(t),y(t),z(t)) where x,y,z are quadratic polys in t (possibly with lc zero), and c(x,y,z) is a plane quartic
int ff_poly_plane_quartic_quadratic_eval (ff_t f[9], ff_t c[15], ff_t x[3], ff_t y[3], ff_t z[3])
{
    ff_t xx[5], xy[5], yy[5], xz[5], yz[5], zz[5], g[9];
    int d = -1;

    // this could be sped up a lot, especially in the expected case where c[0[,...c[5] are zero and x is linear
    ff_poly_square_small (xx, x, 2);  ff_poly_square_small (yy, y, 2);  ff_poly_square_small (zz, z, 2);
    ff_poly_mult_small (xy, x, y, 2); ff_poly_mult_small (xz, x, z, 2); ff_poly_mult_small (yz, y, z, 2);
    if ( c[14] ) { ff_poly_square_small (g, zz, 4);  ff_poly_scalar_mult (f, &d, c[14], g, 8); }
    if ( c[13] ) { ff_poly_mult_small (g, yz, zz, 4);  ff_poly_scalar_mult (g, 0, c[13], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[12] ) { ff_poly_mult_small (g, yy, zz, 4);  ff_poly_scalar_mult (g, 0, c[12], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[11] ) { ff_poly_mult_small (g, yy, yz, 4);  ff_poly_scalar_mult (g, 0, c[11], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[10] ) { ff_poly_square_small (g, yy, 4);    ff_poly_scalar_mult (g, 0, c[10], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[9] ) { ff_poly_mult_small (g, xz, zz, 4);  ff_poly_scalar_mult (g, 0, c[9], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[8] ) { ff_poly_mult_small (g, xy, zz, 4);  ff_poly_scalar_mult (g, 0, c[8], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[7] ) { ff_poly_mult_small (g, xy, yz, 4);  ff_poly_scalar_mult (g, 0, c[7], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[6] ) { ff_poly_mult_small (g, xy, yy, 4);  ff_poly_scalar_mult (g, 0, c[6], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[5] ) { ff_poly_mult_small (g, xx, zz, 4);  ff_poly_scalar_mult (g, 0, c[5], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[4] ) { ff_poly_mult_small (g, xx, yz, 4);  ff_poly_scalar_mult (g, 0, c[4], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[3] ) { ff_poly_mult_small (g, xx, yy, 4);  ff_poly_scalar_mult (g, 0, c[3], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[2] ) { ff_poly_mult_small (g, xx, xz, 4);  ff_poly_scalar_mult (g, 0, c[2], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[1] ) { ff_poly_mult_small (g, xx, xy, 4);  ff_poly_scalar_mult (g, 0, c[1], g, 8);  ff_poly_addto (f, &d, g, 8); }
    if ( c[0] ) { ff_poly_square_small (g, xx, 4);  ff_poly_scalar_mult (g, 0, c[0], g, 8);  ff_poly_addto (f, &d, g, 8); }
    return d;
}

int ff_poly_hyperelliptic_curve_from_conic_cover (ff_t g[9], ff_t c[6], ff_t f[15])
{
    ff_t x0, y0, x[3], y[3], z[3], cc[6], ff[15];
    
    ff_poly_conic_normalize (cc, &x0, &y0, c);
    ff_poly_plane_quartic_translate (ff, f, x0, y0);
    // rational parameterization of conic [c0,c1,c2,c3,c4,0] is (c4*t+c2 : c4*t^2+c2*t : -c3*t^2-c1*t-c0)
    _ff_set (x[0], cc[2]);  _ff_set (x[1], cc[4]);  _ff_set_zero (x[2]);
    _ff_set_zero (y[0]);  _ff_set (y[1], cc[2]);  _ff_set (y[2], cc[4]);
    _ff_neg (z[0], cc[0]);  _ff_neg (z[1], cc[1]);  _ff_neg (z[2], cc[3]);
    return ff_poly_plane_quartic_quadratic_eval (g, ff, x, y, z);
}

// f and W can coincide (W = [a1,a2,a3,a4,a6])
void ff_poly_med_weierstrass (ff_t f[4], ff_t W[5])
{
    register ff_t t1, t2;

//printf ("[%ld,%ld,%ld,%ld,%ld]\n", _ff_get_ui(W[0]), _ff_get_ui(W[1]), _ff_get_ui(W[2]), _ff_get_ui(W[3]), _ff_get_ui(W[4])); 
    _ff_set_one (f[3]);
    _ff_mult (t1, W[0], _ff_half);                  // t1 = a1/2
    _ff_square (t2, t1);                            // c2 = a1^2/4
    _ff_addto (t2, W[1]);                       // f2 = t2 = a2+a1^2/4
    _ff_mult(t1,t1,W[2]);                       // c4 = a1a3/2
    _ff_add(f[1],t1,W[3]);                      // f1 = a4 + a1a3/2
    _ff_mult(t1,W[2], _ff_half);                    // t1 = a3/2
    _ff_square(t1,t1);                          // c6 = a3^2/4
    _ff_add(f[0],t1,W[4]);                      // f0 = a6 + a3^2/4
    _ff_set (f[2],t2);                          // f2 = t2
//printf ("y^2 = x^3 + %ld*x^2 + %ld*x + %ld\n", _ff_get_ui(f[2]), _ff_get_ui(f[1]), _ff_get_ui(f[0]));
    // 3M+2S
}   

// f and W can coincide(W = [a1,a2,a3,a4,a6])
void ff_poly_short_weierstrass (ff_t f[4], ff_t W[5])
{
    register ff_t t1, t2, t3, c2, c4, c6;
    
    _ff_set_one (f[3]);
    _ff_set_zero (f[2]);
    
//printf ("[%ld,%ld,%ld,%ld,%ld]\n", _ff_get_ui(W[0]), _ff_get_ui(W[1]), _ff_get_ui(W[2]), _ff_get_ui(W[3]), _ff_get_ui(W[4]));
    _ff_mult (t1, W[0], _ff_half);                  // t1 = a1/2
    _ff_square (c2, t1);                            // c2 = a1^2/4
    _ff_addto (c2, W[1]);                       // c2 = a2+a1^2/4 is the new a2
    _ff_mult(c4,t1,W[2]);                       // c4 = a1a3/2
    _ff_addto(c4,W[3]);                         // c4 = a4 + a1a3/2 is the new a4
    _ff_mult(t2,W[2], _ff_half);                    // t2 = a3/2
    _ff_square(c6,t2);                          // c6 = a3^2/4
    _ff_addto(c6,W[4]);                         // c6 = a6 + a3^2/4 is the new a6 
//printf ("[%ld,%ld,%ld]\n", _ff_get_ui(c2), _ff_get_ui(c4), _ff_get_ui(c6));
    _ff_mult(t1,c2,_ff_third);                      // t1 = c2/3
    _ff_mult(t2,t1,c2);                         // t2 = c2^2/3
    _ff_sub (f[1], c4, t2);                     // f[1] = c4 - c2^2/3
    _ff_mult(t2,t2,t1);                         // t2 = c2^3/9
    _ff_mult(t3,t2,_ff_third);                      // t3 = c2^3/27
    _ff_add(t2,t3,t3);                          // t2 = 2/27*c2^3
    _ff_mult(t3,t1,c4);                         // t3 = c2c4/3
    _ff_subfrom(t2,t3);                         // t2 = 2/27*c2^3 - c2c4/3
    _ff_add (f[0], c6, t2);                     // f[0] = c6 + 2/27*c2^3 - c2c4/3
//printf ("x^3+%ld*x+%ld\n", _ff_get_ui(f[1]), _ff_get_ui(f[0]));
    // 8M+2S
}   


int ff_poly_parse (ff_t f[], int maxd, char *expr)
{
    int i, d;
    
    d = ui_poly_parse_mod_p (f, maxd, expr, _ff_p);
    for ( i = 0 ; i <= d ; i++ ) _ff_set_ui (f[i], f[i]);
    for ( ; i <= maxd ; i++ ) _ff_set_zero (f[i]);      // pad out with zeroes 
    return d;
}


void ff_poly_print (ff_t f[], int d_f)
{
    char buf[65536];

    ff_poly_sprint (buf, f, d_f);
    out_printf ("%s\n", buf);
}

int ff_poly_sprint (char *s, ff_t f[], int d_f)
{
    char *t;
    int i;
    
    if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
    t = s;
    if ( d_f >= 2 ) {
        if ( _ff_one(f[d_f]) ) t +=  gmp_sprintf (t, "[x^%d", d_f);
        else t += gmp_sprintf (t, "[%lu*x^%d", _ff_get_ui(f[d_f]), d_f);
    } else if ( d_f == 1 ) {
        if ( _ff_one(f[d_f]) ) t +=  gmp_sprintf (t, "[x");
        else t += gmp_sprintf (t, "[%lu*x", _ff_get_ui(f[d_f]));
    } else {
        t += gmp_sprintf (t, "[%lu", _ff_get_ui(f[d_f]));
    }
    for ( i = d_f-1 ; i >= 0 ; i-- ) {
        if ( _ff_nonzero(f[i]) ) {
            if ( i >= 2 ) {
                if ( _ff_one(f[i]) ) t += gmp_sprintf (t, " + x^%d", i);
                else t += gmp_sprintf (t, " + %lu*x^%d", _ff_get_ui(f[i]), i);
            } else if ( i == 1 ) {
                if ( _ff_one(f[i]) ) t += gmp_sprintf (t, " + x");
                else t += gmp_sprintf (t, " + %lu*x", _ff_get_ui(f[i]));
            } else {
                t += gmp_sprintf (t, " + %lu", _ff_get_ui(f[i]));
            }
        }
    }
    *t++ = ']';
    *t= '\0';
    return t-s;
}

/*
    Computes the quadratic twist of y^2=f(x) by a non-residue in F_p
    f and g may coincide
*/
void ff_poly_twist (ff_t g[], ff_t f[], int d)
{
    ff_t a, b;
    int i;

    ff_nonresidue(&b);
    if ( (d&0x1) ) {
        // in the odd degree case we can keep the leading coefficient - means that g is monic if f is.
        _ff_set (g[d], f[d]);                   // leading coefficient unchanged - means twist is monic if f is
        _ff_set (a,b);
        for ( i = d-1 ; i >= 0 ; i-- ) {
            ff_mult (g[i], b, f[i]);                // g[i] = a^{d-i}f[i] for i = d-1 down to 0
            ff_mult (b, b, a);
        }
    } else {
        for ( i = 0 ; i <= d ; i++ ) _ff_mult(g[i],b,f[i]);
    }
    return;
}


// Naive Euclidean division of a by b, computes q and r such that a = qb+r with deg(r) < deg(b).  Overlap is ok.  Does not require any inversions if b is monic
// pd_q, and r optional (if r is specified pd_r must be specified)
void ff_poly_div (ff_t q[], int *pd_q, ff_t r[], int *pd_r, ff_t a[], int d_a, ff_t b[], int d_b)
{
    ff_t *aa;
    register ff_t s, t, x;
    int i, d_q, d_s;
    
    assert ( d_b >= 0 );
    if ( d_a < d_b ) { if ( pd_q ) *pd_q = -1; if ( r ) ff_poly_copy (r, pd_r, a, d_a); return; }

    a = ff_poly_copy (aa = ff_poly_stack_alloc(d_a), 0, a, d_a);
    if ( q == b ) b = ff_poly_copy (ff_poly_stack_alloc(d_b), 0, b, d_b);
    
    d_q = -1;
    // duplicate code for non-monic/monic cases to avoid extra branch
    if ( _ff_one(b[d_b]) ) {
        while ( d_a >= d_b ) {
            _ff_set (s, a[d_a]);
            d_s = d_a - d_b;
            if ( q ) {
                if ( d_s > d_q ) {
                    for ( i = d_q+1 ; i < d_s ; i++ ) _ff_set_zero (q[i]);
                    d_q = d_s;
                    _ff_set (q[d_q], s);
                } else {
                    _ff_addto (q[d_s], s);
                }
            }
            for ( i = 0 ; i < d_b ; i++ ) { _ff_mult (t, s, b[i]); _ff_subfrom (a[i+d_s], t); }
            _ff_set_zero(a[d_a]);
            d_a = ff_poly_degree(a,d_a-1);
        }
    } else {
        _ff_invert (x, b[d_b]);
        while ( d_a >= d_b ) {
            _ff_mult (s, x, a[d_a]);
            d_s = d_a - d_b;
            if ( q ) {
                if ( d_s > d_q ) {
                    for ( i = d_q+1 ; i < d_s ; i++ ) _ff_set_zero (q[i]);
                    d_q = d_s;
                    _ff_set (q[d_q], s);
                } else {
                    _ff_addto (q[d_s], s);
                }
            }
            for ( i = 0 ; i < d_b ; i++ ) { _ff_mult (t, s, b[i]); _ff_subfrom (a[i+d_s], t); }
            _ff_set_zero(a[d_a]);
            d_a = ff_poly_degree(a,d_a-1);
        }
    }
    if ( pd_q ) *pd_q = d_q;
    if ( r ) ff_poly_copy (r, pd_r, a, d_a);
    else if ( pd_r ) *pd_r = d_a;
    ff_poly_stack_pop(aa);
}

// TODO optimize this code
void ff_poly_div_exact (ff_t q[], int *pd_q, ff_t a[], int d_a, ff_t b[], int d_b)
{
    ff_t *aa;
    register ff_t s, t, x;
    int i, d_q, d_s;
    
    if ( d_a < 0 ) { if ( pd_q ) *pd_q = d_a; return; } // a=0 is divisible by everything
    assert ( d_b >= 0 && d_a >= d_b );

    a = ff_poly_copy (aa = ff_poly_stack_alloc (d_a), 0, a, d_a);
    if ( q == b ) b = ff_poly_copy (ff_poly_stack_alloc(d_b), 0, b, d_b);
    
    d_q = -1;
    // duplicate code for non-monic/monic cases to avoid extra branch
    if ( _ff_one(b[d_b]) ) {
        while ( d_a >= d_b ) {
            _ff_set (s, a[d_a]);
            d_s = d_a - d_b;
            if ( q ) {
                if ( d_s > d_q ) {
                    for ( i = d_q+1 ; i < d_s ; i++ ) _ff_set_zero (q[i]);
                    d_q = d_s;
                    _ff_set (q[d_q], s);
                } else {
                    _ff_addto (q[d_s], s);
                }
            }
            for ( i = 0 ; i < d_b ; i++ ) { _ff_mult (t, s, b[i]); _ff_subfrom (a[i+d_s], t); }
            _ff_set_zero(a[d_a]);
            d_a = ff_poly_degree(a,d_a-1);
        }
    } else {
        _ff_invert (x, b[d_b]);
        while ( d_a >= d_b ) {
            _ff_mult (s, x, a[d_a]);
            d_s = d_a - d_b;
            if ( q ) {
                if ( d_s > d_q ) {
                    for ( i = d_q+1 ; i < d_s ; i++ ) _ff_set_zero (q[i]);
                    d_q = d_s;
                    _ff_set (q[d_q], s);
                } else {
                    _ff_addto (q[d_s], s);
                }
            }
            for ( i = 0 ; i < d_b ; i++ ) { _ff_mult (t, s, b[i]); _ff_subfrom (a[i+d_s], t); }
            _ff_set_zero(a[d_a]);
            d_a = ff_poly_degree(a,d_a-1);
        }
    }
    if ( pd_q ) *pd_q = d_q;
//    assert ( d_a < 0 );
    ff_poly_stack_pop(aa);
}

// computes h=af-bg of degree less than max(d_f,d_g) (f and g must both be non-zero).   note that h will not typically be monic, even if f and g are
// assuming d_f >= d_g,  h will be set to lc(g)*f-lc(f)*x^(d_f - d_g)*g. In this case gcd(f,h) = gcd(f,g), but gcd(g,h) = x^k*gcd(f,g) for some k <= d_f-d_g (typically k=0, but not always!)
static inline void ff_poly_reduce (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
    register ff_t t0,t1;
    register int i,j;
    
    assert ( d_f >= 0 && d_g >= 0 );
    if ( ! d_f || ! d_g ) { *pd_h = -1; return; }

    // avoid unnecessary copying, but allow overlap of h with f or g
    if ( d_f > d_g ) {
        _ff_neg(t0,f[d_f]);
        _ff_set(t1,g[d_g]);
        j = d_f-d_g;
        for ( i = d_f-1 ; i >= j ; i-- ) _ff_sum_2_mults(h[i],t0,t1,f[i],g[i-j]);
        for ( ; i >= 0 ; i-- ) ff_mult(h[i],f[i],t1);
        for ( i = d_f-1 ; i >= 0 ; i-- ) if ( ! _ff_zero(h[i]) ) break;
        *pd_h = i;
    } else {
        _ff_set(t0,f[d_f]);
        _ff_neg(t1,g[d_g]);     
        j = d_g-d_f;
        for ( i = d_g-1 ; i >= j ; i-- ) _ff_sum_2_mults(h[i],t0,t1,f[i-j],g[i]);
        for ( ; i >= 0 ; i-- ) ff_mult(h[i],g[i],t0);
        for ( i = d_g-1 ; i >= 0 ; i-- ) if ( ! _ff_zero(h[i]) ) break;
        *pd_h = i;
    }
}

// given nonzero polys s and t, reduces one modulo the other until one of them is zero
static inline void _ff_poly_gcd_reduce (ff_t s[], int *pd_s, ff_t t[], int *pd_t)
{
    while ( *pd_s >= 0 && *pd_t >= 0 ) {
        if ( *pd_s < *pd_t ) {
            ff_poly_reduce (t, pd_t, s, *pd_s, t, *pd_t);
        } else {
            ff_poly_reduce (s, pd_s, s, *pd_s, t, *pd_t);
        }
    }
}

// overlap is ok
int ff_poly_gcd_reduce (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b)
{
    ff_t *s, *t, *x;
    int d_s, d_t, d_x, k;

    if ( d_b > d_a ) { _swap (a, b, s);  _swap (d_a, d_b, d_s); }
    if ( d_b < 0 ) { assert (d_a >= 0 ); ff_poly_copy(g,pd_g, a,d_a); return d_a; }
    
    // remove any factors of x from the gcd
    for ( k = 0 ; _ff_zero(a[k]) && _ff_zero(b[k]) ; k++ );
    a += k;  d_a -= k;  b += k;  d_b -= k;
    if ( ! d_a ) {_ff_set_one(g[k]);  *pd_g = k;  for ( k-- ; k >= 0 ; k-- ) _ff_set_zero(g[k]);  return *pd_g; }
    
    // reduce into s and t rather than copying
    s = ff_poly_stack_alloc (d_a);  t = ff_poly_stack_alloc (d_a);
    ff_poly_reduce (s, &d_s, a, d_a, b, d_b);
    if ( d_s < 0 ) { x = b; d_x = d_b; goto done; }
    ff_poly_reduce (t, &d_t, b, d_b, s, d_s);
    _ff_poly_gcd_reduce (s, &d_s, t, &d_t);
    if ( d_s < 0 ) { x = t; d_x=d_t; } else { x = s; d_x = d_s; }
done:
    // remove any factors of x that may have been introduced in the process of reduction, we know they are not in the gcd, since we removed all factors of x above
    while ( d_x >= 0 && _ff_zero(x[0]) ) { x++; d_x--; };
    ff_poly_copy (g+k,0,x,d_x);  *pd_g = d_x + k;
    for ( k-- ; k >= 0 ; k-- ) _ff_set_zero (g[k]);
    if ( ! *pd_g ) _ff_set_one(g[0]);       // make g monic in the trivial case
    ff_poly_stack_pop(s);
    return *pd_g;
}

// note that g is not made monic unless it is trivial (N.B. even if a and b are monic, g might not be)
int ff_poly_gcd (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b)
{
    ff_t *s, *t, *x;
    int d_s, d_t;

    if ( d_b > d_a ) { _swap (a, b, x);  _swap (d_a, d_b, d_s); }
    if ( d_b < 0 ) { ff_poly_copy (g, pd_g, a, d_a); return d_a; }
    
    s = ff_poly_stack_alloc (d_b);  t = ff_poly_stack_alloc (d_b);
    ff_poly_div (0, 0, s, &d_s, a, d_a, b, d_b);
    if ( d_s < 0 ) { ff_poly_copy (g, pd_g, b, d_b); goto done; }
    d_t = 0;
    ff_poly_div (0, 0, t, &d_t, b, d_b, s, d_s);
    while ( d_s >= 0 && d_t >= 0 ) {
        if ( d_s > d_t ) {
            ff_poly_div (0, 0, s, &d_s, s, d_s, t, d_t);
        } else {
            ff_poly_div (0, 0, t, &d_t, t, d_t, s, d_s);
        }
    }
    if ( d_s < 0 ) ff_poly_copy (g, pd_g, t, d_t); else ff_poly_copy (g, pd_g, s, d_s);
done:
    ff_poly_stack_pop(s);
    return *pd_g;
}


// [CCANT] algorithm 3.2.2, makes output monic.  None of the polynomials can overlap.  g may be null.  degree of gcd is returned in any case
// computes u, v, and g such that g=gcd(a,b) = u*a+v*b
int ff_poly_gcdext (ff_t g[], int *pd_g, ff_t u[], int *pd_u, ff_t v[], int *pd_v, ff_t a[], int d_a, ff_t b[], int d_b)
{
    ff_t *f, *q, *r, *t, *t2, *v1, *v3, w0, w1;
    int d_f, d_q, d_r, d_t, d_t2, d_v1, d_v3;
    int n;

    // GCd(0,b) = b = 0*a+1*b
    if ( d_a < 0 ) { if ( g ) ff_poly_copy (g, pd_g, b, d_b); *pd_u = -1;  _ff_set_one(v[0]); *pd_v = 0; return d_b; }
        
    // GCD(a,0) = a = 1*a+0*b
    if ( d_b < 0 ) { if ( g ) ff_poly_copy (g, pd_g, a, d_a); *pd_u= 0;  _ff_set_one(u[0]); *pd_v = -1; return d_a; }
        
    // GCD(a,c) = 1 = 0*a+1/c*c for c a nonzero constant
    if ( ! d_b ) { if ( g ) { *pd_g = 0;  _ff_set_one(g[0]); } *pd_u = -1; _ff_invert(v[0], b[0]); *pd_v = 0; return 0; }
    
    // GCD(c,b) = 1 = 1/c*c+0*b for c a nonzero constant
    if ( ! d_a ) { if ( g ) { *pd_g = 0;  _ff_set_one(g[0]); } _ff_invert (u[0], a[0]); *pd_u = 0; *pd_v = -1; return 0; }
    
    // for simplicity just set n to the max coeffs needed by any of the polys we are using
    n = ( d_a > d_b ? d_a : d_b );
    f = ff_poly_stack_alloc(n); q = ff_poly_stack_alloc(n); r = ff_poly_stack_alloc(n);
    t = ff_poly_stack_alloc(2*n); t2 = ff_poly_stack_alloc(2*n); v1 = ff_poly_stack_alloc(2*n); v3 = ff_poly_stack_alloc(2*n);

    *pd_u = 0;      _ff_set_one(u[0]);              // u = 1
    ff_poly_copy (f, &d_f, a, d_a);                 // f = a
    d_v1 = -1;                              // v1 = 0
    ff_poly_copy (v3, &d_v3, b, d_b);               // v3 = b   

    while ( d_v3 >= 0 ) {
        ff_poly_div (q, &d_q, r, &d_r, f, d_f, v3, d_v3);
        ff_poly_mult (t2, &d_t2, v1, d_v1, q, d_q);
        ff_poly_sub (t, &d_t, u, *pd_u, t2, d_t2);      // t = u-v1q
        // this is a lot of copying, but it means we only use one multiplication
        ff_poly_copy (u, pd_u,  v1, d_v1);          // u = v1
        ff_poly_copy (v1, &d_v1, t, d_t);           // v1 = t
        ff_poly_copy (f, &d_f, v3, d_v3);           // f = v3
        ff_poly_copy (v3, &d_v3, r, d_r);           // v3 = r
    }
    ff_poly_mult (v3, &d_v3, a, d_a, u, *pd_u);
    ff_poly_sub (v1, &d_v1, f, d_f, v3, d_v3);
    ff_poly_div_exact (v, pd_v, v1, d_v1, b, d_b);   // v = (f-au)/b
    // make GCD monic
    if ( ! _ff_one (f[d_f]) ) {
        _ff_set (w1, f[d_f]);
        ff_poly_monic (f, &d_f, f, d_f);
        _ff_invert (w0, w1);
        ff_poly_scalar_mult (u, pd_u, w0, u, *pd_u);
        ff_poly_scalar_mult (v, pd_v, w0, v, *pd_v);
    }
    if ( g ) ff_poly_copy (g, pd_g, f, d_f);
    ff_poly_stack_pop(f);
    return d_f;
}

int ff_poly_inv_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
    ff_t *v;
    int d_v;
    int n;
    
    v = ff_poly_stack_alloc (d_f);
    n = ff_poly_gcdext (0, 0, h,  pd_h, v, &d_v, f, d_f, g, d_g);
    ff_poly_stack_pop (v);
    return ( n ? 0 : 1);
}


// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e, ff_t f[], int d_f)
{
    ff_poly_modulus_t mod;
    
    assert (d_a < d_f);
    ff_poly_mod_setup (mod, f, d_f);
    ff_poly_pow_modulus (g, pd_g, a, d_a, e, mod);
    ff_poly_mod_clear(mod);
}

// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_modulus (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e,  ff_poly_modulus_t mod)
{
    ff_t *h;
    int d_h;
    int t;

    assert (d_a < mod->d);
    if ( d_a < 0 || !e ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }

    h = ff_poly_stack_alloc (2*mod->d-2);
    ff_poly_copy (h, &d_h, a, d_a);
    t = ui_len (e) - 2;
    while ( t >= 0 ) {              // TODO: use a faster exponentiation algorithm here
        ff_poly_square (h, h, d_h);
        ff_poly_mod_reduce (h, &d_h, h, 2*d_h, mod);
        if ( e&(1UL<<t) ) {
            ff_poly_mult (h, &d_h, h, d_h, a, d_a);
            ff_poly_mod_reduce (h, &d_h, h, d_h, mod);
        }
        t--;
    }
    ff_poly_copy (g, pd_g, h, d_h);
    ff_poly_stack_pop (h);
}

// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_mod_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_t f[], int d_f)
{
    ff_poly_modulus_t mod;
    
    assert (d_a < d_f);
    ff_poly_mod_setup (mod, f, d_f);
    ff_poly_pow_modulus_mpz (g, pd_g, a, d_a, e, mod);
    ff_poly_mod_clear(mod);
}
    

// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_modulus_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e,  ff_poly_modulus_t mod)
{
    ff_t *h;
    int d_h;
    int t;

    assert (d_a < mod->d);
    if ( d_a < 0 || ! mpz_sgn (e) ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }
    
    h = ff_poly_stack_alloc (2*mod->d-2);
    ff_poly_copy (h, &d_h, a, d_a);
    t = mpz_sizeinbase (e, 2) - 2;
    while ( t >= 0 ) {                  // TODO: use a faster exponentiation algorithm here
        ff_poly_square (h, h, d_h);
        ff_poly_mod_reduce (h, &d_h, h, 2*d_h, mod);
        if ( mpz_tstbit (e, t) ) {
            ff_poly_mult (h, &d_h, h, d_h, a, d_a);
            ff_poly_mod_reduce (h, &d_h, h, d_h, mod);
        }
        t--;
    }
    ff_poly_copy (g, pd_g, h, d_h);
    ff_poly_stack_pop (h);
}


// attempts to compute h = sqrt(f) mod g using slow Tonelli-Shanks, where g is assumed to be irreducible and deg(f) < deg(g)
// h and f may overlap
int ff_poly_sqrt_modulus (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_poly_modulus_t mod)
{
    ff_t *a, *ai, *g;
    int d_g, d_a, d_ai;
    mpz_t m, n, e, w;
    int i, s;

//printf ("computing sqrt of "); ff_poly_print (f, d_f); printf (" mod "); ff_poly_print (mod->g, mod->d);
    if ( d_f < 0 ) { *pd_h = -1; return 1; }
    mpz_init (m);  mpz_init (n);
    mpz_ui_pow_ui (m, _ff_p, mod->d);  mpz_sub_ui (m, m, 1);        // m = p^d-1 is the order of the multiplicative group of Fp[x]/(mod)
    mpz_div_2exp (n, m, 1);                                 // n = (p^d-1)/2
    
    // Determine whether f has a square root or not by computing f^((p^d-1)/2)
    g = ff_poly_stack_alloc (2*mod->d-2);
    ff_poly_pow_modulus_mpz (g, &d_g, f, d_f, n, mod);
    if ( ! ff_poly_is_one (g, d_g) ) { if ( h == f ) ff_poly_stack_pop (g); mpz_clear(m);  mpz_clear(n); return 0; }
    
    // Get a random non-residue a
    a = ff_poly_stack_alloc (mod->d-1);  ai = ff_poly_stack_alloc (mod->d-1);
    d_a = mod->d-1;
    do {
        ff_poly_randomize (a, d_a);
        ff_poly_pow_modulus_mpz (g, &d_g, a, d_a, n, mod);
    } while ( ff_poly_is_one (g, d_g) );
    if ( ! ff_poly_inv_modulus (ai, &d_ai, a, d_a, mod) ) { printf ("Unable to compute inverse in ff_poly_mod_sqrt, modulus poly not irreducible!\n");  ff_poly_print (mod->g, mod->d); abort(); }
    
    for ( s = 1 ; !mpz_tstbit (n,0) ; s++ ) mpz_div_2exp (n,n,1);       // m = n*2^s with n odd

    mpz_init (e);  mpz_init (w);
    for ( i = 2 ; i <= s ; i++ ) {                              // this loop computes e such that f=a^e mod the 2-Sylow
        ff_poly_pow_modulus_mpz (g, &d_g, ai, d_ai, e, mod);
        ff_poly_mult (g, &d_g, g, d_g, f, d_f);
        ff_poly_mod_reduce (g, &d_g, g, d_g, mod);
        mpz_div_2exp (w, m, i);                             // w = (p^d-1)/2^i
        ff_poly_pow_modulus_mpz (g, &d_g, g, d_g, w, mod);
        if ( ! ff_poly_is_one (g, d_g) ) mpz_setbit (e, i-1);
    }
    ff_poly_pow_modulus_mpz (g, &d_g, ai, d_ai, e, mod);
    ff_poly_mult (g, &d_g, g, d_g, f, d_f);                     // g = f*a^-e
    ff_poly_mod_reduce (g, &d_g, g, d_g, mod);
    mpz_add_ui (n, n, 1);  mpz_div_2exp (n, n, 1);
    ff_poly_pow_modulus_mpz (g, &d_g, g, d_g, n, mod);          // g = f*a^(-e(n+1)/2) is the square-root of the odd part of f
    mpz_div_2exp (e, e, 1);
    ff_poly_pow_modulus_mpz (ai, &d_ai, a, d_a, e, mod);            // ai = a^(e/2) is the square root of the even part of f
    ff_poly_mult (g, &d_g, g, d_g, ai, d_ai);                       // g = f*a^(-e(n+1)/2)*a^(e/2) = f*a^(-en) is the square root of f
    ff_poly_mod_reduce (h, pd_h, g, d_g, mod);
    ff_poly_mult (g, &d_g, h, *pd_h, h, *pd_h);                 // verify that h^2 = f
    ff_poly_mod_reduce (g, &d_g, g, d_g, mod);
    assert ( ff_poly_equal (g, d_g, f, d_f) );
    mpz_clear (m); mpz_clear (n); mpz_clear (e); mpz_clear (w);
    ff_poly_stack_pop(g);
    return 1;
}

int ff_poly_sqrt_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
    ff_poly_modulus_t mod;
    int ret;
    
    ff_poly_mod_setup (mod, g, d_g);
    ret = ff_poly_sqrt_modulus (h, pd_h, f, d_f, mod);
    ff_poly_mod_clear (mod);
    return ret;
}


// computes x^n mod f for monic f (should be depressed if degree is small)
void ff_poly_xn_mod (ff_t g[], int *pd_g, unsigned long e, ff_t f[], int d_f)
{
    register int i;
    
    assert ( d_f && _ff_one(f[d_f]) );
    if ( d_f == 1 ) { _ff_neg (g[0], f[0]); ff_exp_ui (g, g, e); return; }
    
    if ( ! e ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }
    
    if ( e < d_f ) { *pd_g = e;  _ff_set_one(g[e]);  for ( i = e-1 ; i >= 0 ; i-- ) _ff_set_zero(g[i]);  return; }

    if ( d_f <= FF_POLY_XNMOD_SMALL_DEGREE && _ff_zero(f[d_f-1]) ) {
        ff_t nf[FF_POLY_XNMOD_SMALL_DEGREE+1];

        assert (_ff_zero(f[d_f-1]));    // for small degrees we require the poly to be depressed 
        // for small mod poly operations all assume the modulus is of the form x^n  - c_{n-2}x^{n-2} - ... - c_1*x - c_0, so negate the non-leading coefficients (leading one is implicit)
        ff_poly_neg(nf,0,f,d_f-2);
        switch (d_f) {
        case 2: ff_exp_ui (g,nf,e>>1); if ( e&1 ) { _ff_set(g[1],g[0]); _ff_set_zero(g[0]); } else { _ff_set_zero(g[1]); } break;
        case 3: ff_poly_xn_mod_d3 (g,e,nf); break;
        case 4: ff_poly_xn_mod_d4 (g,e,nf); break;
        case 5: ff_poly_xn_mod_d5 (g,e,nf); break;
        case 6: ff_poly_xn_mod_d6 (g,e,nf); break;
        case 7: ff_poly_xn_mod_d7 (g,e,nf); break;
        case 8: ff_poly_xn_mod_d8 (g,e,nf); break;
        case 10: ff_poly_xn_mod_d10 (g,e,nf); break;
        case 11: ff_poly_xn_mod_d11 (g,e,nf); break;
        case 12: ff_poly_xn_mod_d12 (g,e,nf); break;
        case 13: ff_poly_xn_mod_d13 (g,e,nf); break;
        case 15: ff_poly_xn_mod_d15 (g,e,nf); break;
        case 17: ff_poly_xn_mod_d17 (g,e,nf); break;
        case 19: ff_poly_xn_mod_d19 (g,e,nf); break;
        case 23: ff_poly_xn_mod_d23 (g,e,nf); break;
        case 29: ff_poly_xn_mod_d29 (g,e,nf); break;
        case 31: ff_poly_xn_mod_d31 (g,e,nf); break;
        default: ff_poly_xn_mod_small (g,e,nf,d_f);
        }
        *pd_g = ff_poly_degree(g,d_f-1);
    } else {
        ff_poly_modulus_t mod;
        ff_poly_mod_setup (mod, f, d_f);
        ff_poly_xn_modulus (g, pd_g, e,mod);
        ff_poly_mod_clear(mod);
    }
}

// computes g(x) = x^n mod mod->g
void ff_poly_xn_modulus  (ff_t g[], int *pd_g, unsigned long e, ff_poly_modulus_t mod)
{
    register ff_t *w;
    register int t, d_w;
    int d_g;
    
    if ( mod->d <= FF_POLY_XNMOD_SMALL_DEGREE && _ff_zero(mod->g[mod->d-1]) ) { ff_poly_xn_mod (g, pd_g, e, mod->g, mod->d); return; }
    if ( e < mod->d ) { *pd_g = e;  _ff_set_one(g[e]);  for ( d_g=e-1 ; d_g >= 0 ; d_g-- ) _ff_set_zero(g[d_g]);  return; }

    _ff_set_one(g[1]);  _ff_set_zero(g[0]);  d_g = 1;  w = mod->w1;
    t = ui_len(e)-2;
    while ( t >= 0 ) {
        if ( e&(1UL<<t) ) {
            ff_poly_square (w+1,g,d_g);   _ff_set_zero(w[0]);
            d_w = ff_poly_degree(w,2*d_g+1);            
        } else {
            ff_poly_square (w,g,d_g);   d_w = ff_poly_degree(w,2*d_g);
        }
        ff_poly_mod_reduce (g,&d_g,w,d_w,mod);
        t--;
    }
    *pd_g = d_g;
}

// computes g(x) = x^n mod mod->g
void ff_poly_xn_modulus_mpz (ff_t g[], int *pd_g, mpz_t e, ff_poly_modulus_t mod)
{
    register ff_t *w;
    register int t, d_w;
    int d_g;
    
    t = mpz_sizeinbase(e,2);
    if ( t <= 64 && mod->d <= FF_POLY_XNMOD_SMALL_DEGREE ) { ff_poly_xn_mod (g, pd_g, mpz_get_ui(e), mod->g, mod->d); return; }
    _ff_set_one(g[1]);  _ff_set_zero(g[0]);  d_g = 1;  w = mod->w1;
    t -= 2;
    while ( t >= 0 ) {
        if ( mpz_tstbit(e,t) ) {
            ff_poly_square (w+1,g,d_g);   _ff_set_zero(w[0]);
            d_w = ff_poly_degree(w,2*d_g+1);            
        } else {
            ff_poly_square (w,g,d_g);   d_w = ff_poly_degree(w,2*d_g);
        }
        ff_poly_mod_reduce (g,&d_g,w,d_w,mod);
        t--;
    }
    *pd_g = d_g;
}

void ff_poly_xn_mod_mpz (ff_t g[], int *pd_g, mpz_t e, ff_t f[], int d_f)
{
    ff_poly_modulus_t mod;
    register int i;
    
    if ( mpz_cmp_ui (e, d_f) < 0 )  { *pd_g = (int)mpz_get_ui(e);  _ff_set_one(g[*pd_g]);  for ( i = *pd_g-1 ; i >= 0 ; i-- ) _ff_set_zero(g[i]);  return; }

    ff_poly_mod_setup (mod,f, d_f);
    ff_poly_xn_modulus_mpz (g, pd_g, e,mod);
    ff_poly_mod_clear(mod);
}


// computes (x+a)^n mod f for depressed monic f of deg > 1
void ff_poly_xpan_mod (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_t f[], int d_f)
{
    register int i;

    assert ( d_f > 1 && _ff_one(f[d_f]) && _ff_zero(f[d_f-1]) );
    if ( ! e ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }

    if ( d_f <= FF_POLY_XNMOD_SMALL_DEGREE ) {
        ff_t nf[FF_POLY_XNMOD_SMALL_DEGREE+1];

        // for small mod poly operations all assume the modulus is of the form x^n  - c_{n-2}x^{n-2} - ... - c_1*x - c_0, so negate the non-leading coefficients (leading one is implicit)
        for ( i = 0 ; i < d_f ; i++ ) _ff_neg(nf[i],f[i]);
        switch (d_f) {
        case 2: err_printf ("Degree 2 not supported in ff_poly_xpan_mod_ui, use ff_exp_ui\n"); abort();
        case 3: ff_poly_xpan_mod_d3 (g,a,e,nf); break;
        case 4: ff_poly_xpan_mod_d4 (g,a,e,nf); break;
        case 5: ff_poly_xpan_mod_d5 (g,a,e,nf); break;
        case 6: ff_poly_xpan_mod_d6 (g,a,e,nf); break;
        case 7: ff_poly_xpan_mod_d7 (g,a,e,nf); break;
        case 8: ff_poly_xpan_mod_d8 (g,a,e,nf); break;
        case 10: ff_poly_xpan_mod_d10 (g,a,e,nf); break;
        case 11: ff_poly_xpan_mod_d11 (g,a,e,nf); break;
        case 12: ff_poly_xpan_mod_d12 (g,a,e,nf); break;
        case 13: ff_poly_xpan_mod_d13 (g,a,e,nf); break;
        case 15: ff_poly_xpan_mod_d15 (g,a,e,nf); break;
        case 17: ff_poly_xpan_mod_d17 (g,a,e,nf); break;
        case 19: ff_poly_xpan_mod_d19 (g,a,e,nf); break;
        case 23: ff_poly_xpan_mod_d23 (g,a,e,nf); break;
        case 29: ff_poly_xpan_mod_d29 (g,a,e,nf); break;
        case 31: ff_poly_xpan_mod_d31 (g,a,e,nf); break;
        default: ff_poly_xpan_mod_small (g,a,e,nf,d_f); break;
        }
        *pd_g = ff_poly_degree(g,d_f-1);
        return;
    } else {
        ff_poly_modulus_t mod;

        ff_poly_mod_setup (mod, f, d_f);
        ff_poly_xpan_modulus (g, pd_g, a,e,mod);
        ff_poly_mod_clear(mod);
    }
}

void ff_poly_xpan_modulus  (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_poly_modulus_t mod)
{
    register ff_t *w,t0;
    register int i, t, d_w;
    int d_g;

//printf ("Computing (x+%ld)^%d mod f (p=%ld) ", _ff_get_ui(a), e, _ff_p);  ff_poly_print(f,d_f);
    _ff_set_one(g[1]);  _ff_set(g[0],a); d_g = 1;
    w = mod->w1;
    t = ui_len(e)-2;
    while ( t >= 0 ) {
        if ( e&(1UL<<t) ) {
            ff_poly_square (w+1,g,d_g);   
            _ff_mult(w[0],a,w[1]);
            for ( i = 1 ; i < 2*d_g+1 ; i++ ) { _ff_mult(t0,a,w[i+1]); _ff_addto(w[i],t0); }
            d_w = ff_poly_degree(w,2*d_g+1);
        } else {
            ff_poly_square (w,g,d_g);   d_w = ff_poly_degree(w,2*d_g);
        }
        ff_poly_mod_reduce (g,&d_g,w,d_w,mod);
        t--;
    }
    *pd_g = d_g;
}

// Algorithm IPT from Shoup "A Computational Introduction to Number Theory and Algebra" p. 463
// returns true if f is irreducible, false otherwise.  If pnroots is non-null then *pnroots is set to the number of Fp-roots (whether f is irreducible or not).
int ff_poly_irreducible (ff_t f[], int d_f, int *pnroots)
{
    ff_poly_modulus_t mod;
    ff_t X[2], *h;
    int d_h, d_w;
    int i;

    assert ( d_f > 0 && _ff_one(f[d_f]) );
    if ( d_f == 1 ) {  if ( pnroots ) *pnroots = 1;  return 1;  }
    
    _ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);        // create the poly -x
    h = ff_poly_stack_alloc (d_f-1);
    if ( ! _ff_one(f[d_f]) ) f = ff_poly_monic (ff_poly_stack_alloc(d_f), 0, f, d_f);
    ff_poly_mod_setup (mod, f, d_f);
    ff_poly_xn_modulus (h, &d_h, _ff_p, mod);  i = 1;           // h = x^p mod f
    for ( ;; ) {
        ff_poly_add (mod->w, &d_w, h, d_h, X, 1);           // w = x^(p^i)-x
        ff_poly_gcd (mod->w, &d_w, mod->w, d_w, f, d_f);        // w = gcd(t,f)
        if ( d_w || ++i > d_f/2 ) break;
        ff_poly_pow_modulus (h, &d_h, h, d_h, _ff_p, mod);      // h = x^(p^i) mod f
    }
    ff_poly_stack_pop (h);  ff_poly_mod_clear (mod);
    if ( i > d_f/2 ) return 1;
    if ( pnroots ) { *pnroots = ( i==1 ? d_w : 0 ); }
    return 0;
}

// computes the factorization pattern of a (not necessarily monic) polynomial f, sets counts[d] to the number of irreducible factors of degree d, returns total count
// if root is non-null and poly has at least one linear factor, r will be set to a root
int ff_poly_factorization_pattern_and_root (int counts[], ff_t f[], int d_f, ff_t *root)
{
    ff_poly_modulus_t mod;
    ff_t *g, *h, X[2];
    int d_g, d_h, d_w;
    int i, k, n;

    assert ( d_f > 0 );
    if ( d_f == 1 ) { if ( counts ) counts[1] = 1; if ( root ) ff_poly_find_root (root, f, d_f); return 1; }
    
    if ( counts ) for ( n =0, i = 0 ; i <= d_f ; i++ ) counts[i] = 0;
    _ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);    // create the poly -x, will be reused
    g = ff_poly_stack_alloc (d_f);  h = ff_poly_stack_alloc (d_f-1);
    ff_poly_monic (g, &d_g, f, d_f);                                // g is a working copy of f made monic
    ff_poly_mod_setup (mod, g, d_g);
    i = 1;
    ff_poly_xn_modulus (h, &d_h, _ff_p, mod);                       // h = x^p mod f
    for (;;) {
        ff_poly_add (mod->w, &d_w, h, d_h, X, 1);                   // w = x^{p^i}-x mod g
        for(;;) {
            ff_poly_gcd_reduce (mod->w, &d_w, mod->w, d_w, g, d_g);
            if ( ! d_w ) break;
            k = d_w/i;  assert (k*i==d_w);
            ff_poly_monic(mod->w,0,mod->w,d_w);
            if ( i == 1 && d_w && root && ! ff_poly_find_root (root, mod->w, d_w) ) { err_printf ("Error, couldn't find root of a totally split poly in ff_poly_factorization_pattern_with_root\n"); abort(); }
            ff_poly_div_exact (g, &d_g, g, d_g, mod->w, d_w);
            if ( counts ) counts[i] += k;
            n += k;
        }
        if ( ++i > d_g/2 ) break;
        // we could update the modulus here
        ff_poly_pow_modulus (h, &d_h, h, d_h, _ff_p, mod);         // h = x^(p^i) mod f
    }
    ff_poly_stack_pop(g); ff_poly_mod_clear (mod);
    if ( d_g ) { if ( counts ) counts[d_g]++; n++; }
    return n;
}


// completely factors a square-free monic polynomial f that is the product of d/k irreducible polynomials of degree k, using Cantor-Zassenhaus (Algs 14.8 and 14.10 in GG)
// r must have space for d=(d/k)*k coefficients:  a concatenated list of d/k implicitly monic polys of degree k will be stored in w, each using just d/k coeffs (monic leading coeff is omitted).
// r and f cannot overlap
void ff_poly_equal_degree_factorization (ff_t r[], ff_t f[], int d, int k)
{
    ff_poly_modulus_t mod;
    ff_t *g;
    mpz_t n;
    int d_g;

    assert ( k > 0 && d > 0 && !(d%k) && _ff_one(f[d]) );
    if ( d == k ) { ff_poly_copy (r, 0, f, d-1); return; }
//printf ("EFD: "); ff_poly_print (f,d);
    g = ff_poly_stack_alloc (d);  
    ff_poly_mod_setup (mod, f, d);
    mpz_init (n); mpz_ui_pow_ui (n, _ff_p, k); mpz_sub_ui(n, n, 1); mpz_div_2exp (n, n, 1);     // n = (p^k-1)/2
    do {        
        ff_poly_randomize (g, d-1);  d_g = ff_poly_degree (g, d-1);                             // pick a random poly mod f
        if ( d_g <= 0 ) continue;
        ff_poly_pow_modulus_mpz (g, &d_g, g, d_g, n, mod);                                      // fg= g^n mod
        if ( d_g < 1 ) continue;
        _ff_dec (g[0]);                                                                         // fg= g^n-1 mod 
        ff_poly_gcd_reduce (g, &d_g, g, d_g, f, d);
    } while ( d_g <= 0 );
    mpz_clear (n);
    ff_poly_mod_clear (mod);
    
    assert ( !(d_g%k) );
    
    ff_poly_monic (g, 0, g, d_g);
    ff_poly_equal_degree_factorization (r, g, d_g, k);
    ff_poly_div_exact (g, &d_g, f, d, g, d_g);
    ff_poly_equal_degree_factorization (r+d-d_g, g, d_g, k);
    ff_poly_stack_pop (g);
}


// completely factors a monic polynomial f, using Cantor-Zassenhaus (Algs 14.13 in GG)
// r must have space for d_f coefficients:  a concatenated list of implicitly monic polys f_1,...,f_t will be written to r, with n[i] set to deg f_i.
// r and f cannot overlap.  the number of factors t is returned
int ff_poly_factors (ff_t r[], int n[], ff_t f[], int d)
{
    ff_poly_modulus_t mod;
    ff_t *g, *h, *s, X[2];
    int d_g, d_h, d_s;
    int i, j, k, o, t;

    assert ( d > 0 && _ff_one(f[d]) && r != f );
//printf ("Factoring monic poly "); ff_poly_print (f, d);

    // first remove any factors of x
    for ( t = o = 0; _ff_zero(f[0]) ; o++ ) { f++;  d--;  _ff_set_zero(r[o]); n[t++] = 1; }
    if ( !d ) return t;
    if ( d == 1 ) { _ff_set(r[o],f[0]); n[t++] = 1; return t; }

    _ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);                // create the poly -x
    g = ff_poly_stack_alloc(d);
    h = ff_poly_stack_alloc(d);
    s = ff_poly_stack_alloc(d);
    ff_poly_copy (g, &d_g, f, d);                                   // g is a working copy of f that will be whittled down
    ff_poly_mod_setup (mod, f, d);
    ff_poly_xn_modulus (h, &d_h, _ff_p, mod);                       // h = x^p mod f
    i = 1;
    for (;;) {
        ff_poly_add (s, &d_s, h, d_h, X, 1);                            // w = x^{p^k}-x mod f
        for(;;) {
            ff_poly_gcd_reduce (s, &d_s, s, d_s, g, d_g);
            if ( ! d_s ) break;
            k = d_s/i;  assert (k*i == d_s);
            ff_poly_monic(s,0,s,d_s);
            ff_poly_div_exact (g, &d_g, g, d_g, s, d_s);
            ff_poly_equal_degree_factorization (r+o, s, d_s, i);
            for ( j = 0 ; j < k ; j++ ) n[t++] = i;
            o += d_s;
        }
        if ( ++i > d_g/2 ) break;
        ff_poly_pow_modulus (h, &d_h, h, d_h, _ff_p, mod);              // h = x^(p^k) mod f
    }
    if ( d_g > 0 ) { ff_poly_copy (r+o, 0, g, d_g-1);  n[t++] = d_g; }
    ff_poly_mod_clear (mod);
    ff_poly_stack_pop (g);
//printf ("Found %d factors of degrees: ", t);  for ( i = 0 ; i < t ; i++ ) printf ("%d ", n[i]); puts ("");
    return t;
}


int ff_poly_distinct_roots (ff_t r[], ff_t f[], int d_f)                        //  only returns distinct roots, f need not be monic
{
    ff_t w[d_f+1];
    int k;

    assert (d_f >= 0);
    if ( d_f > 2 ) f = ff_poly_monic (w, 0, f, d_f);

    switch (d_f) {
    case 0: return 0;
    case 1: return ff_poly_roots_d1(r,f);
    case 2: k=ff_poly_roots_d2(r,f,d_f);  if ( k==2 && _ff_equal(r[0],r[1]) ) k=1;  return k;
    case 3:
        k=ff_poly_roots_d3(r,f);
        if ( k==3 ) {
            if ( _ff_equal(r[0],r[1]) ) { if ( _ff_equal(r[1],r[2]) ) { k=1; } else { _ff_set (r[1],r[2]); k =2; } }
            else if ( _ff_equal(r[0],r[2]) ) k = 2;
            else if ( _ff_equal(r[1],r[2]) ) k = 2;
        }
        return k;
    case 4:
        k=ff_poly_roots_d4(r,f);
        if ( k==2 && _ff_equal(r[0],r[1]) ) k=1;
        if ( k==4 ) {
            if ( _ff_equal(r[0],r[1]) ) {
                if ( _ff_equal(r[1],r[2]) ) {
                    if ( _ff_equal(r[2],r[3]) ) {
                        k= 1;
                    } else {
                        _ff_set(r[1],r[3]);
                        k = 2;
                    }
                } else {
                    if ( _ff_equal(r[1],r[3]) ) {
                        _ff_set(r[1],r[2]);  k = 2;
                    } else {
                        _ff_set(r[1],r[3]); k = (_ff_equal(r[1],r[2])?2:3);
                    }
                }
            } else if ( _ff_equal(r[0],r[2]) ) {
                if ( _ff_equal(r[0],r[3]) || _ff_equal(r[1],r[3]) ) {
                    k = 2;
                } else {
                    _ff_set(r[2],r[3]); k= 3;
                }
            } else if ( _ff_equal(r[0],r[3]) ) {
                k = ( _ff_equal(r[1],r[2]) ? 2 : 3 );
            } else if ( _ff_equal(r[1],r[2]) ) {
                _ff_set(r[2],r[3]);
                k = ( _ff_equal(r[1],r[2]) ? 2 : 3 );
            } else if ( _ff_equal(r[1],r[3]) ) {
                k = 3;
            } else {
                k = ( _ff_equal(r[2],r[3]) ? 3 : 4 );
            }
        }
        return k;
    }
    k = _ff_poly_distinct_roots(r,f,d_f,0);
    return k;
}

// this function is not particularly optimized - an area for future work.  requires f monic
// split flag indicates that f is a product of linear polynomials (after the initial call, all recursive calls will have this flag set)
// FF_POLY_ONE_ROOT is used to indicate that one *randomly chosen* root is desired -- in this case r will only have space for one entry
// returns the number of distinct roots (even when FF_POLY_ONE_ROOT is set).
int _ff_poly_distinct_roots  (ff_t *r, ff_t f[], int d_f, int flags)
{
    ff_t *g, *h, X[2], c, c2;
    int d_g, d_h;
    int i, n;

    assert ( d_f >= 0 );
    if ( ! d_f ) return 0;

    if ( d_f < 5 ) {
        if ( (flags&FF_POLY_ONE_ROOT) ) {                                                   // AVS bug fix 03/06/2013, don't go past first entry of r when FF_POLY_ONE_ROOT is set
            ff_t rr[4];
            n = ff_poly_distinct_roots(rr,f,d_f);
            if ( n ) r[0] = rr[ff_randomm_ui(n)];
            return n;
        }
        return ff_poly_distinct_roots(r,f,d_f);
    }
    g = ff_poly_stack_alloc (d_f);  h = ff_poly_stack_alloc (d_f);
    if ( ! _ff_one(f[d_f]) ) f = ff_poly_monic (ff_poly_stack_alloc(d_f), 0, f, d_f);
    if ( (d_f < _ff_p) ) {
        ff_poly_depress_monic(&c,g,f,d_f);  d_g=d_f;                                        // depress f  to g(x) = f(x-c) to speed up mod f operations, if possible
    } else {
        ff_poly_copy (g, &d_g, f, d_f);  _ff_set_zero (c);
    }
    if ( ! (flags&FF_POLY_SPLIT) ) {
        ff_poly_xn_mod (h, &d_h, _ff_p, g, d_g);                                            // h = x^p mod f(x-c)
        _ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);
        ff_poly_addto (h, &d_h, X, 1);                                                      // h = x^p-x mod f(x-c)  (let addto handle the addition of -x to get the degree right)
        ff_poly_gcd_reduce (g, &d_g, h, d_h, g, d_g);                                       // g=gcd(g,h)
        ff_poly_monic (g, 0, g, d_g);
        if ( d_g < 5 ) { n = _ff_poly_distinct_roots(r,g,d_g,flags|FF_POLY_SPLIT);  goto translate_roots; }
        if ( (d_g%_ff_p) ) {
            ff_poly_depress_monic_inplace (&c2,g,d_g);
            _ff_addto(c,c2);
        }
    }
    n = d_g;
    do {
        ff_random(X);                                                                       // pick a random a
        ff_poly_xpan_mod (h, &d_h, X[0], _ff_p>>1, g, d_g);                                 // h = (x+a)^{(p-1)/2} mod g
        if ( d_h < 1 ) continue;
        _ff_dec(h[0]);                                                                      // h = (x+a)^{(p-1)/2}-1 mod g
        ff_poly_gcd_reduce (h, &d_h, h, d_h, g, d_g);                                       // h=gcd(g,h)
    } while ( d_h <= 0 );
    ff_poly_monic(h,0,h,d_h);
    if ( (flags&FF_POLY_ONE_ROOT) ) {
        if ( d_h > d_g-d_h ) ff_poly_div_exact (h, &d_h, g, d_g, h, d_h);
        _ff_poly_distinct_roots(r, h, d_h, flags|FF_POLY_SPLIT);
    } else {
        _ff_poly_distinct_roots(r, h, d_h, flags|FF_POLY_SPLIT);
        ff_poly_div_exact (h, &d_h, g, d_g, h, d_h);
        _ff_poly_distinct_roots(r+d_g-d_h, h, d_h, flags|FF_POLY_SPLIT);
    }
translate_roots:    // if we get here, n is not zero (in fact n is at least 5)
    if ( ! _ff_zero(c) ) {
        if ( n ) _ff_subfrom(r[0],c);
        if ( ! (flags&FF_POLY_ONE_ROOT) ) for ( i = 1 ; i < n ; i++ ) _ff_subfrom(r[i],c);
    }
    ff_poly_stack_pop (g);
    return n;
}

// returns all Fp-roots of f, with multiplicity.  f need not be monic.  Repeated roots will be adjacent to eachother in the returned list.
int ff_poly_roots (ff_t r[], ff_t f[], int d_f)
{
    ff_t *g, t0;
    int i, k, n, d_g;

    assert ( d_f >= 0 );
    if ( d_f <= 2 ) return ff_poly_roots_d2 (r, f, d_f);
    if ( d_f <= 4 ) {
        ff_t w[d_f+1];
        if ( ! _ff_one(f[d_f]) ) f = ff_poly_monic (w, 0, f, d_f);
        return d_f == 3 ? ff_poly_roots_d3(r,f) : ff_poly_roots_d4(r,f);
    }
    k = _ff_poly_distinct_roots(r,f,d_f,0);
    if ( !k || k == d_f ) return k;
    g = ff_poly_stack_alloc (d_f);  ff_poly_copy (g, &d_g, f, d_f);
    for ( i = 0, n = k ; i < k ; i++ ) {
        ff_poly_remove_root (g, g, d_g, r+i);  d_g--;
        for (;;) {
            ff_poly_eval (&t0, g, d_g, r+i);
            if ( ! _ff_zero(t0) ) break;
            ff_poly_remove_root (g, g, d_g, r+i);  d_g--;
            _ff_set (r[n],r[i]); n++;
        }
    }
    ff_organize (r, n); // return roots in canonical order
    ff_poly_stack_pop (g);
    return n;
}

// returns the degree of gcd(f,x^p-x) which will be the number of distinct roots.  does not require f to be monic
int ff_poly_count_distinct_roots (ff_t f[], int d_f)
{
    ff_t *g, *h;
    int d_g, d_h;

    if ( d_f <= 1 ) return d_f;
    if ( d_f < 5 ) {
        ff_t r[4], w[5];
        if ( ! _ff_one(f[d_f]) ) f = ff_poly_monic (w, 0, f, d_f);
        return ff_poly_distinct_roots (r, f, d_f);
    }

    g = ff_poly_stack_alloc (d_f);  h = ff_poly_stack_alloc (d_f-1);
    ff_poly_monic (g, &d_g, f, d_f);
    if ( (d_f%_ff_p) ) ff_poly_depress_monic_inplace (0, g, d_g);
    ff_poly_xn_mod (h, &d_h, _ff_p, g, d_g);                    // h = x^p mod g
    ff_poly_sub_x_from (h, &d_h);                               // h = x^p-x mod g
    ff_poly_gcd_reduce (g, &d_g, g, d_g, h, d_h);
    ff_poly_stack_pop (g);
    return d_g;
}


// returns the total number of Fp-roots of f, with multiplicity.
int ff_poly_count_roots (ff_t f[], int d_f)
{
    ff_t *g, *h;
    int d_g, d_h, n;

    if ( d_f <= 1 ) return d_f;
    if ( d_f < 5 ) {
        ff_t r[4], w[5];
        if ( ! _ff_one(f[d_f]) ) f = ff_poly_monic (w, 0, f, d_f);
        return ff_poly_roots (r, f, d_f);
    }

    g = ff_poly_stack_alloc (d_f);  h = ff_poly_stack_alloc (d_f);
    ff_poly_monic (g, &d_g, f, d_f);
    if ( (d_f%_ff_p) ) ff_poly_depress_monic_inplace (0, g, d_g);
    ff_poly_xn_mod (h, &d_h, _ff_p, g, d_g);                    // h = x^p mod g
    ff_poly_sub_x_from (h, &d_h);                               // h = x^p-x mod g
    ff_poly_gcd_reduce (h, &d_h, g, d_g, h, d_h);               // h = gcd(g,h)
    for ( n = d_h ; d_h ; n += d_h ) {
        ff_poly_div_exact (g, &d_g, g, d_g, h, d_h);
        ff_poly_gcd_reduce (h, &d_h, g, d_g, h, d_h);
    }
    ff_poly_stack_pop (g);
    return n;
}


int ff_poly_square_free_part (ff_t g[], int *pd_g, ff_t f[], int d_f)
{
    ff_t *h, *s, *t, *w, k;
    int d_g, d_h, d_r, d_s, d_t, d_w;
    
    assert ( d_f > 0 );
    h = ff_poly_stack_alloc(d_f); s = ff_poly_stack_alloc(d_f); t = ff_poly_stack_alloc (d_f); w = ff_poly_stack_alloc (d_f);
    ff_poly_monic (g, &d_g, f, d_f);
    ff_poly_derivative (h, &d_h, g, d_g);
    ff_poly_gcd_reduce (s, &d_s, g, d_g, h, d_h);
    ff_poly_monic (s, 0, s, d_s);
    ff_poly_div (t, &d_t, 0, &d_r, g, d_g, s, d_s);
    _ff_set_one(g[0]); d_g = 0;
    _ff_set_one(k);
    while ( d_t > 0 ) {
        if ( ! _ff_zero(k) ) {
            ff_poly_gcd_reduce (h, &d_h, s, d_s, t, d_t);
            ff_poly_monic (h, 0, h, d_h);
            if ( d_h < d_t ) {
                ff_poly_div (w, &d_w, 0, &d_r, t, d_t, h, d_h);  assert (d_r < 0);
                ff_poly_mult (g, &d_g, g, d_g, w, d_w);
            }
            ff_poly_copy (t, &d_t, h, d_h);
        }
        ff_poly_div (s, &d_s, 0, &d_r, s, d_s, h, d_h); assert (d_r < 0);
         _ff_inc(k);
    }
    *pd_g = d_g;
    if ( d_s > 0 ) {
        assert ( d_s % _ff_p == 0 );  d_s /= _ff_p;
        for ( int i = 0 ; i <= d_s ; i++ ) _ff_set(s[i],s[_ff_p*i]);
        ff_poly_square_free_part (t, &d_t, s, d_s);
        ff_poly_mult (g, pd_g, g, d_g, t, d_t);
    }
    ff_poly_stack_pop(h);
    return *pd_g;
}

// computes g = GCD(x^(p^n)-1,f) = product of all distinct factors of f of degree <= n
// f need not be monic, it will be made monic if necessary, g and f may overlap
// intended for small d and large p, does not use modulus
int ff_poly_gcd_xpn (ff_t g[], int *pd_g, int n, ff_t f[], int d_f, int comp)
{
    ff_t *h, *w, t;
    int i, d_h;

    assert (d_f > 1 && n >= 1);
    
    ff_poly_monic (g, 0, f, d_f);
    if ( ! _ff_zero(g[d_f-1]) && (d_f%_ff_p) ) ff_poly_depress_monic_inplace (&t, g, d_f); else _ff_set_zero(t);
    
    h = ff_poly_stack_alloc (d_f-1);
    if ( comp && n > 1 ) {
        ff_poly_xn_mod (h, &d_h, _ff_p, g, d_f);                    // h = x^p mod g
        if ( n > 1 ) {
            ff_poly_zpad (h, d_h, d_f-1);                           // zpad mod g
            w = ff_poly_stack_alloc (d_f-1);
            ff_poly_copy (w, 0, h, d_f-1);                          // save copy for repeated composition
            for ( i = 1 ; i < n ; i++ )                             // note that compose_mod_small zero pads result but computes d_h correctly
                ff_poly_compose_mod_small (h, &d_h, h, d_f-1, w, d_f-1, g, d_f);
        }
    } else {
        if ( n * ui_len(_ff_p) <= FF_BITS ) {
            unsigned long q; for ( q = _ff_p, i = 1 ; i < n ; i++ ) q *= _ff_p;
            ff_poly_xn_mod (h, &d_h, q, g, d_f);
        } else {
            mpz_t q;  mpz_init (q);  mpz_ui_pow_ui (q, _ff_p, n);
            ff_poly_xn_mod_mpz (h, &d_h, q, g, d_f);
            mpz_clear (q);
        }
    }
    ff_poly_sub_x_from (h, &d_h);                               // h = x^(p^n)-x mod g
    ff_poly_gcd_reduce (g, pd_g, g, d_f, h, d_h);               // g = GCD(g,x^(p^n)-x)
    ff_poly_stack_pop (h);
    if ( ! _ff_zero(t) ) ff_poly_translate (g, g, *pd_g, t);
    return *pd_g;
}

// sets cnts[i] to the number of distinct facotrs of f(x) of degree i for i in [1,maxd], returns total
// if d_f is significantly larger thatn log(p), using composition is sub-optimal
int ff_poly_count_distinct_factors (int cnts[], ff_t f[], int d_f, int maxd)
{
    ff_t *g, *h, *w, *t;
    int i, j, k, n, d_g, d_h, d_t, s;

    if ( maxd < 1 ) return 0;
    if ( d_f <= 1 ) return d_f;
    memset (cnts, 0, (maxd+1)*sizeof(*cnts));
    if ( d_f == 1 ) { cnts[1] = 1; return 1; }
    if ( d_f == 2 ) {  ff_t w[2]; cnts[1] = ff_poly_distinct_roots(w, f, d_f); cnts[2] = ( cnts[1] ? 0 : 1 ); return cnts[1]+cnts[2]; }
    if ( maxd == 1 ) { cnts[1] = ff_poly_count_distinct_roots(f, d_f); return cnts[1]; }
    if ( d_f <= 4 ) {
        ff_t g[5], w[4];

        ff_poly_monic (g, &d_g, f, d_f);
        if ( ! _ff_zero(g[d_g-1]) && (d_g%_ff_p) ) ff_poly_depress_monic_inplace (0, g, d_g);
        
        cnts[1] = ff_poly_distinct_roots (w, g, d_g);
        s = ff_poly_disc_legendre_small (g, d_g);
        if ( d_f == 3 ) {
            if ( !s ) return cnts[1];
            if ( ! cnts[1] ) { if ( maxd > 2 ) cnts[3] = 1; return cnts[3]; }
            if ( cnts[1] == 1 ) { cnts[2] = 1; return 2; }
            return cnts[1];
        }
        if ( !s ) {
            if ( ! cnts[1] ) { cnts[2] = 1; return 1; }
            if ( cnts[1] > 1 ) { cnts[2] = 0; return cnts[1]; }
            // need to determine whether we have a double or quadruple root
            ff_poly_remove_root_d4 (g, g, w); ff_poly_remove_root_d3 (g, g, w);
            ff_poly_eval(w,g,2,w);
            cnts[2] = ( _ff_zero(w[0]) ? 0 : 1 );
            return cnts[1]+cnts[2];
        }
        if ( cnts[1] ) { cnts[4-cnts[1]] = 1; cnts[0] = 0; return cnts[1]+cnts[2]+cnts[3]; }
        if ( s > 0 ) cnts[2] = 2; else cnts[4] = 1;
        return cnts[2] + cnts[4];
        // we could also handle degree 5 using the discriminant -- this could be well worth it when p is large
    }

    g = ff_poly_stack_alloc (d_f);  h = ff_poly_stack_alloc (d_f); w = ff_poly_stack_alloc(d_f); t = ff_poly_stack_alloc(d_f);
    
    ff_poly_monic (g, &d_g, f, d_f);
    if ( ! _ff_zero(g[d_g-1]) && (d_g%_ff_p) ) ff_poly_depress_monic_inplace (0, g, d_g);

    n = ( maxd > d_g/2 ? d_g/2 : maxd );
    
    ff_poly_xn_mod (h, &d_h, _ff_p, g, d_g);                    // h = x^p mod g
    ff_poly_zpad (h, d_h, d_g-1);                               // zpad mod g
    ff_poly_copy (w, 0, h, d_g-1);                              // save copy for repeated composition
    ff_poly_copy (t, &d_t, h, d_h);
    ff_poly_sub_x_from (t, &d_t);                           // t = x^(p^i)-x mod g
    ff_poly_gcd_reduce (t, &d_t, g, d_g, t, d_t);           // t = GCD(g,x^(p^i)-x)
    cnts[1] = d_t;
    for ( i = 2 ; i <= n ; i++ ) {
        ff_poly_compose_mod_small (h, &d_h, h, d_h, w, d_g-1, g, d_g);
        ff_poly_copy (t, &d_t, h, d_h);
        ff_poly_sub_x_from (t, &d_t);                           // t = x^(p^i)-x mod g
        ff_poly_gcd_reduce (t, &d_t, g, d_g, t, d_t);           // t = GCD(g,x^(p^i)-x)
        for ( k = 0, j = 1 ; j < i ; j++ ) if ( ! (i%j) ) k += j*cnts[j];
        assert ( (d_t - k) % i == 0 );
        cnts[i] = (d_t - k) / i;
    }

    // if maxd > d_f/2 use what we know to determine if there is a factor of degree greater than d_f/2 and less than or equal to maxd (note that f might not be square free)
    if ( n < maxd ) {
        for ( k = d_f, i = 1 ; i <= maxd ; i++ ) k -= i*cnts[i];
        if ( k >  n ) {
            if ( ff_poly_squarefree (g, d_g) ) {
                if ( k <= maxd ) cnts[k] = 1;
            } else {
                ff_poly_square_free_part (g, &d_g, g, d_g);
                k -= (d_f-d_g);
                if (  k > n && k <= maxd ) cnts[k] = 1;
            }
        }
    }
    for ( k = 0, i = 1 ; i <= maxd ; i++ ) k += cnts[i];
    ff_poly_stack_pop(g);
    return k;
}

// very fast code to test whether gcd(f,g)=1 or not
int ff_poly_trivial_gcd (ff_t f[], int d_f, ff_t g[], int d_g)
{
    ff_t *s, *t;
    register ff_t a, b;
    register ff_t *sp, *se, *tp, *te, *xp;
    register int i;
    
    s = ff_poly_stack_alloc (d_f); t = ff_poly_stack_alloc (d_g);

    se = s+d_f;
    for ( sp = se, tp=f+d_f ; sp >= s ; ) *sp-- = *tp--;
    te = t+d_g;
    for ( tp = te, sp=g+d_g ; tp >= t ; ) *tp-- = *sp--;

    while ( se > s && te > t ) {
        i = (se-s) - (te-t);
        if ( i >= 0 ) {
            a = *te;  b = *se;  xp = s+i;
            ff_negate(b);
            for ( sp = s ; sp < xp ; ) { _ff_mult (*sp,*sp,a); sp++; }
            for ( tp = t ; tp <= te ; tp++ ) { _ff_sum_2_mults (*sp, a, b, *tp, *sp); sp++; }
            for ( se-- ; se >= s && !*se ; se-- );
        } else {
            a = *se;  b = *te;  xp = t-i;
            ff_negate(b);
            for ( tp = t ; tp < xp ; ) { _ff_mult (*tp,*tp,a); tp++; }
            for ( sp = s ; sp <= se ; sp++ ) { _ff_sum_2_mults (*tp, a, b, *sp, *tp); tp++;  }
            for ( te-- ; te >= t && !*te ; te-- );
        }
    }
    ff_poly_stack_pop(s);
    if ( se < s || te < t ) return 0;
    return 1;
}

// tests whether discriminant is zero by computing gcd(f,f')
int ff_poly_squarefree (ff_t f[], int deg)
{
    ff_t *df;
    int b, d_df;

    assert ( deg > 0 );
    
    // check for easy special cases first
    if ( _ff_one(f[deg]) && deg  <= 3 ) { ff_t D; ff_poly_disc_small (&D, f, deg); return ( _ff_zero(D) ? 0 : 1 ); }
    
    df = ff_poly_derivative (ff_poly_stack_alloc (deg-1), &d_df, f, deg);
    b = ff_poly_trivial_gcd (f, deg, df, d_df);
    ff_poly_stack_pop(df);
    return b;
}


void ff_poly_disc (ff_t D[1], ff_t f[], int d)
{
    ff_t *df, *M, t, w;
    register int i, k;

    assert (d > 0 && ! _ff_zero(f[d]) );
    switch (d) {
    case 1: _ff_set_one(D[0]); return;
    case 2: if ( _ff_one(f[d]) ) { ff_poly_disc_2 (D, f); return; } else break;
    case 3: if ( _ff_one(f[d]) ) { ff_poly_disc_3 (D, f); return; } else break;
    case 4: if ( _ff_one(f[d]) ) { ff_poly_disc_4 (D, f); return; } else break;
    }

    // we compute the determinant of the d x d matrix with rows f', x*f' mod f, x^2*f' mod f, ...,  x^(d-1)*f' mod f, where coeffs are store from high to low degree
    // this amounts to row-reducing the Sylvester matrix and reversing the order of the rows (which eliminates the need to fuss with the sign at the end)
    M = ff_poly_stack_alloc (d*d-1);
    df = ff_poly_derivative (ff_poly_stack_alloc(d-1),0,f,d);
    _ff_set_ui(t,d);
    
    for ( i = 0 ; i < d ; i++ ) _ff_set(M[i],df[d-i-1]);                                  // first row is f'
    for ( i = 0 ; i < d-1 ; i++ ) { _ff_mult(w,t,f[d-i-1]); _ff_sub(M[d+i],M[i+1],w); }   // second row is x*f'-d*f
    _ff_mult(w,t,f[0]); _ff_neg (M[2*d-1],w);
    for ( k = 2 ; k < d ; k++ ) {
        for ( i = 0 ; i < d-1 ; i++ ) {
            _ff_diff_2_mults (M[k*d+i],f[d],f[d-i-1],M[(k-1)*d],M[(k-1)*d+i+1]);          // kth row is lc(f)*x*(row k-1) - lc(row k-1)*f
        }
        _ff_mult (w,M[(k-1)*d],f[0]);  _ff_neg (M[k*d+d-1],w);                            // last entry of kth row
    }
    ff_matrix_destructive_determinant(M,d);                                               // compute determinant using ~ 2*d^3/3 mults
    ff_exp_ui (&t,f+d,(d-2)*(d-3)/2);
    ff_invert(t,t);
    _ff_mult(D[0],t,M[0]);
    ff_poly_stack_pop(M);
}

void ff_poly_resultant (ff_t D[1], ff_t f[], int df, ff_t g[], int dg)
{
    ff_t *M;
    register int i,j;
    
    if ( df < 0 || dg < 0 ) { _ff_set_zero(D[0]); return; }
    if ( !df ) { ff_exp_ui(D,f,dg); return; }
    if ( !dg ) { ff_exp_ui(D,g,df); return; }
    M = ff_poly_stack_alloc((df+dg)*(df+dg)-1);
    memset(M,0,(df+dg)*(df+dg)*sizeof(M[0]));
    for ( i = 0 ; i < dg ; i++ ) for ( j = 0 ; j <= df ; j++ ) _ff_set(M[i*(df+dg)+i+df-j],f[j]);
    for ( ; i < df+dg ; i++ ) for ( j = 0 ; j <= dg ; j++ ) _ff_set(M[i*(df+dg)+i-dg+dg-j],g[j]);
    ff_matrix_destructive_determinant(M,df+dg);
    _ff_set(D[0],M[0]);
    ff_poly_stack_pop(M);
}


/*
    Given monic f(x) of degree d, computes linear translation g(x)=(x-f[d-1]/d) which will be a monic poly with d-1 coefficient zero.
    The parameter s is optional, if specified it will be set to f[d-1]/d so that g(x)=f(x-s) (and g(x+s)=f(x)).  Replaces f by g.
*/
void _ff_poly_depress_monic_inplace (ff_t s[1], ff_t f[], int d_f)
{
    ff_t z;
    register ff_t t,c;
    register int i, j;

    assert (d_f > 0);
    ff_invert_small_int(&z,d_f);
    _ff_mult(c,f[d_f-1],z);
    if ( s ) _ff_set(*s,c);                             // g(x) = f(x-s) so that g(x+s)=f(x)    consistent with ff_poly_depress_cubic
    ff_negate(c);                                   // c = -f[d-1]/d, we will replace f by f(x+c)
    for ( i = d_f ; i ; i-- ) {
        for ( j = i ; j < d_f ; j++ ) { _ff_mult(t,c,f[j]); _ff_addto(f[j-1],t); }
        _ff_addto(f[j-1],c);
    }
}


/*
    For monic f, substitute x <- (x - f_{d-1}/[d*f_d]) to make f_{d-1} term zero if necessary
    Returns 1 if poly modified, 0 if not
*/
int mpz_poly_depress_monic_inplace (mpz_t f[], int d)
{
    mpz_t c, s, t, w;
    int i, j;
    
    if ( d < 1 ) return 0;
    if ( ! mpz_sgn(f[d-1]) ) return 0;
    if ( mpz_cmp_ui (f[d],1) != 0 ) return 0;

    // dynamically allocate mpz's here, we don't expect to be called more than once per curve
    mpz_init (c); mpz_init (s); mpz_init (t); mpz_init (w);
    
    // set s to -d*f[d-1]
    mpz_mul_ui (s, f[d-1], d);  mpz_neg (s, s);

    // compute d^(2d) * f ( (x+sd) / d^2 ) = sum_{i=0}^d sum_{k=i}^d binom(k,i) * d^(2*d-2*j) * s^(j-i) * f[j] x^i
    for ( i = 0 ; i < d-1 ; i++ ) {
        mpz_set_ui (w, 1);
        mpz_set_ui (c, 0);
        for ( j = i ; j <= d ; j++ ) {
            mpz_ui_pow_ui (t, d, 2*(d-j));
            mpz_mul (t, t, f[j]);
            mpz_mul (t, t, w);
            mpz_mul_ui (t, t, ui_binomial (j, i));
            mpz_add (c, c, t);
            mpz_mul (w, w, s);
        }
        mpz_set (f[i], c);
    }
    mpz_set_ui (f[d-1], 0);
    
    mpz_clear (c); mpz_clear (s); mpz_clear (t); mpz_clear (w);

    return 1;
}


// uses recursive divide-and conquer algorithm
void ff_poly_from_roots (ff_t f[], ff_t r[], int d)
{
    int df, dg, dh;
    ff_t *g, *h;
    if ( d <= FF_POLY_FROMROOTS_SMALL_DEGREE ) { ff_poly_from_roots_small (f,r,d); return; }
    dg = (d+1)/2;  dh = d/2;
    g = ff_poly_stack_alloc(dg); h = ff_poly_stack_alloc(dh);
    ff_poly_from_roots(g,r,dg); ff_poly_from_roots(h,r+dg,dh);
    ff_poly_mult(f,&df,g,dg,h,dh);
    ff_poly_stack_pop(g);
}


/*
    Use hardwired fast schoolbook for multiplying small polys.
    Uses ff_poly_mult_big if we are above the crossover defined by FF_POLY_MULT_SMALL_DEGREE
    We prefer d_f <= d_g but will swap if needed
    aliasing is OK
*/
void ff_poly_mult (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
    ff_t *x, *y;
    register int i, d;
    if ( f == g && d_f == d_g ) { ff_poly_square (h, f, d_f); if ( pd_h ) *pd_h = ff_poly_degree(h,2*d_f); return; }
    if ( d_f < 0 || d_g < 0 ) { *pd_h = -1; return; }
    if ( d_f > d_g ) { _swap (d_g, d_f, d); _swap (f, g, x); }
    if ( d_g < 10 ) {
        ff_t t[10], u[20];
        if ( d_f < d_g ) {
            for ( i = 0 ; i <= d_f ; i++ ) _ff_set(t[i],f[i]);
            while ( i <= d_g ) { _ff_set_zero(t[i]); i++; }
            x = t;  y = u;
        } else {
            x = f;  y = h;
        }
        switch ( d_g ) {
        case 0: _ff_mult(h[0],f[0],g[0]); break;
        case 1: ff_poly_mult_2_2 (y, x, g); break;
        case 2: ff_poly_mult_3_3 (y, x, g); break;
        case 3: ff_poly_mult_4_4 (y, x, g); break;
        case 4: ff_poly_mult_5_5 (y, x, g); break;
        case 5: ff_poly_mult_6_6 (y, x, g); break;
        case 6: ff_poly_mult_7_7 (y, x, g); break;
        case 7: ff_poly_mult_8_8 (y, x, g); break;
        case 8: ff_poly_mult_9_9 (y, x, g); break;
        case 9: ff_poly_mult_10_10 (y, x, g); break;
        }
    } else if ( d_g <= FF_POLY_MULT_SMALL_DEGREE ) {
        ff_t t[FF_POLY_MULT_SMALL_DEGREE+1], u[2*FF_POLY_MULT_SMALL_DEGREE+1];
        // pad f out to the same size as g if needed
        if ( d_f < d_g ) {
            for ( i = 0 ; i <= d_f ; i++ ) _ff_set(t[i],f[i]);
            while ( i <= d_g ) { _ff_set_zero(t[i]); i++; }
            x = t;  y = u;
        } else {
            x = f;  y = h;
        }
        ff_poly_mult_small(y,x,g,d_g);
    } else {
        y = h;
        ff_poly_mult_big(h,f,d_f,g,d_g);
    }
    if ( pd_h ) *pd_h = ff_poly_degree (y,d_f+d_g);
    if ( y != h ) ff_poly_copy (h, 0, y, d_f+d_g);  // copy everything to make sure we zero pad to d_f+d_g
}


/*
    Compute the inverse of f mod x^n using Algorithm 9.3 in von zur Gathen & Gerhard.
    Incorporates optimization to efficiently treat n not a power of 2.
    g and f cannot overlap
*/
void  ff_poly_inverse_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n)
{
    ff_t *h;
    register ff_t t0;
    register int i, j, k, d, e, r;
    int d_h;
    
    if ( ! n ) { *d_g = -1; return; }
    if ( _ff_zero(f[0]) ) { err_printf ("Uninvertible polynomial in ff_poly_inverse_mod_xn!\n"); ff_poly_print(f,d_f); abort(); }
    if ( _ff_one(f[0]) ) _ff_set_one(g[0]); else _ff_invert(g[0],f[0]);
    if ( n==1 ) { *d_g = 0; return; }
    h = ff_poly_stack_alloc(2*n);
    d = 1; r = ui_len(n);  if ( (1<<(r-1)) == n ) r--;
    for (r--;;r--) {
        e = ui_ceil_ratio (n, (1<<r));
        ff_poly_square(h,g,d-1);
        j = ( 2*d-2 < e ? 2*d-2 : e-1 );
        k = ( d_f < e ? d_f : e-1 );
        ff_poly_mult(h,&d_h,h,j,f,k);
        for ( i = 0 ; i < d ; i++ ) {
            _ff_add(t0,g[i],g[i]);
            _ff_sub(g[i],t0,h[i]);
        }
        if ( e > n ) e = n;
        for ( ; i < e ; i++ ) _ff_neg(g[i],h[i]);
        if ( e == n ) break;
        d = e;
    }
    ff_poly_stack_pop(h);
    *d_g = ff_poly_degree(g,n-1);
}

// setup up data structure to efficiently reduce polys of degree up to max modulo g, see GG Alg 9.5.  max defaults to 2*d_g-1
void ff_poly_mod_setup_max (ff_poly_modulus_t mod, ff_t *g, int d_g, int max)
{
    register int i;
    int d, dw;
    
    if ( max < 2*d_g-1 ) max = 2*d_g-1;
    mod->d = d_g;  mod->e = max-d_g;
    mod->g = ff_poly_alloc (mod->d);
    mod->rgi = ff_poly_alloc(mod->e);
    mod->w = ff_poly_alloc(max);
    mod->w1 = ff_poly_alloc(max);
    ff_poly_monic (mod->g, 0, g, d_g);
    ff_poly_reverse (mod->w, &dw, mod->g, mod->d);
    ff_poly_inverse_mod_xn (mod->rgi, &d, mod->w, dw, mod->e+1);                    // compute the inverse of rev_d(g) mod x^{e+1}
    for ( i = d+1 ; i <= mod->e ; i++ ) _ff_set_zero(mod->rgi[i]);                  // zero-pad out to expected degree e
}

void ff_poly_mod_clear (ff_poly_modulus_t mod) 
    { ff_poly_free (mod->g); ff_poly_free (mod->rgi); ff_poly_free (mod->w); ff_poly_free (mod->w1);  mod->g = 0; }

void ff_poly_mod_reduce (ff_t *h, int *d_h, ff_t *f, int d_f, ff_poly_modulus_t mod)
{
    register ff_t t0;
    register int i, n;
    int d;
    
//printf ("Reducing "); ff_poly_print(f,d_f);
//printf ("Modulo "); ff_poly_print(mod->g,mod->d);
    
    // if f is already reduced, just copy it
    if ( d_f < mod->d ) { *d_h = ff_poly_degree(f,d_f); if ( h != f ) for ( i = *d_h ; i >= 0 ; i-- ) _ff_set(h[i],f[i]);  return; }

    n = mod->d+mod->e;
    assert ( d_f <= n );
    
    // compute rev_n(f), zero-padding f as required
//printf ("f = "); ff_poly_print(f,d_f);
    for ( i = 0 ; i <= n ; i++ ) if ( n-i <= d_f ) _ff_set(mod->w[i],f[n-i]); else _ff_set_zero(mod->w[i]);
//printf ("rev_%d(f) = ", n); ff_poly_print(mod->w,n);
//printf ("(rev_%d(g))^-1 mod x^{%d+1} = ", mod->d, mod->e); ff_poly_print(mod->rgi, mod->e);
    
    // compute q* = rev_n(f) * (rev_d g)^-1 (which we only need mod x^{e+1})
    ff_poly_mult (mod->w,&d,mod->rgi,mod->e,mod->w,mod->e);
//printf ("rev_%d(f) * (rev_%d(g))^-1 mod x^{%d+1} = ", n, mod->d, mod->e); ff_poly_print(mod->w,mod->e);

    // compute the quotient q = rev_e(q*)
    for ( d++ ; d <= mod->e ; d++ ) _ff_set_zero(mod->w[d]);
    for ( i = 0 ; i <= mod->e/2 ; i++ ) { _ff_set(t0, mod->w[i]); _ff_set(mod->w[i],mod->w[mod->e-i]); _ff_set(mod->w[mod->e-i],t0); }
//printf ("q = "); ff_poly_print(mod->w, mod->e);

    // compute the remainder as f - qg, which we know has degree < mod->d
    n = ( mod->e < mod->d ? mod->e : mod->d-1 );
    ff_poly_mult(mod->w,&d,mod->w,n,mod->g,mod->d-1); 
    while ( d < mod->d-1 ) { d++; _ff_set_zero(mod->w[d]); }        // zero pad out to expected degree for subtraction
    for ( i = 0 ; i < mod->d ; i++ ) ff_sub(h[i],f[i],mod->w[i]);
    *d_h = ff_poly_degree(h,mod->d-1);
//printf("f-qg = "); ff_poly_print(h,*d_h);
}

void ff_poly_mod_big (ff_t *h, int *d_h, ff_t *f, int d_f, ff_t *g, int d_g)
{
    ff_poly_modulus_t mod;
    
    ff_poly_mod_setup_max(mod,g,d_g,d_f);
    ff_poly_mod_reduce(h,d_h,f,d_f,mod);
    ff_poly_mod_clear(mod);
}

void ff_poly_compose_mod_big (ff_t g[], int *pd_g, ff_t h1[], int d_h1, ff_t h2[], int d_h2, ff_t f[], int d_f)
{
    ff_poly_modulus_t mod;
    ff_t *w;
    int i, d_w;

    assert ( d_f > 0 );
    if ( d_h1 < 1 ) { ff_poly_copy (g, pd_g, h1, d_h1); return; }
    ff_poly_mod_setup (mod, f, d_f);
    w = ff_poly_stack_alloc (2*d_f);
    ff_poly_scalar_mult (w, 0, h1[d_h1], h2, d_h2); d_w = d_h2;
    ff_poly_addto (w, &d_w, h1+d_h1-1, 0);
    for ( i = d_h1-2 ; i >= 0 ; i-- ) { 
        ff_poly_mult (w, &d_w, w, d_w, h2, d_h2);
        ff_poly_addto (w, &d_w, h1+i, 0);
        ff_poly_mod_reduce (w, &d_w, w, d_w, mod);
    }
    ff_poly_mod_clear(mod);
    ff_poly_copy (g, pd_g, w, d_w);
    ff_poly_stack_pop (w);
}

// g(x)=f(x+c), works in place
// Uses (d*(d-1)/2 + 1) mults (and reductions), d*(d-1)/2+d adds.  Bit complexity is O(d^2*M(log p))
void ff_poly_translate_small (ff_t g[],ff_t f[], int d, ff_t c)
{
    register ff_t t,w;
    register int j, k;
    
    if ( g != f ) ff_poly_copy (g,0,f,d);
    _ff_mult(w,g[d],c);
    for ( j = d-1 ; j >= 0 ; j-- ) { for ( k = j ; k < d-1 ; k++ ) { _ff_mult (t, g[k+1], c); _ff_addto (g[k], t); } _ff_addto (g[k],w); }
}

#define FF_POLY_TRANSLATE_BASE  4

// g(x)=f(x+c) using divide+conquer algorithm adapted from https://hal.inria.fr/ensl-00546102/document
// Bit complexity is O(M(d*log p)log d)
void ff_poly_translate_rec (ff_t g[], ff_t f[], int d, ff_t c)
{
    ff_t *h, *G[2], *w;
    int *split;
    register int i, ii, j, k;
    int d0, d1, e;

    if ( d <= FF_POLY_TRANSLATE_BASE ) { ff_poly_translate_small (g,f,d,c); return; }
    split = malloc (d*sizeof(*split));  // more space than we need but no harm
    for ( split[0] = d+1, k = 1 ; split[k-1] > FF_POLY_TRANSLATE_BASE ; k*=2 )
        for ( j = k-1 ; j >= 0 ; j-- ) { split[2*j+1] = (split[j]+1)/2; split[2*j] = split[j]/2; }

    h = ff_poly_stack_alloc (d);
    G[0] = ff_poly_stack_alloc ((d+1)/2); G[1] = ff_poly_stack_alloc((d+2)/2);
    w = ff_poly_stack_alloc (d);
    ff_poly_copy (h,0,f,d);
    // translate base sub-pols
    for ( i = j = 0 ; i < k ; j+=split[i++] ) ff_poly_translate (h+j, h+j, split[i]-1, c);
    d0 = split[0]; d1 = split[k-1];
    // compute (x+c)^d0 and (x+c)^d1
    _ff_set_one(G[0][d0]); for ( i = 0 ; i < d0 ; i++ ) _ff_set_zero(G[0][i]);
    ff_poly_translate_small (G[0],G[0],d0,c);
    if ( d0 < d1 ) ff_poly_mult_xpa (G[1],G[0],d0,c);

    for (;;) {
        k = k/2;
        for ( i = j = 0 ; j < k ; j++ ) {
            //assert ( split[2*j] == d0 || split[2*j] == d1 );
            ff_poly_mult (w, &e, G[split[2*j]-d0], split[2*j], h+i+split[2*j], split[2*j+1]-1);
            for ( ii = 0 ; ii < split[2*j] ; i++,ii++ ) _ff_addto (h[i],w[ii]);
            for ( ; ii <= e ; i++,ii++ ) _ff_set (h[i],w[ii]);
            split[j] = split[2*j]+split[2*j+1];
            for ( ; ii < split[j] ; i++,ii++ ) _ff_set_zero (h[i]);
        }
        if ( k == 1 ) break;
        if ( split[0] == d0+d0 ) ff_poly_square (G[0],G[0],d0);
        else ff_poly_mult(G[0],&d0,G[0],d0,G[1],d1);
        d0 = split[0]; d1 = split[k-1];
        if ( d0 < d1 ) ff_poly_mult_xpa (G[1],G[0],d0,c);
    }
    ff_poly_copy (g, 0, h, d);
    ff_poly_stack_pop (h);
    free (split);
}
 
// g(x)=f(x+c) using a single polynomial mult plus O(n) field mults and one inversion
// Uses Method F in GG-1997 "Fast algorithms for Taylor shifts and certain difference equations" (http://dl.acm.org/citation.cfm?doid=258726.258745), based on Aho et al (1975).
// Bit complexity is M(d*log p)+O((d+log(p)*M(log p)))
void ff_poly_translate_big (ff_t g[], ff_t f[], int d, ff_t c)
{
    register ff_t *u, *v, *w, *t, x, y;
    register int i, j;
    int e;
    
    if ( d < 2 ) { ff_poly_translate_small(g,f,d,c); return; }
    if ( d >= _ff_p )  { ff_poly_translate_rec(g,f,d,c); return; }
    u = ff_poly_stack_alloc (d);  v = ff_poly_stack_alloc(d); w = ff_poly_stack_alloc(2*d); t = ff_poly_stack_alloc(d);
    
    //set t[i] = 1/i! for i in [0..d]
    ff_inverted_factorials(t,d);
    
    // set u = sum_i (d-i)!*f[n-i]*x^i
    _ff_set(u[d],f[0]); _ff_set(u[d-1],f[1]); _ff_set_one(x); _ff_x2(x); _ff_set(y,x); _ff_mult(u[d-2],f[2],y);
    for ( i = 3, j = d-i ; i <= d ; i++,j-- ) { _ff_inc(x); _ff_multby(y,x); _ff_mult(u[j],f[i],y); }
        
    // set v = sum_i a^i/i!*x^i
    _ff_set_one(v[0]);  _ff_set(v[1],c); _ff_set(x,c);
    for ( i = 2 ; i <= d ; i++ ) { _ff_multby(x,c); _ff_mult(v[i],x,t[i]); }
        
    ff_poly_mult (w,&e,u,d,v,d);
    _ff_set(g[0],w[d]); _ff_set(g[1],w[d-1]);
    for ( i = 2, j = d-i ; i <= d ; i++,j-- ) _ff_mult(g[i],w[j],t[i]);
    
    ff_poly_stack_pop (u);
}

void ff_poly_translate (ff_t g[], ff_t f[], int d, ff_t c)
{
    ff_t t,w;

    if ( d < 0 ) return;
    switch ( d ) {
    case 0: _ff_set(g[0],f[0]); return;
    case 1: _ff_set(g[1],f[1]); _ff_mult(w,c,f[1]); _ff_add(g[0],f[0],w); return;
    case 2: _ff_set(g[2],f[2]); _ff_mult(w,c,f[2]); _ff_add(t,w,f[1]); _ff_x2(w); _ff_add(g[1],f[1],w); _ff_mult(t,t,c); _ff_add (g[0],f[0],t); return;
    case 3: _ff_set(g[3],f[3]); _ff_mult(w,c,f[3]); _ff_add(g[2],f[2],w); _ff_mult(t,g[2],c); _ff_add(g[1],f[1],t); _ff_addto(g[2],w); _ff_mult(t,g[1],c); _ff_add(g[0],f[0],t); _ff_mult(t,g[2],c); _ff_addto(g[1],t); _ff_addto(g[2],w); return;
    }
    if ( d <= FF_POLY_TRANSLATE_SMALL_DEGREE ) { ff_poly_translate_small (g,f,d,c); return; }
    if ( d < FF_POLY_TRANSLATE_BIG_DEGREE ) { ff_poly_translate_rec (g,f,d,c); return; }
    ff_poly_translate_big(g,f,d,c);
}

/*
    Given polys a, b, and c, solves the first-order linear differential equation af'+bf = c mod x^n,
    using the initial condition f(0)=0.  Based on the algorithm of Brent & Kung in Section 2.3 of BSMS2008
    f can alias any of the inputs.
*/

