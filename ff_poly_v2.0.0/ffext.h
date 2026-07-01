#ifndef _FFEXT_INCLUDE_
#define _FFEXT_INCLUDE_

/*
    Copyright 2008-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "ff.h"
#include "ffpoly.h"
#include "ffpolyalloc.h"

#ifdef __cplusplus
extern "C" {
#endif

// subset of finite field arithmetic supported in degree 2 and 3 extension fields
// elements are represented as polys over F_p[x]/(g(x)) where g(x)=x^2-_ff_2g for degree 2
// For degree 3 minpoly is stored in _ff3_f as a depressed monic cubic
// For p = 1 mod 3 we use x^3 - c, where c=_ff_3g is the generator of the 3-Sylow subgroup of Fp
// for p != 1 mod 3 we use x^3 - x - c where c is chosen in ff3_setup() and stored in _ff3_c

extern int _ff2_cbrt_setup;
extern ff_t _ff2_cbrt_unity[2];
extern ff_t _ff3_zp[3];
extern ff_t _ff3_z2p[3];
extern ff_t _ff3_f[4];
extern ff_t _ff3_c;     // negation of constant coeff of _ff3_f, equal to _ff_3g when p is 1 mod 3

void ff_ext_setup(void);
void _ff3_setup(void);
static inline void ff3_setup(void) { if ( ! _ff3_f[3] ) _ff3_setup(); }
static inline void ff2_setup(void) { if ( ! _ff_2g ) ff_setup_2g(); }


// for convenience
static inline ff_t *ff2_random(ff_t o[2]) { _ff_random(o[0]);  _ff_random_nz(o[1]); return o; }
static inline ff_t *ff3_random(ff_t o[3]) { do { _ff_random(o[0]);  _ff_random(o[1]);  _ff_random(o[2]); } while ( _ff_zero(o[1]) && _ff_zero(o[2]) ); return o; }
static inline ff_t *ff2_set_zero(ff_t o[2]) { _ff_set_zero(o[0]); _ff_set_zero(o[1]); return o; }
static inline ff_t *ff3_set_zero(ff_t o[3]) { _ff_set_zero(o[0]); _ff_set_zero(o[1]); _ff_set_zero(o[2]); return o; }
static inline ff_t *ff2_set_one(ff_t o[2]) { _ff_set_one(o[0]); _ff_set_zero(o[1]); return o; }
static inline ff_t *ff3_set_one(ff_t o[3]) { _ff_set_one(o[0]); _ff_set_zero(o[1]); _ff_set_zero(o[2]); return o; }
static inline ff_t *ff2_set_ui (ff_t o[2], unsigned long x) { _ff_set_ui(o[0],x); _ff_set_zero(o[1]); return o; }
static inline ff_t *ff3_set_ui (ff_t o[3], unsigned long x) { _ff_set_ui(o[0],x); _ff_set_zero(o[1]); _ff_set_zero (o[2]); return o; }
static inline int ff2_is_zero(ff_t o[2]) { return ( _ff_zero(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_is_zero(ff_t o[3]) { return ( _ff_zero(o[0]) && _ff_zero(o[1]) &&  _ff_zero(o[2]) ); }
static inline int ff2_is_one(ff_t o[2]) { return ( _ff_one(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_is_one(ff_t o[3]) { return ( _ff_one(o[0]) && _ff_zero(o[1]) && _ff_zero(o[2]) ); }
static inline int ff2_neg_one(ff_t o[2]) { return ( _ff_neg_one(o[0]) && _ff_zero(o[1]) ); }
static inline int ff2_pm_one(ff_t o[2]) { return ( _ff_pm_one(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_neg_one(ff_t o[3]) { return ( _ff_neg_one(o[0]) && _ff_zero(o[1]) &&  _ff_zero(o[2]) ); }
static inline int ff3_pm_one(ff_t o[3]) { return ( _ff_pm_one(o[0]) && _ff_zero(o[1]) && _ff_zero(o[2]) ); }
static inline void ff2_dec(ff_t o[2]) { _ff_dec(o[0]); }
static inline void ff3_dec(ff_t o[2]) { _ff_dec(o[0]); }
static inline void ff2_inc(ff_t o[2]) { _ff_inc(o[0]); }
static inline void ff3_inc(ff_t o[2]) { _ff_inc(o[0]); }
static inline void ff2_next(ff_t o[2]) { _ff_inc(o[0]); if ( _ff_zero(o[0]) ) _ff_inc(o[1]); }
static inline void ff3_next(ff_t o[3]) { _ff_inc(o[0]); if ( _ff_zero(o[0]) ) { _ff_inc(o[1]); if ( _ff_zero(o[1])  ) _ff_inc(o[2]); } }
static inline ff_t *ff2_set(ff_t o[2], ff_t a[2]) { _ff_set(o[0],a[0]);  _ff_set(o[1],a[1]); return o; }
static inline ff_t *ff3_set(ff_t o[3], ff_t a[3]) { _ff_set(o[0],a[0]);  _ff_set(o[1],a[1]);  _ff_set(o[2],a[2]); return o; }
static inline ff_t *ff2_set_ff(ff_t o[2], ff_t a) { _ff_set(o[0],a);  _ff_set_zero(o[1]); return o; }
static inline ff_t *ff3_set_ff(ff_t o[3], ff_t a) { _ff_set(o[0],a);  _ff_set_zero(o[1]);  _ff_set_zero(o[2]); return o; }
static inline int ff2_equal(ff_t o[2], ff_t a[2]) { return (_ff_equal(a[0],o[0])&&_ff_equal(a[1],o[1])); }
static inline int ff3_equal(ff_t o[3], ff_t a[3]) { return (_ff_equal(a[0],o[0])&&_ff_equal(a[1],o[1])&&_ff_equal(a[2],o[2])); }

static inline ff_t *ff2_negate(ff_t a[2]) { ff_negate(a[0]); ff_negate(a[1]); return a;}
static inline ff_t *ff3_negate(ff_t a[3]) { ff_negate(a[0]); ff_negate(a[1]); ff_negate(a[2]); return a;}
static inline ff_t *ff2_neg(ff_t o[2], ff_t a[2]) { if ( o == a ) ff2_negate(o); else {_ff_neg(o[0],a[0]); _ff_neg(o[1],a[1]); } return o; }
static inline ff_t *ff3_neg(ff_t o[3], ff_t a[3]) { if ( o == a ) ff3_negate(o); else {_ff_neg(o[0],a[0]); _ff_neg(o[1],a[1]); _ff_neg(o[2],a[2]); } return o; }

// overlap ok (we avoid _ff_add(x,y,x) and _ff_sub(x,y,x)
static inline ff_t *ff2_add(ff_t o[2], ff_t a[2], ff_t b[2]) { if ( o == b ) { _ff_addto (o[0],a[0]); _ff_addto(o[1],a[1]); } else {_ff_add(o[0],a[0],b[0]);  _ff_add(o[1],a[1],b[1]); } return o; }
static inline ff_t *ff3_add(ff_t o[3], ff_t a[3], ff_t b[3]) { if ( o == b ) {_ff_addto (o[0],a[0]); _ff_addto(o[1],a[1]); _ff_addto(o[2],a[2]); } else { _ff_add(o[0],a[0],b[0]);  _ff_add(o[1],a[1],b[1]);  _ff_add(o[2],a[2],b[2]); } return o; }
static inline ff_t *ff2_sub(ff_t o[2], ff_t a[2], ff_t b[2]) { if ( o == b ) { _ff_subfrom(o[0],a[0]); _ff_subfrom(o[1],a[1]); ff2_negate(o); } else { _ff_sub(o[0],a[0],b[0]);  _ff_sub(o[1],a[1],b[1]); } return o; }
static inline ff_t *ff3_sub(ff_t o[3], ff_t a[3], ff_t b[3]) { if ( o == b ) { _ff_subfrom(o[0],a[0]); _ff_subfrom(o[1],a[1]); _ff_subfrom(o[2],a[2]);  ff3_negate(o); } else { _ff_sub(o[0],a[0],b[0]);  _ff_sub(o[1],a[1],b[1]);  _ff_sub(o[2],a[2],b[2]); } return o; }

static inline ff_t *ff2_dbl(ff_t o[2], ff_t a[2]) { return ff2_add(o,a,a); }
static inline ff_t *ff3_dbl(ff_t o[2], ff_t a[2]) { return ff3_add(o,a,a); }

// compute N(a) in the field Fp[z]/(z^2-s) (typically s=_ff_2g is our 2-Sylow generator)  3M+2A
static inline void ff2_norm_s (ff_t o[1], ff_t a[2], ff_t s)
{
    register ff_t t1,t2;
    
    _ff_square(t1,a[1]);
    _ff_neg(t2,s);                              // it is worth the negation in order to use sum_2_mults
    _ff_sum_2_mults (o[0],a[0],t1,t2,a[0]);     // N(a) = a[0]^2 - s*a[1]^2
}
static inline void ff2_norm (ff_t o[1], ff_t a[2]) { ff2_norm_s(o,a,_ff_2g); }

static inline void ff2_trace (ff_t o[1], ff_t a[2])
    { ff_add(o[0],a[0],a[0]); }
    
static inline ff_t *ff2_exp_p (ff_t o[2], ff_t a[2])
    { _ff_set (o[0], a[0]); _ff_neg (o[1], a[1]); return o; }
static inline ff_t *ff2_bar (ff_t o[2], ff_t a[2])
    { _ff_set (o[0], a[0]); _ff_neg (o[1], a[1]); return o; }

static inline ff_t *ff2_scalar_mult (ff_t o[2], ff_t c, ff_t a[2]) { ff_mult(o[0],c,a[0]); ff_mult(o[1],c,a[1]); return o; }

// multiplication in Fp[z]/(z^2-s) 5M+2A
static inline ff_t *ff2_mult_s (ff_t o[2], ff_t a[2], ff_t b[2], ff_t s) // compute (a[1]z+a[0])(b[1]z+b[0]) mod (z^2-s)
{
    register ff_t t1;

    _ff_mult(t1,a[1],b[1]);
    _ff_sum_2_mults(o[1],a[0],a[1],b[0],b[1]);          // o[1] = a[0]b[1]+a[1]b[0]
    _ff_sum_2_mults(o[0],a[0],s,t1,b[0]);               // o[0] = a[0]^2+a[1]b[1]s
    return o;
}
static inline ff_t *ff2_mult (ff_t o[2], ff_t a[2], ff_t b[2]) { return ff2_mult_s(o,a,b,_ff_2g); }

// o = tr(a*bbar) 3M+3A
static inline void ff2_trace_pair (ff_t o[1], ff_t a[2], ff_t b[2])
{
    ff_t t1,t2;
    
    _ff_neg(t2,b[1]);
    _ff_mult (t1, t2, _ff_2g);
    _ff_sum_2_mults (o[0], a[0], a[1], t1, b[0]);
    _ff_x2(o[0]);
}

// multiplies (b[1]z+b[0])*(z+a) mod (z^2-s)  3M+2A
static inline ff_t *ff2_mult_zpa_s (ff_t o[2], ff_t b[2], ff_t a, ff_t s)
{
    register ff_t t1;
    
    _ff_mult(t1,a,b[1]); _ff_addto(t1,b[0]);
    _ff_sum_2_mults(o[0],b[0],b[1],s,a);
    _ff_set(o[1],t1);
    return o;
}

// squaring in F[z]/(z^2-s)  4M+2A
static inline ff_t *ff2_square_s (ff_t o[2], ff_t a[2], ff_t s)
{
    register ff_t t1;

    _ff_square(t1,a[1]);
    _ff_dbl_mult(o[1],a[0],a[1]);
    _ff_sum_2_mults(o[0],a[0],s,t1,a[0]);
    return o;
}

static inline ff_t *ff2_square (ff_t o[2], ff_t a[2]) { return ff2_square_s(o,a,_ff_2g); }

// inversion in F[z]/(z^2-s) I+5M+3A
static inline ff_t *ff2_invert_s (ff_t o[2], ff_t a[2], ff_t s)
{
    ff_t t;

    ff2_norm_s (&t, a, s);
    _ff_invert (t, t);
    ff2_scalar_mult (o, t, a);
    ff_negate (o[1]);
    return o;
}
static inline ff_t *ff2_invert (ff_t o[2], ff_t a[2]) { return ff2_invert_s(o,a,_ff_2g); }

void ff2_exp_ui_s (ff_t o[2], ff_t a[2], unsigned long e, ff_t s);  // computes in F_p[z]/(z^2-s)
void ff2_exp_ui_s_x (ff_t o[2], ff_t a[2], unsigned long e, ff_t s);
static inline void ff2_exp_ui (ff_t o[2], ff_t a[2], unsigned long e) { ff2_exp_ui_s(o,a,e,_ff_2g); }


// Cube roots in F_p^2 are only relevant when p=2mod3, since otherwise the 3-Sylow of F_p^2 is the 3-Sylow of F_p.
void _ff2_setup_cbrt(void);
static inline void ff2_setup_cbrt(void) { if ( ! _ff2_cbrt_setup ) _ff2_setup_cbrt(); }
int ff2_3Sylow_invcbrt (ff_t o[2], ff_t a[2]);

// g(x)=f(x+c), works in place, currently d_f must be <= FF_BINOMIAL_MAX
void ff2_poly_translate (ff_t g[], int *pd_g, ff_t f[], int d_f, ff_t c[2]);


// Note ff_setup_3g() must be called before using any of the ff3 functions below
// ff3_poly_eval, ff3_sqrt, and ff3_exp_ui do this automatically

static inline ff_t *ff3_scalar_mult (ff_t o[3], ff_t c, ff_t a[3]) { ff_mult(o[0],c,a[0]); ff_mult(o[1],c,a[1]); ff_mult(o[2],c,a[2]); return o; }

// overlap ok, 11M + 6A (8A)
static inline ff_t *ff3_mult (ff_t o[3], ff_t a[3], ff_t b[3])
{
    register ff_t t0,t1,t2,w1,w2,w3;

    _ff_sum_2_mults(w1,a[1],a[2],b[1],b[2]);
    _ff_sum_2_mults(t0,a[0],_ff3_c,w1,b[0]);
    _ff_mult(w2,a[2],b[2]);
    _ff_sum_2_mults(w3,a[0],a[1],b[0],b[1]);
    _ff_mult(t2,w2,_ff3_c); 
    _ff_add(t1,t2,w3);
    _ff_sum_2_mults(w3,a[0],a[2],b[0],b[2]);
    _ff_mult(t2,a[1],b[1]);
    _ff_add(o[2],w3,t2);
    _ff_set(o[1],t1);
    _ff_set(o[0],t0);
    if ( ! _ff_p1mod3 ) { _ff_addto(o[1],w1); _ff_addto(o[2],w2); }
    return o;
}


// squares mod x^3-s, 7M + 5A
static inline void _ff3_square_mod_0s (ff_t o[3], ff_t a[3], ff_t s)
{
    register ff_t s1,s2,t0,t1,t2,w1,w2,w3;

    _ff_add(s1,a[1],a[1]);                      // 2a1
    _ff_mult(w1,s1,a[0]);                       // 2a0a1
    _ff_mult(s2,_ff3_c,a[2]);                   // a2s
    _ff_mult(w2,s1,s2);                         // 2a1a2s
    _ff_mult(t2,s2,a[2]);                       // a2^2s
    _ff_square(t0,a[0]);  _ff_square(t1,a[1]);  // a0^2, a1^2
    _ff_mult(w3,a[0],a[2]);                     // a0a2
    _ff_x2(w3);                                 // 2a0a2
    _ff_add(o[0],t0,w2);                        // a0^2 + 2a1a2s
    _ff_add(o[1],w1,t2);                        // 2a0a1 + a2^2s
    _ff_add(o[2],w3,t1);                        // 2a0a2 + a1^2
}

// squares mod x^3-x-s, 7M + 6A
static inline void _ff3_square_mod_1s (ff_t o[3], ff_t a[3], ff_t s)
{
    register ff_t s1,s2,t0,t1,t2,w1,w2,w3;

    _ff_add(s1,a[0],a[2]);                      // a0+a2
    _ff_square(w1,s1);                          // a0^2 + 2a0a2 + a2^2
    _ff_add(s2,a[1],a[1]);                      // 2a1
    _ff_mult(w3,s1,s2);                         // 2a0a1 + 2a1a2
    _ff_mult(w2,a[2],s);                        // a2s
    ff_mult(s2,s2,w2);                          // 2a1a2s
    _ff_mult(t2,w2,a[2]);                       // a2^2s
    _ff_square(t0,a[0]); _ff_square(t1,a[1]);   // a0^2, a1^2
    _ff_subfrom(w1,t0);                         // 2a0a2 + a2^2
    _ff_add(o[0],t0,s2);                        // a0^2 + 2a1a2s
    _ff_add(o[1],w3,t2);                        // 2a0a1 + 2a1a2 + a2^2s
    _ff_add(o[2],w1,t1);                        // 2a0a2 + a2^2 + a1^2
}

// overlap ok, 7M + 5A (6A)
static inline ff_t *ff3_square (ff_t o[3], ff_t a[3])
{
    register ff_t t0,t1,t2,t3,t4;
    
    _ff_set(t0,a[0]);                           // t0 = a0
    if ( ! _ff_p1mod3 ) _ff_addto(t0,a[2]);     // t0 = a0 or a0 + a2
    _ff_add (t1, a[1], a[1]);                   // t1 = 2a1
    _ff_mult (t2, a[2], _ff3_c);                // t2 = a2s
    _ff_set (t3,a[0]);                          // copy a0 to handle overlap
    _ff_sum_2_mults(o[0],a[0],t1,t2,a[0]);      // a0^2 + 2a1a2s
    _ff_set (t4,a[1]);                          // copy a1 to handle overlap
    _ff_sum_2_mults(o[1],t0,t2,a[2],t1);        // 2a1a0 + a2^2s or 2a1(a0+a2) + a2^2s
    _ff_addto(t0,t3);                           // t0 = 2a0 or 2a0 + a2 
    _ff_sum_2_mults(o[2],t0,t4,t4,a[2]);        // 2a0a2 + a1^2 or (2a0+a2)a2 + a1^2
    return o;
}

// exponentiating by p (the Frobenius map) is very fast, potentially only 2M and at most 6M+4A, faster than ff3_mult by a lot
// computes a^p = a[0]+a[1]z^p+a[2]z^{2p} using 2M or 6M+4A, overlap ok
static inline ff_t *ff3_exp_p (ff_t o[3], ff_t a[3])
{
    register ff_t w0,w1,w2;
    
    if ( _ff_p1mod3 ) {
        // if p=1mod3 we know z^p is a multiple of z and z^2p is a multiple of z^2
        _ff_set(o[0],a[0]);
        ff_mult(o[1],a[1],_ff3_zp[1]);
        ff_mult(o[2],a[2],_ff3_z2p[2]);
    } else {
        _ff_sum_2_mults(w0,a[1],a[2],_ff3_z2p[0],_ff3_zp[0]);
        _ff_sum_2_mults(w1,a[1],a[2],_ff3_z2p[1],_ff3_zp[1]);
        _ff_sum_2_mults(w2,a[1],a[2],_ff3_z2p[2],_ff3_zp[2]);
        _ff_add(o[0],a[0],w0);
        _ff_set(o[1],w1);
        _ff_set(o[2],w2);
    }
    return o;
}


extern int _ff3_trace_z2;

// tr(a[0]+a[1]z+a[2]z^2 = 3a[0]+a[1]tr(z)+a[2]tr(z^2) = 3a[0]+a[2]tr(z^2) since tr(z)=0 for z^3-z-s=0 and z^3-s=0
static inline void ff3_trace (ff_t o[1], ff_t a[3])
{
    register ff_t t1,t2;
    
    _ff_add(t1,a[0],a[0]);
    _ff_add(o[0],a[0],t1);
    if ( ! _ff3_trace_z2 ) return;
    _ff_add(t2,a[2],a[2]);          // we rely on the fact that tr(z^2) is either 0 or 2
    _ff_addto(o[0],t2); 
}

// computes a*a^p*a^(p^2) using 9M+6A or 9M+9A
static inline void ff3_norm (ff_t o[1], ff_t a[3])
{
    ff_t t1,t2,t3,t4,t5,t6;
    
    if ( _ff_p1mod3 ) { // Fp^3 = F_p[x]/(x^3-c) where c=_ff3_c
        _ff_add(t1,a[0],a[0]); _ff_addto (t1,a[0]); ff_negate(t1);  // t1 = -3*a0
        _ff_sum_2_mults (t3,a[1],a[2],t1,a[1]);                     // t3 = a1^2-3*a0*a2
        _ff_mult (t1,_ff3_c, a[2]); _ff_square (t2, a[2]);          // t1 = c*a2, t2 = a2^2
        _ff_sum_2_mults (t4,t3,t2,t1,a[1]);                         // t4 = c*a2^3 + a1^3 - 3*a0*a1*a2
        _ff_square (t1, a[0]);                                      // t1 = a0^2
        _ff_sum_2_mults (o[0],t1,_ff3_c,t4,a[0]);                   // o = a0^3 + c*(c*a2^3 + a1^3 - 3*a0*a1*a2)
    } else { // Fp^3 = Fp[x]/(x^3-x-c) where c = _ff3_c
        _ff_square(t1,a[1]); _ff_square(t2,a[2]);                   // t1 = a1^2, t2 = a2^2
        _ff_subfrom(t1,t2);                                         // t1 = a1^2-a2^2
        _ff_mult(t4,_ff3_c,a[2]);                                   // t4 = c*a2
        _ff_add(t3,t4,t4); _ff_addto(t3,t4); ff_negate(t3);         // t3 = -3*c*a2
        _ff_add(t6,a[0],a[2]); _ff_addto(t6,a[2]);                  // t6 = a0+2*a2
        _ff_sum_2_mults(t5,a[0],a[1],t3,t6);                        // t5 = a0*(a0+2*a2) - 3*c*a1*a2
        _ff_subfrom(t5,t1);                                         // t5 = a0*(a0+2*a2) - (a1^2-a2^2) - 3*c*a1*a2
        _ff_sum_2_mults(t3,a[1],t2,t4,t1);                          // t3 = a1*(a1^2-a2^2) + c*a2^3
        _ff_sum_2_mults (o[0],a[0],_ff3_c,t3,t5);                   // o = a0*(a0*(a0+2*a2)-(a1^2-a2^2)-3*c*a1*a2) + c*(a1*(a1^2-a2^2) + c*a2^3)
    }
}

// I+27M+12A or I+35M+25A
static inline ff_t *ff3_invert (ff_t o[3], ff_t a[3])
{
    ff_t t, b[3];
    
    ff3_norm (&t, a);
    _ff_invert(t,t);
    ff3_exp_p (b,a); 
    ff3_exp_p (o,b);
    ff3_mult (b,b,o);
    ff3_scalar_mult (o,t,b);
    return o;
}


void ff3_minpoly (ff_t f[4], ff_t a[3]);

// these functions evaluate a poly defined over F_p at a point in F_p^2 or F_p^3
void ff2_poly_eval_ff (ff_t y[2], ff_t f[], int d, ff_t x[2]);
void ff3_poly_eval_ff (ff_t y[3], ff_t f[], int d, ff_t x[3]);

// g(x)=f(x+c), works in place, currently d_f must be <= FF_BINOMIAL_MAX
void ff2_poly_translate (ff_t g[], int *pd_g, ff_t f[], int d_f, ff_t c[2]);


// these functions evaluate a poly defined over F_p^2 or F_p^3, computing y = f(x)
// the array f holds all the coefficients, ordered by degree; a total of (2*(d+1) entries in Fp^2, 3*(d+1) entries in Fp^3)
void ff2_poly_eval (ff_t y[2], ff_t f[], int d, ff_t x[2]);
void ff3_poly_eval (ff_t y[3], ff_t f[], int d, ff_t x[3]);

void ff2_nonresidue (ff_t nr[2]);
int ff2_sqrt (ff_t o[2], ff_t a[2]);
int ff3_sqrt (ff_t o[3], ff_t a[3]);

static inline int ff2_residue (ff_t a[2])
    { ff_t t;  if ( _ff_zero(a[1]) ) return 1;  ff2_norm(&t,a); return ff_residue (t); }
static inline int ff3_residue (ff_t a[3])
    { ff_t t;  ff3_norm(&t,a); return ff_residue (t); }
static inline int ff2_legendre (ff_t a[2])
    { ff_t t;  if ( ff2_is_zero(a) ) return 0; if ( _ff_zero(a[1]) ) return 1;  ff2_norm(&t,a); return ff_residue (t) ? 1 : -1; }
static inline int ff3_legendre (ff_t a[3])
    { ff_t t;  if ( ff3_is_zero(a) ) return 0; ff3_norm(&t,a); return ff_residue (t) ? 1 : -1; }
    
static inline int ff2_is_cube (ff_t a[2])
{
    ff_t b[2];
    
    if ( _ff_p == 3 || ff2_is_zero(a) ) return 1;
    if ( _ff_p1mod3 ) { ff2_norm(b,a); ff_exp_ui(b,b,(_ff_p-1)/3); return _ff_one(b[0]); } else { ff2_exp_ui(b,a,(_ff_p+1)/3); return _ff_zero(b[1]); }
}

static inline int ff3_is_cube (ff_t a[3])
{
    ff_t b[1];
    
    if ( ! _ff_p1mod3 || ff3_is_zero(a) ) return 1;
    ff3_norm(b,a); ff_exp_ui(b,b,(_ff_p-1)/3); return _ff_one(b[0]);
}

// finds all roots of a monic quadratic, returns number of roots found (0 or 2)
// f and r may overlap, note that each element of Fp^2 occupies 2 positions in these arrays (so first root is r[1]*z+r[0] and second is r[3]*z+r[2])
static inline int ff2_poly_roots_d2 (ff_t r[4], ff_t f[4])
{
    ff_t w[2];
    
    if ( ff2_is_zero (f+2) ) { ff2_neg (w,f);  if ( ! ff2_sqrt (r, w) ) return 0;  ff2_neg (r+2,r); return 2; }
    // w = f[1]^2 - 4*f[0]
    ff2_add (w, f, f); ff2_add (w, w, w); 
    ff2_square (r, f+2);  ff2_sub (w, r, w);
    if ( ! ff2_sqrt (w, w) ) return 0;
    ff2_sub (r, w, f+2);
    ff2_scalar_mult (r, _ff_half, r);
    ff2_sub (r+2, r, w);
    return 2;
}

// f need not be monic
static inline int ff2_poly_count_roots_d2 (ff_t f[6])
{
    ff_t r[2], w[2];
    
    ff2_add (w, f, f); ff2_add (w, w, w); ff2_mult (w, w, f+4);
    ff2_square (r, f+2);  ff2_sub (w, r, w);
    return ( ff2_residue (w) ? 2 : 0 );
}

// f need not be monic
static inline int ff2_poly_count_distinct_roots_d2 (ff_t f[6])
{
    ff_t r[2], w[2];
    
    ff2_add (w, f, f); ff2_add (w, w, w); ff2_mult (w, w, f+4);
    ff2_square (r, f+2);  ff2_sub (w, r, w);
    return 1 + ff2_legendre(w);
}

int ff2_poly_count_distinct_roots_d3 (ff_t f[8]);

// finds all roots of a monic quadratic, returns number of roots found (0 or 2)
// f and r may overlap, note that each element of Fp^3 occupies 3 positions in these arrays
static inline int ff3_poly_roots_d2 (ff_t r[6], ff_t f[6])
{
    ff_t w[3];
    
    if ( ff3_is_zero (f+3) ) { ff3_neg (w,f);  if ( ! ff3_sqrt (r, w) ) return 0;  ff3_neg (r+3,r); return 2; }
    // w = f[1]^2 - 4*f[0]
    ff3_add (w, f, f); ff3_add (w, w, w); 
    ff3_square (r, f+3);  ff3_sub (w, r, w);
    if ( ! ff3_sqrt (w, w) ) return 0;
    ff3_sub (r, w, f+3);
    ff3_scalar_mult (r, _ff_half, r);
    ff3_sub (r+3, r, w);
    return 2;
}

// f need not be monic
static inline int ff3_poly_count_roots_d2 (ff_t f[9])
{
    ff_t r[3], w[3];
    
    ff3_add (w, f, f); ff3_add (w, w, w); ff3_mult (w, w, f+6);
    ff3_square (r, f+3);  ff3_sub (w, r, w);
    return ( ff2_residue (w) ? 2 : 0 );
}

// f need not be monic
static inline int ff3_poly_count_distinct_roots_d2 (ff_t f[9])
{
    ff_t r[3], w[3];
    
    ff3_add (w, f, f); ff3_add (w, w, w); ff3_mult (w, w, f+6);
    ff3_square (r, f+3);  ff3_sub (w, r, w);
    return 1 + ff3_legendre(w);
}


static inline ff_t *ff2_poly_alloc (int d) { return ff_poly_alloc(2*d+1); } // this will allocate 2*d+2 ff_t's, which is what we want
static inline ff_t *ff3_poly_alloc (int d) { return ff_poly_alloc(3*d+2); } // this will allocate 3*d+3 ff_t's, which is what we want
static inline ff_t *ff2_poly_stack_alloc (int d) { return ff_poly_stack_alloc(2*d+1); } // this will allocate 2*d+2 ff_t's, which is what we want
static inline ff_t *ff3_poly_stack_alloc (int d) { return ff_poly_stack_alloc(3*d+2); } // this will allocate 3*d+3 ff_t's, which is what we want

static inline int ff2_poly_degree (ff_t f[], int d) { while ( d >= 0 && ff2_is_zero(f+2*d) ) d--; return d; }
static inline int ff3_poly_degree (ff_t f[], int d) { while ( d >= 0 && ff3_is_zero(f+3*d) ) d--; return d; }

static inline int ff2_poly_from_ui (ff_t f[], unsigned long F[], int d) { for ( int i = 0 ; i < 2*d+2 ; i++ ) _ff_set_ui (f[i],F[i]); return ff2_poly_degree(f,d); }
static inline int ff2_poly_from_si (ff_t f[], long F[], int d) { for ( int i = 0 ; i < 2*d+2 ; i++ ) _ff_set_si (f[i],F[i]); return ff2_poly_degree(f,d); }
static inline int ff2_poly_from_mpz (ff_t f[], mpz_t F[], int d) { for ( int i = 0 ; i < 2*d+2 ; i++ ) _ff_set_mpz (f[i],F[i]); return ff2_poly_degree(f,d); }
static inline int ff3_poly_from_ui (ff_t f[], unsigned long F[], int d) { for ( int i = 0 ; i < 3*d+3 ; i++ ) _ff_set_ui (f[i],F[i]); return ff3_poly_degree(f,d); }
static inline int ff3_poly_from_si (ff_t f[], long F[], int d) { for ( int i = 0 ; i < 3*d+3 ; i++ ) _ff_set_si (f[i],F[i]); return ff3_poly_degree(f,d); }
static inline int ff3_poly_from_mpz (ff_t f[], mpz_t F[], int d) { for ( int i = 0 ; i < 3*d+3 ; i++ ) _ff_set_mpz (f[i],F[i]); return ff3_poly_degree(f,d); }


static inline ff_t *ff2_poly_copy (ff_t g[], int *dg, ff_t f[], int df)
    { memcpy (g, f, 2*(df+1)*sizeof(ff_t)); if ( dg ) *dg = df; return g; }
static inline ff_t *ff3_poly_copy (ff_t g[], int *dg, ff_t f[], int df)
    { memcpy (g, f, 3*(df+1)*sizeof(ff_t)); if ( dg ) *dg = df; return g; }

static inline int ff2_poly_random (ff_t f[], int d)
    { register int i;  for ( i = 0 ; i <= d ; i++ ) ff2_random(f+2*i); return ff2_poly_degree(f,d); }
static inline void ff2_poly_random_monic (ff_t f[], int d)   // d must be >= 0
    { register int i; assert(d>=0);  for ( i = 0 ; i < d ; i++ ) ff2_random(f+2*i); ff2_set_one(f+2*d); return; }
static inline void ff2_poly_random_depressed_monic (ff_t f[], int d) // d must be >= 1
    { register int i;  assert(d>0); for ( i = 0 ; i < d-1 ; i++ ) ff2_random(f+2*i); ff2_set_one(f+2*d); ff2_set_zero(f+2*(d-1)); return; }

static inline int ff3_poly_random (ff_t f[], int d)
    { register int i;  for ( i = 0 ; i <= d ; i++ ) ff3_random(f+3*i); return ff3_poly_degree(f,d); }
static inline void ff3_poly_random_monic (ff_t f[], int d)   // d must be >= 0
    { register int i; assert(d>=0);  for ( i = 0 ; i < d ; i++ ) ff3_random(f+3*i); ff3_set_one(f+3*d); return; }
static inline void ff3_poly_random_depressed_monic (ff_t f[], int d) // d must be >= 1
    { register int i;  assert(d>0); for ( i = 0 ; i < d-1 ; i++ ) ff3_random(f+3*i); ff3_set_one(f+3*d); ff2_set_zero(f+3*(d-1)); return; }

static inline ff_t *ff2_poly_make_monic (ff_t f[], int d)
{
    ff_t c[2];

    if ( ff2_is_one(f+2*d) ) return f;
    ff2_invert (c, f+2*d);  ff2_set_one (f+2*d);
    for ( d-- ; d >= 0 ; d-- ) ff2_mult (f+2*d, f+2*d, c);
    return f;
}

static inline ff_t *ff3_poly_make_monic (ff_t f[], int d)
{
    ff_t c[3];

    if ( ff3_is_one(f+3*d) ) return f;
    ff3_invert (c, f+3*d);  ff3_set_one (f+3*d);
    for ( d-- ; d >= 0 ; d-- ) ff3_mult (f+3*d, f+3*d, c);
    return f;
}

// computes the discriminant -27*f0^2*f3^2 + 18*f0*f1*f2*f3 - 4*f0*f2^3 - 4*f1^3*f3 + f1^2*f2^2 of a general cubic
// = 9*f0*f3(-3*f0*f3+2*f1*f2) + (-4*f0*f2)*f2*2 + f1^2*(f2^2-4*f1*f3)
static inline ff_t *ff2_cubic_disc (ff_t d[2], ff_t f[8])
{
    ff_t three; _ff_triple(three, _ff_mont_R);
    ff_t a1[2]; ff2_mult(a1,f,f+3*2); ff2_scalar_mult(a1,three,a1);             // 3*f0*f3
    ff_t b1[2]; ff2_mult(b1,f+1*2,f+2*2); ff2_dbl(b1,b1); ff2_sub(b1,b1,a1);    // -3*f0*f3+2*f1*f2
    ff2_scalar_mult(a1,three,a1);                                               // 9*f0*f3
    ff_t a2[2]; ff2_negate(ff2_mult(a2,f,f+2*2)); ff2_dbl(a2,ff2_dbl(a2,a2));   // -4*f0*f2
    ff_t b2[2]; ff2_square(b2,f+2*2);                                           // f2^2
    ff_t a3[2]; ff2_square(a3,f+1*2);                                           // f1^2
    ff_t b3[2]; ff2_mult(b3,f+1*2,f+3*2); ff2_dbl(b3,ff2_dbl(b3,b3));           // 4*f1*f3
    ff2_sub(b3,b2,b3);                                                          // f2^2-4*f1*f3
    ff2_mult(a1,a1,b1); ff2_mult(a2,a2,b2); ff2_mult(a3,a3,b3);
    ff2_add(d,a1,ff2_add(a2,a2,a3));
    return d;
}

void _ff2_poly_depress_monic_inplace (ff_t s[2], ff_t f[], int d_f);
static inline void ff2_poly_depress_monic_inplace (ff_t s[2], ff_t f[], int d_f)
    { if ( d_f > 0 && ff2_is_zero(f+2*d_f-2) ) { if ( s ) ff2_set_zero(s); return; } _ff2_poly_depress_monic_inplace (s,f,d_f); }
static inline ff_t *ff2_poly_depress_monic (ff_t s[1], ff_t g[], ff_t f[], int d_f)
    { ff2_poly_depress_monic_inplace(s,ff2_poly_copy(g,0,f,d_f),d_f); return g; }
void ff2_poly_depress_monic_cubic (ff_t g[4], ff_t f[8]);

void _ff3_poly_depress_monic_inplace (ff_t s[3], ff_t f[], int d_f);
static inline void ff3_poly_depress_monic_inplace (ff_t s[3], ff_t f[], int d_f)
    { if ( d_f > 0 && ff3_is_zero(f+3*d_f-3) ) { if ( s ) ff3_set_zero(s); return; } _ff3_poly_depress_monic_inplace (s,f,d_f); }
static inline ff_t *ff3_poly_depress_monic (ff_t s[1], ff_t g[], ff_t f[], int d_f)
    { ff3_poly_depress_monic_inplace(s,ff3_poly_copy(g,0,f,d_f),d_f); return g; }
void ff3_poly_depress_monic_cubic (ff_t g[6], ff_t f[12]);

int ff2_poly_count_distinct_roots_d2 (ff_t f[6]);
int ff3_poly_count_distinct_roots_d2 (ff_t f[9]);
int ff2_poly_count_distinct_roots_d3 (ff_t f[8]);
int ff3_poly_count_distinct_roots_d3 (ff_t f[12]);
int ff2_poly_count_distinct_roots (ff_t f[], int d); // f should have 2*(d+1) entries
int ff3_poly_count_distinct_roots (ff_t f[], int d); // f should have 3*(d+1) entries

int ff3_trsqrt_zpa_mod_rs (ff_t o[1], ff_t a, ff_t r, ff_t s);          // Computes tr(sqrt(z)) in F_p^3=F_p[z]/(z^3-rz-s).  This is a support function for factoring quartics.

void ff3_exp_ui (ff_t o[3], ff_t a[3], unsigned long e);
void ff3_exp_ui_rs (ff_t o[3], ff_t a[3], unsigned long e, ff_t r, ff_t s);
// inversion not currently supported/needed in degree 3

// computes z^n mod f, where z is our generator for F_p^3
void ff3_zn_mod (ff_t o[3], unsigned long n, ff_t f[2]);            // only looks at f[0] and f[1], assumes f[2]=0 and f[3]=1


// overlap not allowed
static inline void ff2_poly_mult (ff_t h[], int *dh, ff_t f[], int df, ff_t g[], int dg)
{
    *dh = df+dg;  memset (h, 0, 2*(*dh+1)*sizeof(h[0]));
    for ( int i = 0 ; i <= df ; i++ ) for ( int j = 0 ; j <= dg ; j++ ) { ff_t w[2]; ff2_mult (w, f+2*i, g+2*j);  ff2_add (h+2*(i+j),h+2*(i+j), w); }
}

// overlap not allowed
static inline void ff3_poly_mult (ff_t h[], int *dh, ff_t f[], int df, ff_t g[], int dg)
{
    *dh = df+dg;  memset (h, 0, 3*(*dh+1)*sizeof(h[0]));
    for ( int i = 0 ; i <= df ; i++ ) for ( int j = 0 ; j <= dg ; j++ ) { ff_t w[3]; ff3_mult (w, f+3*i, g+3*j);  ff3_add (h+3*(i+j),h+3*(i+j), w); }
}

static inline int ff2_poly_squarefree (ff_t f[], int d)
    { err_printf ("function ff2_poly_squarefree is not yet implemented\n"); abort(); }

char *ff2_poly_sprint (char *buf, ff_t f[], int d);
void ff2_poly_print (ff_t f[], int d);

char *ff3_poly_sprint (char *buf, ff_t f[], int d);
void ff3_poly_print (ff_t f[], int d);

#ifdef __cplusplus
}
#endif

#endif
