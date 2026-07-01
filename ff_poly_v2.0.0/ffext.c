#include <stdlib.h>
#include <stdio.h>
#include "ff.h"
#include "ffext.h"
#include "ffpoly.h"
#include "ffpolysmall.h"
#include "ffpolyalloc.h"
#include "cstd.h"

/*
    Copyright 2008-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

ff_t _ff3_c;
ff_t _ff3_f[4];                                     // irred minimal poly of z - polynomial basis {1,z,z^2}  of the form x^3-s or x^3-x-s
ff_t _ff3_zp[3];                                    // if p=1mod3 this is just a multiple of z
ff_t _ff3_z2p[3];                                   // z^{2p} cached because it used by ff3_exp_p
int _ff3_trace_z2;                                  // trace of z^2 is 0 for p=1mod and 2 o.w., note that trace of z is always 0, store this as in int rather than an ff_t

static ff_t _ff2_3Sylow_tab[42][2][2];              // These values depend on s (the quadratic non-residue passed in) and are not reused, having a static table is simple a convenience
ff_t _ff2_cbrt_unity[2];
int _ff2_cbrt_setup;
int _ff2_nr_setup;
static ff_t _ff2_nr_a;                              // for p=1mod4, a=_ff2_nr_a is an alement of  F_p s.t. z+a is not a quadratic residue in Fp^2

void ff_ext_setup(void) { _ff3_f[3] = 0; _ff2_nr_setup = _ff2_cbrt_setup = 0;  }    // don't zero everything, just enough to detect uninitialized cases

// standard 4-ary exponentiation (fixed 2-bit window)  TODO: include x^p = bar(x) optimization
void ff2_exp_ui_s (ff_t o[2], ff_t a[2], unsigned long e, ff_t s)
{
    register int i;
    register unsigned long j, m;
    ff_t b[4][2], c[2];
    
//printf ("exp(%ldz+%ld, %ld)\n", _ff_get_ui(a[1]), _ff_get_ui(a[0]), e);
    switch (e) {
    case 0:  ff2_set_one (o);  return;
    case 1:  ff2_set(o,a); return;
    case 2:  ff2_square_s(o,a,s); return;
    }
    i = ui_highbit(e);
    if ( i&1 ) i--;
    m = 3UL<<i;
    ff2_set (b[1], a);
    ff2_square_s (b[2],b[1],s);
    ff2_mult_s(b[3],b[2],b[1],s);
    ff2_set (c, b[(m&e)>>i]);
    for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
        ff2_square_s(c,c,s);  ff2_square_s(c,c,s);
        j = (m&e)>>i;
        if ( j ) ff2_mult_s(c,c,b[j],s);
    }
    ff2_set (o, c);
//printf("expresult=%ldz+%ld\n", _ff_get_ui(o[1]), _ff_get_ui(o[0]));
}

void ff2_nonresidue (ff_t o[2])
{
    _ff_set_one (o[1]);
    if ( (_ff_p&3)==1 ) _ff_set(o[0],0);                            // if p=1 mod 4 then z is a non-residue (since N(z)=-s is non-residue in Fp)
    else if ( _ff2_nr_setup ) _ff_set (o[0], _ff2_nr_a);            // use known non-residue z+a if we have already computed it
    else {
        register ff_t t;
        
        ff_setup_2g();
        do {
            ff_random (&_ff2_nr_a);
            _ff_square (t, _ff2_nr_a);  _ff_subfrom (t, _ff_2g);    // compute N(z+a) = a^2-z^2 = a^2-_ff_2g
        } while ( ff_residue (t) );
        _ff2_nr_setup = 1;
        _ff_set (o[0], _ff2_nr_a);
    }
}


/*
    sqrt algorithm over F_p^2, reduces problem to two sqrts in F_p

    We are given an element of the form a1z+a0 where z^2=s=_ff_2g (a non-residue in F_p)
    If a1 is zero, we just compute sqrt(a0) in F_p, and if it doesn't exist, sqrt(a0)=sqrt(a0/s)*sqrt(s).

    We now assume a1!=0.  If (b1z+b0)^2 = (a1z+a0) then

           b0^2+sb1^2 = a0  and  2b0b1=a1

    and we know that both b0 and b1 are non-zero.  We then obtain

           b0^4 - a0b0^2 + a1^2s/4 = 0 =>  b0^2 = (a0 +/- sqrt(N(a)))/2
 
    and similarly

          sb1^4 - a0b1^2 + a1^2/4 = 0 => b1^2 = (a0 -/+ sqrt(N(a)))/(2s)

    Note that b0^2*b1^2 = a1^2/4 is a QR, but 1/s is not a QR, so exactly one of (-a0+sqrt(N(a)))/2 and (-a0-sqrt(N(a)))/2 is a QR

*/
int ff2_sqrt (ff_t o[2], ff_t a[2])
{
    ff_t x, y, t;
    
    ff_setup_2g();
    if ( _ff_zero(a[1]) ) {
        if ( ff_invsqrt(&t, a, 1) ) {
            _ff_mult (o[0], t, a[0]);
            _ff_set_zero(o[1]);
        } else {
            _ff_mult (o[1],t,a[0]);
            _ff_set_zero(o[0]);
        }
        return 1;
    }
    ff2_norm(&t,a);
    if ( ! ff_sqrt (&x, &t) ) return 0;  // a is a QR in F_p^2 iff N(a) is a QR in F_p
    _ff_add (t,a[0],x);
    _ff_div2(y,t);
    if ( ! ff_invsqrt (&t, &y, 1) ) {   
        _ff_mult(x,y,t);                // if r=sqrt(y) is not in F_p then b1 = r/sqrt(s)=r/z, we have r=y*tz, so b1=y*t
        ff_mult(t,t,_ff_2g);            // b0=a1/(2b1)=a1/(2yt)=a1z/(2ytz)=a1z/(2r)=a1zr/(2y)=a1ytz^2/(2y)=a1ts/2
        ff_mult(t,t,a[1]);
        _ff_div2(y,t);                  // got b0
    } else {
        ff_mult(y,y,t);                 // got b0
        ff_mult(t,t,a[1]);              // b1=a1/(2b0)=a1/(2yt)=a1t/2
        _ff_div2(x,t);                  // got b1
    }
#if ! FF_FAST
{
    ff_t b[2];
    _ff_set(b[0],y); _ff_set(b[1],x);
    ff2_square(b,b);
    if ( ! ff2_equal(b,a) ) { printf ("ff2_sqrt failed, (%lu+%luz)^2 != %lu+%luz\n",_ff_get_ui(b[0]), _ff_get_ui(b[1]), _ff_get_ui(a[0]), _ff_get_ui(a[1])); abort(); }
    }
#endif
    _ff_set(o[1],x);
    _ff_set(o[0],y);
    return 1;
}

void _ff2_setup_cbrt(void)
{
    register int i;
    register unsigned long n;
    ff_t t0,t1;
    ff_t r[2],s[2],t[2],w[2];
    
    if ( _ff2_cbrt_setup ) return;
    ff_setup_2g();
    if ( _ff_p1mod3 || _ff_p==3 ) { puts ("_ff2_setup_cbrt called when p!= 2mod 3, this should never happen.");  abort(); }
    _ff2_cbrt_setup = 1;
    // Note that p+1 = 3^e*m where m is not divisible by 3, and we have p-1 not div by 3

    if ( _ff_p3_e == 1 ) {
        // life is slightly easer in this situation, since (-1+sqrt(-3))/2 generates the Sylow 3-subgroup of F_p^2
        // unfortunately we still need to compute sqrt(-3) in the standard basis (-3 is not necessarily _ff_2g)
        _ff_set_i(t0,-3);
        if ( ff_sqrt_ext (&t1,&t0) ) { printf ("Impossible, -3 is a QR mod %ld with p=2mod3\n", _ff_p);  abort(); }
        _ff_div2(_ff2_3Sylow_tab[0][0][0], _ff_negone);
        _ff_div2(_ff2_3Sylow_tab[0][0][1], t1);
        ff2_square (_ff2_3Sylow_tab[0][1],_ff2_3Sylow_tab[0][0]);
        ff2_set (_ff2_cbrt_unity,_ff2_3Sylow_tab[0][0]);
//printf ("cube root of unity is %ldz+%ld, 2g=%d\n", _ff_get_ui(_ff2_cbrt_unity[1]), _ff_get_ui(_ff2_cbrt_unity[0]), _ff_get_ui(_ff_2g));
        return;
    }

    // To find a cubic non-residue, it suffices to find r^(m*3^(e-1)) not in F_p.   (when e is 1 this is just r^m)
    // To get into the 3-Sylow, we need to exponentiate by (p-1)m = (p+1)(m-1) + 2m(3^(e-1)-1), so we exponentiate by m-1 first.
    n = (_ff_p+1)/(3*_ff_p3_m)-1;                   // n = 3^(e-1)-1
    _ff_set_one(t[1]); _ff_set_one(t[0]);
    for(;;) { 
        ff2_exp_ui (r, t, _ff_p3_m-1);              // r = t^(m-1)
        ff2_mult (s, r, t);                         // s = t^m
        if ( ! _ff_zero(s[1]) ) {                   // if s is in F_p, t is a cubic residue
            ff2_exp_ui (w,s,n);                     // w = t^(m*(3^(e-1)-1))
            _ff_mult(t0,w[0],s[1]); _ff_mult(t1,w[1],s[0]);
            _ff_addto(t1,t0);                       // t1 = z coeff of t^(m*3^(e-1)) = t^((p+1)/3)
            if ( ! _ff_zero(t1) ) break;            // t is  a cubic residue iff t^((p+1)/3) is in F_p iff t1 is zero
        }
        _ff_inc(t[0]);
        if ( _ff_zero(t[0]) ) { _ff_inc(t[1]); _ff_set_one(t[0]); }
    }
    // t is a non CR, r=t^(m-1), s=t^m,  w=t^(m*(3^(e-1)-1))
    ff2_norm(&t0,r);                                // t0 = t^((p+1)(m-1))
    ff2_mult(s,s,w);                                // s = t^(3^(e-1)m)
    ff2_square(w,w);                                // w = t^(2*(3^(e-1)-1)m)
    ff2_mult(s,s,w);                                // s = t^((3^(e-1)+2(3^(e-1)-1))m) = t^((3^e-2)m)
    ff2_scalar_mult(r,t0,s);                        // r = t^((p+1)(m-1)+(3^e-2)m) = t^((p-1)m) generates the 3-Sylow of F_p^2  
    ff2_set (_ff2_3Sylow_tab[0][0], r);
    ff2_square (_ff2_3Sylow_tab[0][1],_ff2_3Sylow_tab[0][0]);
    for ( i = 1 ; i < _ff_p3_e ; i++ ) {
        ff2_mult(_ff2_3Sylow_tab[i][0],_ff2_3Sylow_tab[i-1][0],_ff2_3Sylow_tab[i-1][1]);    // tab[i][0] = tab[i-1][0]^3
        ff2_square(_ff2_3Sylow_tab[i][1], _ff2_3Sylow_tab[i][0]);                           // tab[i][1] = tab[i][0]^2
    }
    ff2_set (_ff2_cbrt_unity,_ff2_3Sylow_tab[_ff_p3_e-1][0]);
//printf ("cube root of unity is %ldz+%ld, 2g=%d\n", _ff_get_ui(_ff2_cbrt_unity[1]), _ff_get_ui(_ff2_cbrt_unity[0]), _ff_get_ui(_ff_2g));
}


// computes a^{-1/3} for a in the 3-Sylow subgroup, returns 0 if not a quadratic residue
// uses precomputed table of 2e powers of the 3-Sylow generator to reduce running time by a factor of 2 over standard Tonelli-Shanks (still O(e^2)).
int ff2_3Sylow_invcbrt (ff_t o[2], ff_t a[2])
{
    ff_t b[2], q[2], q1[2], t[2], w1[2], w2[2];
    register int i, j, k;
    
    ff2_setup_cbrt();

    // handle easy cases without setting up
    if ( _ff_one(a[0]) ) { _ff_set_one(o[0]);  return 1; }      // use 1 as the cube root of 1
    if ( _ff_p3_e == 1 ) return 0;
    
    // set w1 and w2 to the two elements of order 3 in the 3-Sylow (i.e. the two non-trivial cube roots of unity)
    ff2_set (w1, _ff2_3Sylow_tab[_ff_p3_e-1][0]);
    ff2_set (w2, _ff2_3Sylow_tab[_ff_p3_e-1][1]);
    ff2_set_one (t);
            
    ff2_set (b,a);
    j  = 0;     // avoid compiler warning
    do {
        ff2_set (q, b);
        for ( i = 1 ; i <= _ff_p3_e ; i++ ) {       // i<=e is just a safety check in case a isn't in the 3-Sylow, this could be removed
            j=1;
            if ( ff2_equal(q,w1) ) break;
            j=0;
            if ( ff2_equal(q,w2) ) break;
            ff2_set(q1,q);
            ff2_square (q,q);  ff2_mult(q,q,q1);
        }
        k = _ff_p3_e-i;
#if ! FF_FAST
        if ( k < 0 ) { printf ("Unexpected result: k<0 in ff2_3Sylow_invsqrt?!  a = %luz+%lu, p = %lu\n", _ff_get_ui(a[1]), _ff_get_ui(a[0]), _ff_p);  abort(); }
#endif
        if ( k <= 0 ) return 0;
        ff2_mult (b, b, _ff2_3Sylow_tab[k][j]);         // the product of all the elements S[k] we multiply into b here is b^{-1}, since we terminate with b=1
        ff2_mult (t, t, _ff2_3Sylow_tab[k-1][j]);       // multiply t by S[k-1]=cbrt(S[k]), the product of all these will be b^{-1/3}
    } while ( ! ff2_is_one(b) );
#if ! FF_FAST
    ff2_square (b, t); ff2_mult(b,b,t);
    ff2_mult(b,b,a);
    if ( ! ff2_is_one(b) ) { printf ("ff2_3Sylow_invcbrt failed, (%luz+%lu)^3 *  (%luz+%lu) != 1\n", _ff_get_ui(t[1]), _ff_get_ui(t[0]), _ff_get_ui(a[1]), _ff_get_ui(a[0])); abort(); }
#endif
    ff2_set(o,t);
    return 1;
}

void _ff3_setup (void)
{
    _ff_set_one(_ff3_f[3]);  _ff_set_zero(_ff3_f[2]);
    if ( _ff_p1mod3 ) {
        ff_setup_3g();
        _ff3_c = _ff_3g;
        ff3_set_zero(_ff3_zp);
        ff_exp_ui(_ff3_zp+1,&_ff3_c,(_ff_p-1)/3);
        _ff_set_zero(_ff3_f[1]);
        _ff_neg (_ff3_f[0], _ff3_c);
        _ff3_trace_z2 = 0;
    } else {
        _ff_set_one(_ff3_f[1]);
        ff_negate (_ff3_f[1]);
        _ff_set(_ff3_f[0],_ff3_f[1]);
        for(;;) {
            if ( ! ff_poly_roots_d3 (0, _ff3_f) ) break;
            _ff_dec(_ff3_f[0]);
        }
        _ff_neg (_ff3_c, _ff3_f[0]);
        ff3_zn_mod (_ff3_zp, _ff_p, _ff3_f);
        _ff3_trace_z2 = 2;
    }
    ff3_square (_ff3_z2p, _ff3_zp);
    // Note that the norm of z is the product of the roots of f, which is -f[0] = __ff3_c
    // the trace of z is -f[2] = 0, trace of z^2 is 0 for p=1mod3 and 2 o.w.
}


// multiplies modulo (z^3-rz-s), overlap ok (duplicates code in ff3_mult, then adjusts for r)
static inline void _ff3_mult_mod_rs (ff_t o[3], ff_t a[3], ff_t b[3], ff_t r, ff_t s)
{
    register ff_t s1,s2,s3,t0,t1,t2,w1,w2,w3;
    
    _ff_add(t1,a[0],a[1]);   _ff_add(t2,b[0],b[1]);  _ff_mult(s1,t1,t2);
    _ff_add(t1,a[0],a[2]);   _ff_add(t2,b[0],b[2]);  _ff_mult(s2,t1,t2);
    _ff_add(t1,a[1],a[2]);   _ff_add(t2,b[1],b[2]);  _ff_mult(s3,t1,t2);
    _ff_mult(t0,a[0],b[0]);  _ff_mult(t1,a[1],b[1]); _ff_mult(t2,a[2],b[2]);
    _ff_sub(w1,s3,t1);  _ff_subfrom(w1,t2); // w1 = (a1+a2)(b1+b2)-a1b1-a2b2 = a1b2+b2b1
    _ff_mult(w3,w1,s);                      // w3 = s(a1b2+b2b1)
    _ff_add(o[0],t0,w3);                    // o[0] = a0b0 + (a1b2+b2b1)s
    _ff_mult(w2,t2,s);                      // w2 = a2b2s
    _ff_sub(w3,s1,t0);  _ff_subfrom(w3,t1); // w3 = (a0+a1)(b0+b1)-a0b0-a1b1 = a0b1+b1a0
    _ff_add(o[1],w2,w3);                    // o[1] = a0b1+b1a0+a2b2s
    _ff_sub(w3,s2,t0);  _ff_subfrom(w3,t2); // w3 = (a0+a2)(b0+b2)-a0b0-a2b2 = a0b2+b2a0
    _ff_add(o[2],w3,t1);                    // o[2] = a0b2+b2a0+a1b1
    // z^3=rz+s so we need to add rz[(a1b2+b2b1)+a2b2z] - don't optimize for r=0 or 1 here, we assume general case applies
    _ff_mult(t0,r,w1);
    _ff_addto(o[1],t0);                     // o[1] = a0b1+b1a0+a2b2s + r(a1b2+b2b1)
    _ff_mult(t1,r,t2);
    _ff_addto(o[2],t1);                     // o[2] = a0b2+b2a0+a1b1 + ra2b2
    // 10M + 17A
}

// squares mod x^3-rx-s
static inline void _ff3_square_mod_rs (ff_t o[3], ff_t a[3], ff_t r, ff_t s)
{
    register ff_t s1,s2,t0,t1,t2,w1,w2,w3;

    _ff_mult(w1,a[2],r);                        // a2r
    _ff_add(s1,a[0],w1);                        // a0+a2r
    _ff_add(s2,a[1],a[1]);                      // 2a1
    _ff_mult(w3,s1,s2);                         // 2a0a1+2a1a2r
    _ff_mult(w2,a[2],s);                        // a2s
    ff_mult(s2,s2,w2);                          // 2a1a2s
    _ff_mult(t2,w2,a[2]);                       // a2^2s
    _ff_square(t0,a[0]); _ff_square(t1,a[1]);   // a0^2, a1^2
    _ff_addto(s1,a[0]);                         // 2a0+a2r
    _ff_mult(w1,s1,a[2]);                       // 2a0a2+a2^2r
    _ff_add(o[0],t0,s2);                        // a0^2+2a1a2s
    _ff_add(o[1],w3,t2);                        // 2a0a1+2a1a2+a2^2s
    _ff_add(o[2],w1,t1);                        // 2a0a2 + a2^2r + a1^2
    // 8M + 6A
}

// compute o=f(x) where x is in Fp^2 and f is in Fp[x]
void ff2_poly_eval_ff (ff_t o[2], ff_t f[], int d, ff_t x[2])
{
    ff_t t[2], y[2];
    register int i;

    if ( d < 0 ) { ff2_set_zero(o); return; }
    ff2_set_ff (y,f[d]);
    ff2_set(t,x);
    for ( i = d-1 ; i >= 0 ; i-- ) { ff2_mult(y,y,t);  _ff_addto(y[0], f[i]); }
    ff2_set(o,y);
    return;
}


// compute o=f(x) where x is in Fp^3 and f is in Fp[x]
void ff3_poly_eval_ff (ff_t o[3], ff_t f[], int d, ff_t x[3])
{
    ff_t t[3], y[3];
    register int i;

    ff_setup_3g();
    if ( d < 0 ) { ff3_set_zero (o); return; }
    ff3_set_ff (y, f[d]);
    ff3_set (t,x);
    for ( i = d-1 ; i >= 0 ; i-- ) { ff3_mult(y,y,t);  _ff_addto(y[0],f[i]); }
    ff3_set(o,y);
    return;
}


// compute o=f(x) where x is in Fp^2 and f is in Fp^2[x]
void ff2_poly_eval (ff_t o[2], ff_t f[], int d, ff_t x[2])
{
    ff_t t[2], y[2];
    register int i;

    if ( d < 0 ) { ff2_set_zero(o); return; }
    ff2_set (y, f+2*d);
    ff2_set(t,x);
    for ( i = d-1 ; i >= 0 ; i-- ) { ff2_mult(y,y,t);  ff2_add (y, y, f+2*i); }
    ff2_set(o,y);
    return;
}


// compute o=f(x) where x is in Fp^2 and f is in Fp^2[x]
void ff3_poly_eval (ff_t o[3], ff_t f[], int d, ff_t x[3])
{
    ff_t t[3], y[3];
    register int i;

    if ( d < 0 ) { ff3_set_zero(o); return; }
    ff3_set (y, f+3*d);
    ff3_set(t,x);
    for ( i = d-1 ; i >= 0 ; i-- ) { ff3_mult(y,y,t);  ff3_add (y, y, f+3*i); }
    ff3_set(o,y);
    return;
}


void ff3_minpoly (ff_t f[4], ff_t a[3])
{
    ff_t x[3];
    ff_t t;
    
    ff3_setup();
    _ff_set_one(f[3]);
    ff3_trace(&t,a);  _ff_neg(f[2],t);  // f[2] = -tr(a) = -(a+a^p+a^(p^2))
    ff3_exp_p (x,a);
    ff3_mult (x,x,a);               // x = a*a^p
    ff3_trace(f+1,x);               // f[1] = (a*a^p+a^p*a^(p^2)+a^(p^2)*a)
    ff3_norm(&t,a);  _ff_neg(f[0],t);   // f[0] = -N(a) = -a*a^p*a^(p^2)
}


// computes z^n mod f=z^3-az-b, optimized for a=1 case
// standard 4-ary exponentiation (fixed 2-bit window)
void ff3_zn_mod (ff_t o[3], unsigned long n, ff_t f[2])         // only looks at f[0] and f[1], implicitly assumes f[2]=0 and f[3]=1
{
    register int i;
    register unsigned long j, m, e;
    register ff_t a, b, w1,w2, w3;
    ff_t t[3];

    if ( ! n ) { ff3_set_one(o); return; }
    e=n;
    i = ui_highbit(e);
    if ( i&1 ) i--;
    m = 3UL<<i;
    j = (m&e)>>i;
    _ff_neg(a,f[1]);
    _ff_neg(b,f[0]);
//printf("a=%ld, f[1]=%ld, b=%ld, f[0]=%ld\n", _ff_get_ui(a), _ff_get_ui(f[1]), _ff_get_ui(b), _ff_get_ui(f[0]));
    _ff_set_zero(t[0]);  _ff_set_zero(t[1]);  _ff_set_zero(t[2]);
    switch(j){
    case 1: _ff_set_one(t[1]);  break;              // t = z
    case 2:  _ff_set_one(t[2]); break;              // t = z^2
    case 3: _ff_set(t[0],b);  _ff_set(t[1],a); break;       // t = z^3=az+b
    default:
        printf ("reached unreachable code in ff3_zn_mod\n");  abort();
    }
//printf("j=%d, %ldz^2+%ldz+%ld!\n", j, _ff_get_ui(t[2]), _ff_get_ui(t[1]), _ff_get_ui(t[0]));  
    if ( _ff_one(a) ) {
        for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
            _ff3_square_mod_1s(t,t,b);  _ff3_square_mod_1s(t,t,b);
//printf("i=%d, %ldz^2+%ldz+%ld!\n", i, _ff_get_ui(t[2]), _ff_get_ui(t[1]), _ff_get_ui(t[0]));  
            j = (m&e)>>i;
            switch(j) {
            case 1:
                _ff_mult(w1,t[2],b);
                _ff_set(w2,t[1]);
                _ff_add(t[1],t[0],t[2]);
                _ff_set(t[2],w2);
                _ff_set(t[0],w1);
                break;
            case 2:
                _ff_mult(w1,t[2],b);
                _ff_addto(t[2],t[0]);
                _ff_mult(w2,t[1],b);
                _ff_addto(t[1],w1);
                _ff_set(t[0],w2);
                break;
            case 3:
                _ff_mult(w1,t[2],b);
                _ff_mult(w2,t[1],b);
                _ff_addto(w2,t[0]);
                _ff_add(w3,t[0],t[2]);
                _ff_mult(t[0],w3,b);
                _ff_set(w3,t[2]);
                _ff_add(t[2],t[1],w1);
                _ff_add(t[1],w3,w2);
                break;
            }
            // 1.5M+1.75A on average, every 4 bits
        }
    } else {        
        for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
            _ff3_square_mod_rs(t,t,a,b);  _ff3_square_mod_rs(t,t,a,b);
//printf("i=%d, %ldz^2+%ldz+%ld!\n", i, _ff_get_ui(t[2]), _ff_get_ui(t[1]), _ff_get_ui(t[0]));  
            j = (m&e)>>i;
            switch(j) {
            case 1:
                _ff_mult(w1,t[2],b);
                _ff_mult(w2,t[2],a);
                _ff_set(t[2],t[1]);
                _ff_add(t[1],w2,t[0]);
                _ff_set(t[0],w1);
                break;
            case 2:
                _ff_mult(w1,t[2],b);
                _ff_mult(w2,t[2],a);
                _ff_add(t[2],t[0],w2);
                _ff_mult(w2,t[1],a);
                _ff_mult(t[0],t[1],b);
                _ff_add(t[1],w1,w2);
                break;
            case 3:
                _ff_mult(w1,t[2],b);
                _ff_mult(w2,t[2],a);
                _ff_add(w3,t[0],w2);
                _ff_mult(t[0],w3,b);
                _ff_mult(w2,w3,a);
                _ff_mult(w3,t[1],b);
                _ff_mult(t[2],t[1],a);
                _ff_add(t[1],w2,w3);
                _ff_addto(t[2],w1);
                break;
            }
            // 3M+1.5A on average, every 4 bits
//printf("j=%d, %ldz^2+%ldz+%ld!\n", j, _ff_get_ui(t[2]), _ff_get_ui(t[1]), _ff_get_ui(t[0]));  
        }
    }
    ff3_set (o, t);
}


// standard 4-ary exponentiation (fixed 2-bit window)
void ff3_exp_ui (ff_t o[3], ff_t a[3], unsigned long e)
{
    register int i;
    register unsigned long j, m;
    ff_t b[4][3], c[3];
    
    ff_setup_3g();
    switch (e) {
    case 0:  ff3_set_one (o);  return;
    case 1:  ff3_set(o,a); return;
    case 2:  ff3_square(o,a); return;
    }
    i = ui_highbit(e);
    if ( i&1 ) i--;
    m = 3UL<<i;
    ff3_set (b[1], a);
    ff3_square (b[2],b[1]);
    ff3_mult(b[3],b[2],b[1]);
    ff3_set (c, b[(m&e)>>i]);
    for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
        ff3_square(c,c);
        ff3_square(c,c);
        j = (m&e)>>i;
        if ( j ) ff3_mult(c,c,b[j]);
    }
    ff3_set (o, c);
}

// standard 4-ary exponentiation (fixed 2-bit window)
void ff3_exp_ui_rs (ff_t o[3], ff_t a[3], unsigned long e, ff_t r, ff_t s)
{
    register int i;
    register unsigned long j, m;
    ff_t b[4][3], c[3];
    
    ff_setup_3g();
    switch (e) {
    case 0:  ff3_set_one (o);  return;
    case 1:  ff3_set(o,a); return;
    case 2:  _ff3_square_mod_rs (o,a,r,s); return;
    }
    i = ui_highbit(e);
    if ( i&1 ) i--;
    m = 3UL<<i;
    ff3_set (b[1], a);
    _ff3_square_mod_rs (b[2],b[1],r,s);
    _ff3_mult_mod_rs(b[3],b[2],b[1],r,s);
    ff3_set (c, b[(m&e)>>i]);
    for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
        _ff3_square_mod_rs(c,c,r,s);  _ff3_square_mod_rs(c,c,r,s);
        j = (m&e)>>i;
        if ( j ) _ff3_mult_mod_rs(c,c,b[j],r,s);
    }
    ff3_set (o, c);
}

/*
   This algorithm is deterministic (given 2Sylow generator in F_p) - caller should flip coin to randomize roots
*/
int ff3_sqrt (ff_t o[3], ff_t a[3])
{
    ff_t x[3];
    ff_t b, c, d;

    if ( ff3_is_zero(a) ) { ff3_set_zero (o);  return 1; }
    if ( ff3_is_one(a) ) { ff3_set (o,a);  return 1; }

    ff3_setup();

    /*
        We need to compute a^n and a^(n+1)/2, where n is the odd part of p^3-1, which is m*(p^2+p+1), where m=2k+1 is the odd part of p-1
        The element a^n=N(a)^m is in the 2-sylow subgroup, and is also in Fp (Fp and Fp^3 have the same 2-Sylow subgroup).
        
        We can compute N(a) very efficiently, and then exponentiate by k in Fp, followed by a square and a multiply to obtain N(a)^(2k+1)=a^n.
        We use a standard Fp Tonelli-Shanks algorithm to compute the inverse of the square root of a^n, if it exists, and if it doesn't we detect this
        for less than the cost of a Legendre computation.
    
        Once we have a^(-n/2) in hand, we then compute a^((n+1)/2) to obtain a^(1/2).
    
        (n+1)/2 = ((2k+1)(p^2+p+1)+1)/2 = k(p^2+p+1) + p(p+1)/2 + 1
    
        We already have a^(k(p^2+p+1)) from above, so we need only compute a^((p+1)/2) followed by a Frobenius exponentiation by p, and two mults.
    */
    ff3_norm(&b,a);                                 // b = N(a) = a^(p^2+p+1)
    ff_exp_ui(&c,&b,(_ff_p2_m-1)/2);                // c = N(a)^k (save for later)
    _ff_square(d,c);                                // d = N(a)^(2k)
    ff_mult(b,b,d);                                 // b = N(a)^(2k+1) = N(a)^m is in the 2-Sylow
    if ( ! ff_2Sylow_invsqrt(&d,&b,0) ) return 0;   // d = 1/sqrt(b) = a^{-n/2)
    ff3_exp_ui(x,a,(_ff_p+1)/2);                    // x = a^((p+1)/2)
    ff3_exp_p(x,x);                                 // x = a^(p(p+1)/2)
    ff3_mult (x, x, a);                             // x = a^(p(p+1)/2+1)
    _ff_mult (b, d, c);                             // b = a^(k(p^2+p+1)-n/2)       (do Fp mults before scalar_mult)
    ff3_scalar_mult (x, b, x);                      // x = a^(k(p^2+p+1)+p(p+1)/2+1-n/2) = a^((n+1)/2-n/2) = a^(1/2)
    ff3_set (o, x); 
    return 1;
}

// computes (z+a)*v in F_p[z]/(z^3-rz-s), overlap ok
static inline void _ff3_mult_zpa_mod_rs (ff_t o[3], ff_t v[3], ff_t a, ff_t r, ff_t s)
{
    register ff_t t2,w1,w2;
    
    _ff_set(t2,v[2]);
    _ff_mult(w1,t2,a);
    _ff_add(o[2],v[1],w1);
    _ff_mult(w1,v[1],a);
    _ff_mult(w2,t2,r);
    _ff_addto(w1,w2);
    _ff_add(o[1],v[0],w1);
    _ff_mult(w1,v[0],a);
    _ff_mult(w2,t2,s);
    _ff_add(o[0],w1,w2);
    // 5M+4A
}


// Computes tr(sqrt(z)) in F_p^3=F_p[z]/(z^3-rz-s).  This is a support function for factoring quartics.
int ff3_trsqrt_zpa_mod_rs (ff_t o[1], ff_t a, ff_t r, ff_t s)
{
    ff_t u[3], v[3], w[3], x[3], f[4];
    ff_t b, c, d, t1, t2;

//printf("Computing tr(sqrt(z+a)) for F_p[z]/(z^3-rz-s) with p=%ld, a=%ld, r=%ld, s=%ld\n", _ff_p, _ff_get_ui(a), _ff_get_ui(r), _ff_get_ui(s));
    
    /*
        We use effectively the same algorithm as ff3_sqrt to compute sqrt(z), except we work in the user supplied basis
            and our input is (z+a) which simplifies a few things.
    */
    // N(z+a) = (z+a)(z+a)^p(z+a)^(p^2) = (z+a)(z^p+a)(z^(p^2)+a) = N(z)+a*tr(z^p*z)+a^2tr(z)+a^3 = s-ar+a^3  (since N(z)=s, tr(z)=0 and tr(z*z^p)=-r)
    _ff_mult(t1,a,r);                       // t1=ar
    _ff_sub(b,s,t1);                        // b=s-ar
    _ff_square(t2,a); ff_mult(t1,t2,a);         // t1=a^3
    _ff_addto(b,t1);                        // b=s-ar+a^3=N(z+a)
    ff_exp_ui(&c,&b,(_ff_p2_m-1)/2);            // c = N(z+a)^k (save for later)
    _ff_square(d,c);                        // d = N(z+a)^(2k)
    ff_mult(b,b,d);                         // b = N(z+a)^(2k+1) = N(z)^m is in the 2-Sylow
    if ( ! ff_2Sylow_invsqrt(&d,&b,0) ) return 0;   // d = 1/sqrt(b) = z^{-m(p^2+p+1)/2)

    // set modulus for ff3 mults
    _ff_set_one(f[3]);  _ff_set_zero(f[2]);
    _ff_set(f[1],r);  _ff_set(f[0],s);
    // combine computation of (z+a)^((p+1)/2) and z^p
    ff_poly_xpan_mod_d3 (w,a,(_ff_p-1)/2,f);        // w = (z+a)^((p-1)/2)
    _ff3_mult_zpa_mod_rs (x,w,a,r,s);           // x = (z+a)^((p+1)/2)
    _ff3_mult_mod_rs(u,w,x,r,s);                // u = (z+a)^p mod f
    _ff_subfrom(u[0],a);                        // u = (z+a)^p - a = z^p+a^p-a = z^p

    // x^p = x0+x1*z^p+x2*z^(2p) = x0+x1*u+x2*u^2
    ff3_scalar_mult(w,x[1],u);              // w = x1u
    _ff3_square_mod_rs(v,u,r,s);                // v = u^2
    ff3_scalar_mult(v,x[2],v);                  // v = x2u^2
    _ff_set(t1,x[0]);
    ff3_add(x,w,v);                         // x = x1u+x2u^2
    _ff_addto(x[0],t1);                     // x = x0+x1u+x2u^2 = (z^((p+1)/2))^p
    _ff3_mult_zpa_mod_rs (x,x,a,r,s);           // x = z^(p(p+1)/2+1)
    _ff_mult (b, d, c);                     // b = z^(k(p^2+p+1)-n/2)       (do F_p mults before scalar_mult)  recall that n=m(p^2+p+1) where m is odd part of p-1
    ff3_scalar_mult (x, b, x);                  // x = z^(k(p^2+p+1)+p(p+1)/2+1-n/2) = z^((n+1)/2-n/2) = z^(1/2)
#if ! FF_FAST
    _ff3_square_mod_rs (w, x,r,s);
    if ( !_ff_zero(w[0])||!_ff_one(w[1])||!_ff_zero(w[2]) ) { printf ("ff3_sqrt failed, (%luz^2+%luz+%lu)^2 = ", _ff_get_ui(x[2]), _ff_get_ui(x[1]), _ff_get_ui(x[0]));
                        printf ("%luz^^2+%luz+%luz != z",   _ff_get_ui(w[2]), _ff_get_ui(w[1]), _ff_get_ui(w[0])); }
#endif
    // now compute tr(x)=x0*tr(1)+x1*tr(z)+x2*(tr(z^2)) = 3x0 + 0 + 2rx2
    _ff_mult(t2,r,x[2]);
    _ff_add(t1,t2,t2);                      // t1=2rx2
    _ff_add(t2,x[0],x[0]);
    _ff_addto(t2,x[0]);                     // t2=3x0
    _ff_add(o[0],t1,t2);                        // o=tr(x)
                        
//printf("tr(z^(1/2)) = %ld\n", _ff_get_ui(o[0]));
    return 1;
}

static inline char *ff_monomial_str (char *buf, int d)
{
    if ( d <= 0 ) { buf[0] = '\0'; return buf; }
    if ( d == 1 ) { buf[0] = 'x'; buf[1] = '\0'; return buf; }
    sprintf (buf, "x^%d", d);
    return buf;
}

static inline char *ff2_poly_coeff_str (char *buf, ff_t a[2])
{
    if ( ff2_is_one(a) ) { buf[0] = '\0'; }
    else if ( _ff_zero(a[1]) ) { sprintf(buf, "%ld*", _ff_get_ui(a[0])); }
    else if ( _ff_zero(a[0]) ) { sprintf(buf, "%ld*a*", _ff_get_ui(a[1])); }
    else if ( _ff_one(a[1]) ) { sprintf(buf, "(a+%ld)*", _ff_get_ui(a[0])); }
    else { sprintf(buf, "(%ld*a+%ld)*", _ff_get_ui(a[1]), _ff_get_ui(a[0])); }
    return buf;
}

char *ff2_poly_sprint (char *buf, ff_t f[], int d)
{
    char *s;

    d = ff2_poly_degree (f, d);
    if ( d < 0 ) { strcpy (buf, "[0]"); return buf; }
    buf[0] = '[';
    s = buf+1;
    ff2_poly_coeff_str(s, f+2*d); s += strlen(s);
    ff_monomial_str(s, d); s += strlen(s);
    for ( int i = d-1 ; i >= 0 ; i-- ) {
        if ( ! ff2_is_zero(f+2*i) ) {
            memcpy (s, " + ", 3); s+= 3;
            ff2_poly_coeff_str(s, f+2*i); s += strlen(s);
            ff_monomial_str(s, i); s += strlen(s);
        }
    }
    if ( *(s-1) == '*' ) s--;
    else if ( *(s-1) == ' ' || *(s-1) == '[' ) *s++ = '1';
    sprintf (s, "] / (a^2 - %ld, %ld)", _ff_get_ui(_ff_2g), _ff_p);
    return buf;
}

void ff2_poly_print (ff_t f[], int d)
{
    char *buf;
    
    if ( d < 0 ) d = -1;
    buf = malloc (256+64*(d+1));
    puts (ff2_poly_sprint(buf,f,d));
    free (buf);
}

static inline char *ff3_poly_coeff_str (char *buf, ff_t a[3])
{
    if ( ff3_is_one(a) ) { buf[0] = '\0'; }
    else if ( _ff_zero(a[2]) && _ff_zero(a[1]) ) { sprintf(buf, "%ld*", _ff_get_ui(a[0])); }
    else if ( _ff_zero(a[2]) && _ff_zero(a[0]) ) { sprintf(buf, "%ld*a*", _ff_get_ui(a[1])); }
    else if ( _ff_zero(a[1]) && _ff_zero(a[0]) ) { sprintf(buf, "%ld*a^2*", _ff_get_ui(a[2])); }
    else if ( _ff_zero(a[2]) ) { sprintf(buf, "(%ld*a+%ld)*", _ff_get_ui(a[1]), _ff_get_ui(a[0])); }
    else if ( _ff_zero(a[1]) ) { sprintf(buf, "(%ld*a^2+%ld)*", _ff_get_ui(a[2]), _ff_get_ui(a[0])); }
    else if ( _ff_zero(a[0]) ) { sprintf(buf, "(%ld*a^2+%ld*a)*", _ff_get_ui(a[2]), _ff_get_ui(a[1])); }
    else { sprintf(buf, "(%ld*a^2+%ld*a+%ld)*", _ff_get_ui(a[2]), _ff_get_ui(a[1]), _ff_get_ui(a[0])); }
    return buf;
}

char *ff3_poly_sprint (char *buf, ff_t f[], int d)
{
    char *s;

    d = ff3_poly_degree (f, d);
    if ( d < 0 ) { strcpy (buf, "[0]"); return buf; }
    buf[0] = '[';
    s = buf+1;
    ff3_poly_coeff_str(s, f+3*d); s += strlen(s);
    ff_monomial_str(s, d); s += strlen(s);
    for ( int i = d-1 ; i >= 0 ; i-- ) {
        if ( ! ff3_is_zero(f+3*i) ) {
            memcpy (s, " + ", 3); s += 3;
            ff3_poly_coeff_str(s, f+3*i); s += strlen(s);
            ff_monomial_str(s, i); s += strlen(s);
        }
    }
    if ( *(s-1) == '*' ) s--;
    else if ( *(s-1) == ' ' || *(s-1) == '[' ) *s++ = '1';
    if ( _ff_p1mod3 ) 
        sprintf (s, "] / (a^3 - %ld, %ld)", _ff_get_ui(_ff3_c), _ff_p);
    else
        sprintf (s, "] / (a^3 - a - %ld, %ld)", _ff_get_ui(_ff3_c), _ff_p);
    return buf;
}

void ff3_poly_print (ff_t f[], int d)
{
    char *buf;
    
    if ( d < 0 ) d = -1;
    buf = malloc (256+96*(d+1));
    puts (ff3_poly_sprint(buf,f,d));
    free (buf);
}

// puts monic cubic in the form x^3+g1*x+g0 (note g2=0 and g3=1 are implicit, they are not set)
// overlap is ok
void ff2_poly_depress_monic_cubic (ff_t g[4], ff_t f[8])
{
    ff_t a[2], w0[2], w1[2];
    
    assert (_ff_p > 3);
    
    // compute g(x) = f(x+a) where a = - f2 / 3.
    ff2_scalar_mult (a, _ff_third, f+4);  ff2_negate (a);
    ff2_mult (w1, a, f+4); 
    ff2_square (w0,a);  ff2_add(w0,w0,w0);  ff2_sub(w0,f+2,w0);  ff2_mult (w0,a,w0);
    ff2_add (g+2, f+2, w1);                                             // g1 = f1 + a*f2
    ff2_add (g, f, w0);                                                 // g0 = f0 + a*(f1 - 2*a^2)
}

// puts monic cubic in the form x^3+g1*x+g0 (note g2=0 and g3=1 are implicit, they are not set)
// overlap is ok
void ff3_poly_depress_monic_cubic (ff_t g[6], ff_t f[12])
{
    ff_t a[3], w0[3], w1[3];

    if ( ff3_is_zero(f+6) ) { memmove(g,f,6*sizeof(ff_t)); return; }

    assert (_ff_p > 3);

    // compute g(x) = f(x+a) where a = - f2 / 3.
    ff3_scalar_mult (a, _ff_third, f+6);  ff3_negate (a);
    ff3_mult (w1, a, f+6); 
    ff3_square (w0,a);  ff3_add(w0,w0,w0);  ff3_sub(w0,f+3,w0);  ff3_mult (w0,a,w0);
    ff3_add (g+3, f+3, w1);                                             // g1 = f1 + a*f2
    ff3_add (g, f, w0);                                                 // g0 = f0 + a*(f1 - 2*a^2)
}

// Replaces monic f(x) with f(x-s) where s is chosen so that the coeff of x^(d-1) is zero.  Requires d to be nonzero mod p
void _ff2_poly_depress_monic_inplace (ff_t s[2], ff_t f[], int d)
{
    ff_t z, t[2], c[2];
    register int i, j;

    assert (d > 0);
    ff_invert_small_int(&z,d);
    ff_negate(z);
    _ff_mult(c[0],f[2*d-2],z);  _ff_mult(c[1],f[2*d-1],z);
    if ( s ) ff2_set(s,c);
    for ( i = d ; i ; i-- ) {
        for ( j = i ; j < d ; j++ ) { ff2_mult(t,c,f+2*j); ff2_add(f+2*(j-1),f+2*(j-1),t); }
        ff2_add(f+2*j-2,f+2*j-2,c);
    }
}

// Replaces monic f(x) with f(x-s) where s is chosen so that the coeff of x^(d-1) is zero.  Requires d to be nonzero mod p
void _ff3_poly_depress_monic_inplace (ff_t s[3], ff_t f[], int d)
{
    ff_t z, t[3], c[3];
    register int i, j;

    assert (d > 0);
    ff_invert_small_int(&z,d);
    ff_negate(z);
    _ff_mult(c[0],f[3*d-3],z);  _ff_mult(c[1],f[3*d-2],z);  _ff_mult(c[2],f[3*d-1],z);
    if ( s ) ff3_set(s,c);
    for ( i = d ; i ; i-- ) {
        for ( j = i ; j < d ; j++ ) { ff3_mult(t,c,f+3*j); ff3_add(f+3*(j-1),f+3*(j-1),t); }
        ff3_add(f+3*j-3,f+3*j-3,c);
    }
}

// computes g = f*fbar as a poly in Fp[x] of degree 2d, overlap not allowed
static inline void ff2_poly_monic_norm (ff_t g[], int *dg, ff_t f[], int d)
{
    ff_t t;

    assert (f!=g);
    memset (g, 0, 2*d*sizeof(ff_t));
    for ( int i = 0 ; i <= d ; i++ ) {
        ff2_norm (g+2*i,f+2*i);
        for ( int j = 0 ; j < i ; j++ ) {
            ff2_trace_pair(&t,f+2*i,f+2*j);
            _ff_addto (g[i+j],t);
        }
    }
    ff_poly_monic (g, dg, g, 2*d);
}


// g = f + fbar made monic 
// overlap OK
static inline void ff2_poly_monic_trace (ff_t g[], int *dg, ff_t f[], int d)
{
    for ( int i = 0 ; i <= d ; i++ ) _ff_set (g[i], f[2*i]);
    ff_poly_monic (g, dg, g, d);
}

// g = f - fbar made monic (if f is monic g will have lower degree)
// overlap OK
static inline void ff2_poly_monic_delta (ff_t g[], int *dg, ff_t f[], int d)
{
    for ( int i = 0 ; i <= d ; i++ ) _ff_add (g[i], f[2*i+1], f[2*i+1]);
    d = ff_poly_degree (g,d);
    ff_poly_monic (g,dg,g,d);
}

// table of precomputed values for ff2_poly_count_distinct_roots_d3 when p = 3
static char ff2_poly_roots_d3_mod3_tab[729] = {
1,1,1,1,1,1,1,1,1,3,3,3,0,0,0,0,0,0,3,0,0,3,0,0,3,0,0,3,0,0,0,3,0,0,0,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,0,0,0,0,3,0,3,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
2,3,0,0,1,1,0,1,1,2,3,0,0,1,1,0,1,1,3,0,2,1,1,0,1,1,0,1,1,0,3,0,2,1,1,0,3,0,2,1,1,0,1,1,0,1,0,1,1,0,1,0,2,3,1,1,0,1,1,0,3,0,2,3,0,2,1,1,0,1,1,0,1,0,1,0,2,3,1,0,1,
2,0,3,0,1,1,0,1,1,2,0,3,0,1,1,0,1,1,3,2,0,1,0,1,1,0,1,1,0,1,1,0,1,3,2,0,3,2,0,1,0,1,1,0,1,1,1,0,0,3,2,1,1,0,1,0,1,3,2,0,1,0,1,3,2,0,1,0,1,1,0,1,1,1,0,1,1,0,0,3,2,
2,0,0,0,1,1,3,1,1,3,1,1,2,0,0,0,1,1,2,0,0,0,1,1,3,1,1,1,1,3,0,0,2,1,1,0,1,0,1,1,3,1,0,2,0,3,1,1,2,0,0,0,1,1,1,3,1,0,2,0,1,0,1,1,1,0,1,1,3,0,0,2,3,1,1,2,0,0,0,1,1,
2,1,1,1,0,0,1,3,0,1,0,0,1,3,0,2,1,1,1,2,1,0,1,0,0,1,3,3,0,1,1,1,2,0,0,1,1,3,0,2,1,1,1,0,0,1,1,2,0,0,1,3,0,1,2,1,1,1,0,0,1,3,0,3,0,1,1,1,2,0,0,1,3,0,1,1,1,2,0,0,1,
2,1,1,1,0,0,1,0,3,1,0,0,1,0,3,2,1,1,1,1,2,0,0,1,0,3,1,2,1,1,1,0,0,1,0,3,3,1,0,1,2,1,0,1,0,3,1,0,1,2,1,0,1,0,3,1,0,1,2,1,0,1,0,1,0,3,2,1,1,1,0,0,1,2,1,0,1,0,3,1,0,
2,0,0,3,1,1,0,1,1,3,1,1,0,1,1,2,0,0,2,0,0,3,1,1,0,1,1,1,3,1,1,0,1,0,2,0,1,1,0,0,0,2,1,1,3,3,1,1,0,1,1,2,0,0,1,1,3,1,1,0,0,0,2,1,0,1,0,2,0,1,3,1,3,1,1,0,1,1,2,0,0,
2,1,1,1,3,0,1,0,0,1,0,0,2,1,1,1,3,0,1,2,1,0,1,3,0,1,0,2,1,1,1,3,0,1,0,0,3,0,1,0,0,1,1,1,2,3,0,1,0,0,1,1,1,2,3,0,1,0,0,1,1,1,2,1,3,0,1,0,0,2,1,1,1,1,2,3,0,1,0,0,1,
2,1,1,1,0,3,1,0,0,1,0,0,2,1,1,1,0,3,1,1,2,0,3,1,0,0,1,3,1,0,0,1,0,1,2,1,1,0,3,1,0,0,2,1,1,1,2,1,3,1,0,0,1,0,2,1,1,1,0,3,1,0,0,3,1,0,0,1,0,1,2,1,3,1,0,0,1,0,1,2,1
};

int ff2_poly_count_distinct_roots_d3 (ff_t f[8])
{
    ff_t g[8], t[2], w[2], *a, *b, t9;

    if ( ! ff2_is_one(f+6) ) f = ff2_poly_make_monic(ff2_poly_copy(g,0,f,3),3);
    if ( _ff_p == 3 ) return ff2_poly_roots_d3_mod3_tab[_ff_get_ui(f[0])+3*_ff_get_ui(f[1])+9*_ff_get_ui(f[2])+27*_ff_get_ui(f[3])+81*_ff_get_ui(f[4])+243*_ff_get_ui(f[5])];
    ff2_poly_depress_monic_cubic (g, f);
    
    a = g+2; b = g;     // g(x) = x^3 + ax + b

    // handle a = 0
    if ( ff2_is_zero(a) ) {
        if ( ff2_is_zero(b) ) return 1;
        ff2_negate(b);  return ff2_is_cube (b) ? 3 : 0;
    }
    
    // handle b = 0
    if ( ff2_is_zero(b) ) { ff2_negate(a);  return 2 + ff2_legendre(a); }

    // compute t = -D/3 = (4a^3+27b^2)/3
    ff2_square(t,a);  ff2_mult(t,t,a);  ff2_add(t,t,t);  ff2_add(t,t,t);  ff2_scalar_mult(t,_ff_third,t);   // t = 4a^3/3
    ff2_square(w,b);  _ff_set_ui(t9,9);  ff2_scalar_mult(w,t9,w);  ff2_add(t,t,w);                          // t = (4a^3+27b^2)/3
    
    if ( ff2_is_zero(t) ) return 2; // must be a double root because triple root would have a=b=0 handled above
    
    if ( ! ff2_sqrt (w,t) ) return 1;
    ff2_scalar_mult(w,_ff_third,w);  ff2_sub (w,w,b);  ff2_scalar_mult(w,_ff_half,w);                  // w = sqrt(t)/6 - b/2
    return ff2_is_cube(w) ? 3 : 0;
}

int ff2_poly_count_distinct_roots (ff_t f[], int d)
{
    ff_t *g, *h1, *h2;
    int dg, dh1, dh2;
    
    if ( d <= 1 ) return d;
    if ( d == 2 ) return ff2_poly_count_distinct_roots_d2 (f);
    if ( d == 3 ) return ff2_poly_count_distinct_roots_d3 (f);
    
    g = ff_poly_stack_alloc (2*d);
    h1 = ff_poly_stack_alloc (d);  h2 = ff_poly_stack_alloc(d);
    
    if ( ! ff2_is_one(f+2*d) || ! ff2_is_zero(f+2*d-2) ) {
        f = ff2_poly_make_monic(ff2_poly_copy(ff2_poly_stack_alloc(d),0,f,d),d);
        if ( (d%_ff_p) ) ff2_poly_depress_monic_inplace(0,f,d);
    }
    
    ff2_poly_monic_norm (g, &dg, f, d);
    ff_poly_gcd_xpn (g, &dg, 2, g, dg, (ui_highbit(_ff_p) > d ? 1 : 0));
    
    if ( ! dg ) { ff_poly_stack_pop(g); return 0; }
    
    ff2_poly_monic_trace (h1, &dh1, f, d);
    ff2_poly_monic_delta (h2, &dh2, f, d);
    ff_poly_gcd_reduce (h1, &dh1, h1, dh1, h2, dh2);
    // the irreducible factors of g that divide f are precisely those that divide h1 (with double the multiplicity)
    
    if ( dh1 ) ff_poly_gcd_reduce (h1, &dh1, g, dg, h1, dh1); // we need to do this to handle double roots correclty
    ff_poly_stack_pop(g);
    assert ((dg-dh1)%2 == 0);
    return dh1 + (dg-dh1)/2;
}

// computes g = f*fbar*fbar2 as a poly in Fp[x] of degree 3d, overlap not allowed
static inline void ff3_poly_monic_norm (ff_t g[], int *dg, ff_t f[], int d)
{
    ff_t *fbar, *fbar2;
    ff_t t, w[3];

    assert (f!=g);
    fbar = ff3_poly_stack_alloc(d);  fbar2 = ff3_poly_stack_alloc(d);
    for ( int i = 0 ; i <= d ; i++ ) ff3_exp_p(fbar+3*i,f+3*i);
    for ( int i = 0 ; i <= d ; i++ ) ff3_exp_p(fbar2+3*i,fbar+3*i);
    memset (g, 0, 3*d*sizeof(ff_t));
    for ( int i = 0 ; i <= d ; i++ ) {
        ff3_norm (g+3*i,f+3*i);                                      // iii
        for ( int j = 0 ; j < i ; j++ ) {
            ff3_mult (w,f+3*i,fbar+3*i);  ff3_mult (w,w,fbar2+3*j);  // ijj + jij + jji
            ff3_trace (&t,w);
            _ff_addto (g[2*i+j],t);
            ff3_mult (w,f+3*i,fbar+3*j);  ff3_mult (w,w,fbar2+3*j);  // iij + jii + iji
            ff3_trace (&t,w);
            _ff_addto (g[i+2*j],t);
            for ( int k = 0 ; k < j ; k++ ) {
                ff3_mult(w,f+3*i,fbar+3*j); ff3_mult(w,w,fbar2+3*k); // ijk + jki + kij
                ff3_trace (&t,w);
                _ff_addto(g[i+j+k],t);
                ff3_mult(w,f+3*i,fbar+3*k); ff3_mult(w,w,fbar2+3*j); // ikj + jik + kji
                ff3_trace (&t,w);
                _ff_addto(g[i+j+k],t);
            }
        }
    }
    ff_poly_stack_pop(fbar);
    ff_poly_monic (g, dg, g, 3*d);
}

// g = f + fbar + fbar2 made monic, as a poly in Fp[x] of degree d
// overlap OK
static inline void ff3_poly_monic_trace (ff_t g[], int *dg, ff_t f[], int d)
{
    for ( int i = 0 ; i <= d ; i++ ) ff3_trace(g+i,f+3*i);
    d = ff_poly_degree (g, d);
    ff_poly_monic (g, dg, g, d);
}

// g = f*fbar + fbar*fbar2 + fbar2*fbar made monic, as a poly in Fp[x] of degree d
// overlap OK
static inline void ff3_poly_monic_second_trace (ff_t g[], int *dg, ff_t f[], int d)
{
    ff_t *fbar, *h;
    int dh;
    
    fbar = ff3_poly_stack_alloc(d);
    for ( int i = 0 ; i <= d ; i++ ) ff3_exp_p(fbar+3*i,f+3*i);
    h = ff3_poly_stack_alloc(2*d);
    ff3_poly_mult(h,&dh,f,d,fbar,d);
    ff3_poly_monic_trace (g,dg,h,dh);
    ff_poly_stack_pop(fbar);
}

// currently we do not handle p = 3
int ff3_poly_count_distinct_roots_d3 (ff_t f[12])
{
    ff_t g[12], s[3], t[3], w[3], *a, *b, c;

    assert ( _ff_p > 3 );

    if ( ! ff3_is_one(f+9) ) f = ff3_poly_make_monic(ff3_poly_copy(g,0,f,3),3);
    ff3_poly_depress_monic_cubic (g, f);
    
    a = g+3; b = g;     // g(x) = x^3 + ax + b

    // handle a = 0
    if ( ff3_is_zero(a) ) {
        if ( ff3_is_zero(b) ) return 1;
        if ( ! _ff_p1mod3 ) return 1;
        ff3_negate(b);  return ff3_is_cube (b) ? 3 : 0;
    }
    
    // handle b = 0
    if ( ff3_is_zero(b) ) { ff3_negate(a);  return 2 + ff3_legendre(a); }
    
    // compute t = -D/3 = (4a^3+27b^2)/3
    ff3_square(t,a);  ff3_mult(t,t,a);  ff3_add(t,t,t);  ff3_add(t,t,t);  ff3_scalar_mult(t,_ff_third,t);   // t = 4a^3/3
    ff3_square(w,b);  _ff_set_ui(c,9);  ff3_scalar_mult(w,c,w);  ff3_add(t,t,w);                            // t = (4a^3+27b^2)/3
    if ( ff3_is_zero(t) ) return 2; // must be a double root because triple root would have a=b=0 handled above
    
    if ( _ff_p1mod3 ) {
        if ( ! ff3_sqrt (w,t) ) return 1;                                                                   // t is a square iff D is + Stickelberger
        ff3_scalar_mult(w,_ff_third,w);  ff3_sub (w,w,b);  ff3_scalar_mult(w,_ff_half,w);                   // w = sqrt(t)/6 - b/2
        return ff3_is_cube(w) ? 3 : 0;
    } else {
        ff_t x[3], ss[3], z[2], d, e;
        int sts;
        
        if ( ff3_residue (t) ) return 1;                                                                    // t a square iff b is not + Stickelberger
        _ff_mult(c,_ff_half,_ff_third); _ff_square(c,c); ff3_scalar_mult(s,c,t);                            // s = t/36
        ff3_scalar_mult(t,_ff_half,b);  ff3_negate(t);                                                      // t = t-b/2
        ff3_norm(&c,s);                                                                                     // c = N(s)
        ff3_exp_p(ss,s);  ff3_exp_ui (ss,ss,(_ff_p+1)/2);  ff3_mult(ss,ss,s);                               // ss = s^((p^2+p+2)/2) = sqrt(N(s))*sqrt(s))
        sts = ff_invsqrt (z+1, &c, 1);  assert (!sts);  _ff_mult(z[1],z[1],c);                              // z[1] = sqrt(N(s))
        ff3_norm(&d,t);  _ff_mult(d,d,c);                                                                   // d = N(s)*N(t)
        ff3_exp_p(w,ss); ff3_exp_p(w,w);  ff3_exp_p (x,t);  ff3_mult(w,w,x);  ff3_mult(w,w,ss);             // w = ss*t^p*ss^(p^2), x = t^p
        ff3_trace(&e,w); _ff_addto(d,e); _ff_mult(z[0],c,d);                                                // z[0] = N(s)(T(ss*t^p*ss^(p^2))+N(s)N(t))
        ff3_exp_p(w,x); ff3_mult(w,w,x); ff3_mult(w,w,ss); ff3_trace(&d,w);  _ff_mult(d,d,c);               // d = N(s)*T(ss*t^p*t^(p^2))
        ff3_norm(&e,ss);                                                                                    // e = N(ss)
        _ff_addto(d,e); _ff_mult(z[1],z[1],d);                                                              // z[1] = sqrt(N(s))*(N(ss)+N(s)*T(ss*t^p*t^(p^2)))
        return ff2_is_cube(z) ? 3 : 0;
    }
}


int ff3_poly_count_distinct_roots (ff_t f[], int d)
{
    ff_t *g1,*g2,*g3;
    int d1, d2, d3;
    
    if ( d <= 1 ) return d;
    if ( d == 2 ) return ff3_poly_count_distinct_roots_d2(f);
    if ( d == 3 && _ff_p > 3 ) return ff3_poly_count_distinct_roots_d3(f);

    g3 = ff_poly_stack_alloc(3*d);
    
    if ( ! ff2_is_one(f+3*d) || ! ff3_is_zero(f+3*(d-1)) ) {
        f = ff3_poly_make_monic(ff3_poly_copy(ff3_poly_stack_alloc(d),0,f,d),d);
        if ( (d%_ff_p) ) ff3_poly_depress_monic_inplace(0,f,d);
    }
    
    // Compute g3 as product of minpolys of distinct roots of any of f,fbar,fbar2
    ff3_poly_monic_norm (g3, &d3, f, d);
    ff_poly_gcd_xpn (g3, &d3, 3, g3, d3, (ui_highbit(_ff_p) > d ? 1 : 0));

    if ( ! d3 ) { ff_poly_stack_pop(g3); return 0; }
    
    // Compute g2 as product of minpolys of distinct roots of any two of f,fbar,fbar2 (these are also counted in d1)
    g2 = ff_poly_stack_alloc(2*d);
    ff3_poly_monic_second_trace (g2, &d2, f, d);
    ff_poly_gcd_reduce (g2, &d2, g2, d2, g3, d3);
    if ( ! d2 ) { assert ( d3%3 == 0 ); ff_poly_stack_pop(g3); return d3/3; }
    
    // compute g1 as product of minpolys of distinct roots one of all three of f,fbar,fbar2 (these are also counted in d2 and d1)
    g1 = ff_poly_stack_alloc(d);
    ff3_poly_monic_trace (g1, &d1, f, d);
    ff_poly_gcd_reduce (g1, &d1, g1, d1, g2, d2);

    // factors only in d2 or d3 must have degree 3, degree 1 factors in d1 get triple counted, other factors of d1 have degree 3, so d1+d2+d3 = 0 mod 3
    // factors in d2 and d3 should be double counted (they correspon to roots a and abar with abar2 not a factor), and they are 6/3 = 2
    // degree 3 factors in d1, d2, d3 should be triple counted (they correspond to roots a, abar, abar2), and they are 9/3 = 3
    // degree 1 factors in d1, d2, d3 should only get counted once, and they are 3/3 = 1
    assert ((d1+d2+d3)%3 == 0);
    ff_poly_stack_pop(g3);
    return (d1+d2+d3)/3;
}



