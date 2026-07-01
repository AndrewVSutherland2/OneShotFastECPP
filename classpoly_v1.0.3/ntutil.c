#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include "ntutil.h"
#include "cstd.h"


// standard 4-ary exponentiation (fixed 2-bit window)
// p must fit in 32-bits
static unsigned long ui_powm_ui (unsigned long a, unsigned long e, unsigned long p)
{
    register unsigned long c, m;
    unsigned long b[4];
    register int i, j;
     
    assert ( p <= UINT_MAX );
    if ( ! e ) return 1;
    i = ui_highbit(e);
    if ( i&1 ) i--;
    m = 3UL<<i;
    b[1] = a%p;
    b[2] = (b[1]*b[1])%p;
    b[3] = (b[1]*b[2])%p;
    c = b[(m&e)>>i];
    for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
        c = (c*c)%p; 
        c = (c*c)%p;
        j = (m&e)>>i;
        if ( j ) c = (c*b[j])%p;
    }
    return c;
}

// hardwired table of mod-35 residue classes of odd primes p for which 5 or 7 is a non-residue mod p
static char _nr35_tab[35] =    { 0, 0, 5, 5, 0, 7, 7, 5, 5, 0, 7, 0, 5, 5, 0, 0, 0, 5, 5, 7, 7, 0, 5, 5, 7, 0, 7, 5, 5, 0, 0, 7, 5, 5, 7};

// Algorithm 1.5.1 in [CANT] - standard Tonelli-Shanks n^2 algorithm.  p must be an odd prime and this is not verified.
// Note that ff_sqrt is faster, but this code is useful if you only need one sqrt and p is smallish.  currently only supports p<2^31
// TODO: switch to mm.h to speed up modular reductions and handle p > INT_MAX
long i_sqrt_modp (long a, long p)
{
    register long q, x, y, z, b;
    register int i, r, m;

    assert ( p <= INT_MAX );
    if ( p == 2 ) return a&1;
    a %= p;  if ( a < 0 ) a += p;
    if ( ! a ) return 0;
    if ( a==1 ) return 1;
    // Get a non-residue, use hardwired tests to catch all but 1/32 of the cases
    if ( (p&3)==3 ) {
        x = p-1;                        // -1 is  a non-residue whenever p = 3 mod 4
    } else if ( (p&7)==5 ) {
        x = 2;                          // 2 is a non-residue whenever p = 5 mod 8
    } else if ( (p%3)==2 ) {
        x = p-3;                        // -3 is a non-residue whenever p = 2 mod 3
    } else {
        x = _nr35_tab[p%35];            // use 5 or 7 when possible
        if ( ! x ) for ( x = 11 ; ui_legendre(x,p) >= 0 ; x += 2 );
    }
    q = (p-1)>>1;
    for ( r = 1 ; !(q&1) ; r++ ) q >>= 1;
    y = ui_powm_ui (x,q,p);
    q = (q-1)>>1;
    z = ui_powm_ui (a,q,p);
    b = (a*z)%p;
    x= b;
    b = (b*z)%p;
    while ( b != 1 ) {
        q = b;
        for ( m = 1 ; m <= r; m++ ) {
            q = (q*q)%p;
            if ( q==1 ) break;
        }
        assert ( m <= r );
        if ( m == r ) return -1;
        q = y;
        for ( i = 0 ; i < r-m-1 ; i++ ) q = (q*q)%p;
        y = (q*q)%p;
        r = m;
        x = (x*q)%p;
        b = (b*y)%p;
    }
    // sanity ccheck, cant be removed to improve performance (very slightly)
    assert ( (x*x)%p == a );
    return x;
}

// computes a square-root of a modulo q=p^e provided a != 0 mod p and p > 2
long i_sqrt_modq (long a, long p, long q)
{
    long x, z, x1, z1;
    
    assert ( p != 2 && (a%p) != 0 );
    x = i_sqrt_modp (a, p);
    if ( x < 0 ) return -1;
    if ( q == p ) return x;
    z = ui_inverse (2*x, p);
    if ( q > INT_MAX ) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
        __int128 pp, xx, zz, x1, z1;
#pragma GCC diagnostic pop
        xx = x; zz = z; 
        for ( pp = p ; pp < q ; ) {
            assert ((xx*xx-a)%pp == 0);
            if ( pp > INT_MAX ) {
                while ( pp*pp > q ) pp /= p;
                if ( pp*pp < q ) pp *= p;
                xx %= pp; zz %= pp;
            }
            pp = pp*pp;
            x1 = (xx - (xx*xx-a)*zz) % pp; if ( x1 < 0 ) x1 += pp;
            z1 = (2LL*zz*(1LL - ((x1*zz) % pp))) % pp; if ( z1 < 0 ) z1 += pp;
            xx = x1; zz = z1;
            assert ((xx*xx-a)%pp == 0);
        }
        xx %= q;
        assert ((xx*xx-a)%q == 0);
        return (long)xx;
    } else {
        long pp;
        for ( pp = p ; pp < q ; ) {
            pp = pp*pp;
            x1 = (x - (x*x-a)*z) % pp; if ( x1 < 0 ) x1 += pp;
            z1 = (2*z*(1 - ((x1*z) % pp))) % pp; if ( z1 < 0 ) z1 += pp;
            x = x1; z = z1;
        }
        x %= q;
        assert ((x*x-a)%q == 0);
        return x;
    }
}

// Given q a power of an odd prime p and 0<d<4q coprime ot p and congruent to 0 or 3 mod 4, finds a solution to x^2+dy^2=4q with y!=0 if one exists
int i_norm_equation (long *x, long *y, long d, long p, long q)
{
    int dzm4;
    long dd, r, r0, r1, m, s;

    assert ( (p&1) && !(q%p) && d > 0 && d < 4*q && (d%p) && (!(d&3) || (d&3)==3) );
    dzm4 = !(d&3);
    dd = dzm4 ? d/4 : d;
    r = i_sqrt_modq (-dd, p, q);
    if ( r < 0 ) return 0;
    if ( (r&1) ^ (dd&1) ) r = q-r;  // make parity of r match parity of d (so r^2 = -d mod 4)
    m = dzm4 ? q : 4*q;
    if ( r > (m>>1) ) r = m-r;
    r0 = m; r1 = r;
    while ( r1 > m/r1 ) { r = r0 % r1; r0 = r1; r1 = r; }
    if ( r1*r1 >= m )  { r = r0 % r1; r0 = r1; r1 = r; }
    if ( (m-r1*r1)%dd != 0 ) return 0;
    s = i_sqrt ((m-r1*r1)/dd);
    if ( ! s ) return 0;
    *x = dzm4 ? 2*r1 : r1; *y = s;
    assert ( (*x)*(*x) + (*y)*(*y)*d == 4*q );
    return 1;
}