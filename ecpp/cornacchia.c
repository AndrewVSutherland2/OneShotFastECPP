#include <stdlib.h>
#include <string.h>
#include "cornacchia.h"

// r = a*b mod p, using ctx scratch t1
static inline void mulmod (mpz_t r, const mpz_t a, const mpz_t b, const mpz_t p, mpz_t tmp)
{
    mpz_mul (tmp, a, b);
    mpz_mod (r, tmp, p);
}

void cornacchia_init (cornacchia_ctx *c, const mpz_t p)
{
    memset (c, 0, sizeof(*c));
    mpz_init_set (c->p, p);
    mpz_init (c->p4);   mpz_mul_2exp (c->p4, p, 2);            // 4p
    mpz_init (c->L);    mpz_sqrt (c->L, c->p4);               // floor(sqrt(4p))
    c->pmod8 = (int) mpz_fdiv_ui (p, 8);

    mpz_init (c->Q);    mpz_sub_ui (c->Q, p, 1);              // p-1 = Q*2^S
    c->S = 0;
    while ( mpz_even_p (c->Q) ) { mpz_fdiv_q_2exp (c->Q, c->Q, 1); c->S++; }
    mpz_init (c->e);    mpz_add_ui (c->e, c->Q, 1);  mpz_fdiv_q_2exp (c->e, c->e, 1); // (Q+1)/2
    mpz_init (c->z);
    mpz_init (c->pp1q);
    if ( (c->pmod8 & 3) == 3 ) {                              // p == 3 mod 4
        mpz_add_ui (c->pp1q, p, 1);  mpz_fdiv_q_2exp (c->pp1q, c->pp1q, 2);  // (p+1)/4
    } else {                                                  // p == 1 mod 4: need a non-residue
        mpz_t n;  mpz_init_set_ui (n, 2);
        while ( mpz_jacobi (n, p) != -1 ) mpz_add_ui (n, n, 1);
        mpz_powm (c->z, n, c->Q, p);
        mpz_clear (n);
    }
    mpz_init (c->x0); mpz_init (c->a); mpz_init (c->b); mpz_init (c->r); mpz_init (c->c);
    mpz_init (c->t0); mpz_init (c->t1); mpz_init (c->R); mpz_init (c->cc); mpz_init (c->bb);
    mpz_init (c->two_p);  mpz_mul_2exp (c->two_p, p, 1);
    c->dcap = mpz_size (c->two_p) + 2;
    c->da = malloc (c->dcap * sizeof(mp_limb_t));
    c->db = malloc (c->dcap * sizeof(mp_limb_t));
    c->dr = malloc (c->dcap * sizeof(mp_limb_t));
    c->dq = malloc (c->dcap * sizeof(mp_limb_t));
}

void cornacchia_clear (cornacchia_ctx *c)
{
    mpz_clear (c->p); mpz_clear (c->p4); mpz_clear (c->L); mpz_clear (c->Q);
    mpz_clear (c->e); mpz_clear (c->z); mpz_clear (c->pp1q);
    mpz_clear (c->x0); mpz_clear (c->a); mpz_clear (c->b); mpz_clear (c->r); mpz_clear (c->c);
    mpz_clear (c->t0); mpz_clear (c->t1); mpz_clear (c->R); mpz_clear (c->cc); mpz_clear (c->bb);
    mpz_clear (c->two_p);  free (c->da); free (c->db); free (c->dr); free (c->dq);
}

int cornacchia_sqrtmodp (cornacchia_ctx *c, mpz_t out, const mpz_t a)
{
    c->n_sqrt++;
    if ( mpz_sgn (a) == 0 ) { mpz_set_ui (out, 0); return 1; }

    if ( (c->pmod8 & 3) == 3 ) {                              // p == 3 mod 4
        // compute into scratch (out may alias a, and we must keep a for the check)
        mpz_powm (c->t0, a, c->pp1q, c->p);                  // candidate a^((p+1)/4)
        mulmod (c->R, c->t0, c->t0, c->p, c->t1);
        if ( mpz_cmp (c->R, a) != 0 ) return 0;              // QR iff it squares back to a
        mpz_set (out, c->t0);
        return 1;
    }

    // Tonelli-Shanks for p == 1 mod 4
    int M = c->S;
    mpz_set (c->bb, c->z);                                    // running "c" = z
    mpz_powm (c->R, a, c->Q, c->p);                          // t = a^Q
    mpz_powm (out, a, c->e, c->p);                           // R = a^((Q+1)/2)
    for (;;) {
        if ( mpz_cmp_ui (c->R, 1) == 0 ) return 1;           // t == 1 -> out is the root
        int i = 0;
        mpz_set (c->t0, c->R);
        while ( mpz_cmp_ui (c->t0, 1) != 0 ) {
            mulmod (c->t0, c->t0, c->t0, c->p, c->t1);
            if ( ++i == M ) return 0;                        // order does not divide 2^M -> non-residue
        }
        mpz_set (c->r, c->bb);                               // b = c^(2^(M-i-1))
        for ( int j = 0 ; j < M - i - 1 ; j++ ) mulmod (c->r, c->r, c->r, c->p, c->t1);
        M = i;
        mulmod (c->bb, c->r, c->r, c->p, c->t1);             // c = b^2
        mulmod (c->R, c->R, c->bb, c->p, c->t1);             // t = t*c
        mulmod (out, out, c->r, c->p, c->t1);                // R = R*b
    }
}

// Given c->x0 = a square root of -d mod p, do the parity fix, Euclidean descent
// and final square test.  Returns 1 and sets t,v on success.
static int cornacchia_finish (cornacchia_ctx *c, unsigned long d, mpz_t t, mpz_t v)
{
    if ( mpz_tstbit (c->x0, 0) ^ (int)(d & 1) ) mpz_sub (c->x0, c->p, c->x0);  // parity x0 == d (mod 2)

    c->n_descent++;
    // mpn Euclidean descent: (a,b) = (2p, x0) reduce until b <= L.  Most quotients
    // are 1 or 2, handled by subtraction; only larger ones need a full division.
    mp_limb_t *a = c->da, *b = c->db, *r = c->dr;
    mp_size_t an = c->two_p->_mp_size, bn = c->x0->_mp_size;
    mpn_copyi (a, c->two_p->_mp_d, an);
    mpn_copyi (b, c->x0->_mp_d, bn);
    mp_size_t Ln = c->L->_mp_size;
    const mp_limb_t *Ld = c->L->_mp_d;
    for (;;) {
        if ( bn < Ln || (bn == Ln && mpn_cmp (b, Ld, bn) <= 0) ) break;
        mp_size_t rn;
        if ( an == bn ) {
            mpn_sub_n (r, a, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--;        // q=1?
            if ( ! (rn < bn || (rn == bn && mpn_cmp (r, b, bn) < 0)) ) {
                mpn_sub_n (r, r, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--;    // q=2?
                if ( ! (rn < bn || (rn == bn && mpn_cmp (r, b, bn) < 0)) )          // q>=3: divide
                    { mpn_tdiv_qr (c->dq, r, 0, a, an, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--; }
            }
        } else {
            mpn_tdiv_qr (c->dq, r, 0, a, an, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--;
        }
        mp_limb_t *tmp = a; a = b; b = r; r = tmp;  an = bn; bn = rn;
        if ( bn == 0 ) break;
    }
    if ( bn == 0 ) mpz_set_ui (c->b, 0);
    else mpz_import (c->b, bn, -1, sizeof(mp_limb_t), 0, 0, b);

    mpz_mul (c->cc, c->b, c->b);  mpz_sub (c->r, c->p4, c->cc);   // 4p - b^2
    if ( mpz_sgn (c->r) < 0 ) return 0;
    if ( ! mpz_divisible_ui_p (c->r, d) ) return 0;
    mpz_divexact_ui (c->cc, c->r, d);
    if ( ! mpz_perfect_square_p (c->cc) ) return 0;
    mpz_sqrt (v, c->cc);
    mpz_set (t, c->b);
    return 1;
}

int cornacchia_solve (cornacchia_ctx *c, unsigned long d, mpz_t t, mpz_t v)
{
    int dm4 = (int)(d & 3);
    if ( d == 0 || (dm4 != 0 && dm4 != 3) ) return 0;
    if ( mpz_cmp_ui (c->p4, d) <= 0 ) return 0;              // need d < 4p

    // x0 = -d mod p
    mpz_set_ui (c->t0, d);  mpz_mod (c->t0, c->t0, c->p);
    if ( mpz_sgn (c->t0) == 0 ) return 0;                    // p | d
    mpz_sub (c->x0, c->p, c->t0);

    // Legendre pre-filter: (-d/p) must be 1 (skips a doomed Tonelli for ~half of D)
    if ( mpz_jacobi (c->x0, c->p) != 1 ) return 0;
    if ( ! cornacchia_sqrtmodp (c, c->x0, c->x0) ) return 0; // x0 = sqrt(-d) mod p
    return cornacchia_finish (c, d, t, v);
}

// As cornacchia_solve but with x0 = sqrt(-d) mod p supplied (e.g. computed
// multiplicatively from a factor base).  Caller guarantees (-d/p) = 1.
int cornacchia_solve_x0 (cornacchia_ctx *c, unsigned long d, const mpz_t x0, mpz_t t, mpz_t v)
{
    mpz_set (c->x0, x0);
    return cornacchia_finish (c, d, t, v);
}
