#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

/*
    Micro-benchmark + A/B test harness for the Cornacchia Euclidean descent.

    The descent reduces (a,b) = (2p, x0) by a_i-2 = q_i a_i-1 + a_i until b <= L
    = floor(sqrt(4p)); the final b is the trace candidate t.  This isolates that
    reduction (the ~97% of substep 2) from the rest of the scan so we can measure
    step counts and compare implementations on identical inputs.

    Usage: descent_bench <pbits> <ntrials>
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

static unsigned long g_steps;

// ---- reference: current mpz descent (pointer-rotated to drop the per-step copies) ----
static void descent_mpz (mpz_t outb, const mpz_t a0, const mpz_t b0, const mpz_t L)
{
    static mpz_t s[3]; static int init;
    if (!init) { mpz_init(s[0]); mpz_init(s[1]); mpz_init(s[2]); init=1; }
    mpz_set (s[0], a0); mpz_set (s[1], b0);
    int ai=0, bi=1, ri=2;
    while ( mpz_cmp (s[bi], L) > 0 ) {
        mpz_mod (s[ri], s[ai], s[bi]);
        int t=ai; ai=bi; bi=ri; ri=t;          // rotate: (a,b,r)->(b,r,a)
        g_steps++;
    }
    mpz_set (outb, s[bi]);
}

// ---- mpn descent: fixed limb buffers, mpn_tdiv_qr per step, pointer rotation ----
static void descent_mpn (mpz_t outb, const mpz_t a0, const mpz_t b0, const mpz_t L)
{
    mp_size_t cap = mpz_size (a0) + 1;
    mp_limb_t A[cap], B[cap], R[cap], Q[cap];
    mp_limb_t *a=A,*b=B,*r=R;
    mp_size_t an = mpz_size(a0), bn = mpz_size(b0);
    mpn_copyi (a, a0->_mp_d, an);
    mpn_copyi (b, b0->_mp_d, bn);
    mp_size_t Ln = mpz_size(L);
    const mp_limb_t *Ld = L->_mp_d;

    for (;;) {
        // stop if b <= L
        if ( bn < Ln || (bn == Ln && mpn_cmp (b, Ld, bn) <= 0) ) break;
        // r = a mod b   (a has an limbs >= bn = b's limbs)
        mpn_tdiv_qr (Q, r, 0, a, an, b, bn);
        mp_size_t rn = bn; while ( rn > 0 && r[rn-1] == 0 ) rn--;
        // rotate (a,b,r) -> (b,r,a)
        mp_limb_t *t=a; a=b; b=r; r=t;
        an = bn; bn = rn;
        g_steps++;
        if ( bn == 0 ) break;
    }
    mpz_set_ui (outb, 0);
    if ( bn ) { mpz_import (outb, bn, -1, sizeof(mp_limb_t), 0, 0, b); }
}

// ---- mpn descent with a small-quotient subtractive fast path (q=1,2 via subtraction) ----
static void descent_mpn_fast (mpz_t outb, const mpz_t a0, const mpz_t b0, const mpz_t L)
{
    mp_size_t cap = mpz_size (a0) + 1;
    mp_limb_t A[cap], B[cap], R[cap], Q[cap];
    mp_limb_t *a=A,*b=B,*r=R;
    mp_size_t an = mpz_size(a0), bn = mpz_size(b0);
    mpn_copyi (a, a0->_mp_d, an);
    mpn_copyi (b, b0->_mp_d, bn);
    mp_size_t Ln = mpz_size(L);
    const mp_limb_t *Ld = L->_mp_d;
    #define RLT_B(rn) ((rn) < bn || ((rn) == bn && mpn_cmp (r, b, bn) < 0))   // r < b ?
    for (;;) {
        if ( bn < Ln || (bn == Ln && mpn_cmp (b, Ld, bn) <= 0) ) break;
        mp_size_t rn;
        if ( an == bn ) {                                    // quotient is small
            mpn_sub_n (r, a, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--;     // r = a-b (q=1?)
            if ( ! RLT_B (rn) ) {                            // r >= b -> q >= 2
                mpn_sub_n (r, r, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--; // r -= b (q=2?)
                if ( ! RLT_B (rn) ) {                        // q >= 3 -> full division
                    mpn_tdiv_qr (Q, r, 0, a, an, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--;
                }
            }
        } else {
            mpn_tdiv_qr (Q, r, 0, a, an, b, bn); rn = bn; while ( rn && !r[rn-1] ) rn--;
        }
        mp_limb_t *t=a; a=b; b=r; r=t;  an = bn; bn = rn;  g_steps++;
        if ( bn == 0 ) break;
    }
    mpz_set_ui (outb, 0);
    if ( bn ) mpz_import (outb, bn, -1, sizeof(mp_limb_t), 0, 0, b);
}

int main (int argc, char *argv[])
{
    int bits = argc > 1 ? atoi (argv[1]) : 256;
    long N = argc > 2 ? atol (argv[2]) : 1000000;

    mpz_t p, p2, L, b0, bm, bn, a0;
    mpz_inits (p, p2, L, b0, bm, bn, a0, NULL);
    gmp_randstate_t rs; gmp_randinit_default (rs); gmp_randseed_ui (rs, 1);
    mpz_urandomb (p, rs, bits); mpz_setbit (p, bits-1); mpz_nextprime (p, p);
    mpz_mul_2exp (p2, p, 1);                          // a0 = 2p
    mpz_mul_2exp (L, p, 2); mpz_sqrt (L, L);          // L = floor(sqrt(4p))

    // pre-generate N random starting b0 in [0, p)
    mpz_t *bs = malloc (N * sizeof(mpz_t));
    for (long i=0;i<N;i++){ mpz_init(bs[i]); mpz_urandomm(bs[i], rs, p); }

    // mpz reference
    g_steps = 0; double t0 = wall ();
    for (long i=0;i<N;i++) descent_mpz (bm, p2, bs[i], L);
    double tm = wall () - t0;  unsigned long steps_mpz = g_steps;

    // mpn (naive)
    g_steps = 0; t0 = wall ();
    for (long i=0;i<N;i++) descent_mpn (bn, p2, bs[i], L);
    double tn = wall () - t0;

    // mpn + small-quotient fast path
    mpz_t bf; mpz_init (bf);
    t0 = wall ();
    for (long i=0;i<N;i++) descent_mpn_fast (bf, p2, bs[i], L);
    double tf = wall () - t0;

    // correctness: same final b on a sample
    int bad = 0;
    for (long i=0;i<N && i<200000;i++){
        descent_mpz(bm,p2,bs[i],L); descent_mpn(bn,p2,bs[i],L); descent_mpn_fast(bf,p2,bs[i],L);
        if(mpz_cmp(bm,bn)||mpz_cmp(bm,bf)) bad++;
    }

    printf ("%d-bit p, %ld descents:  avg steps %.1f\n", bits, N, (double)steps_mpz/N);
    printf ("  mpz      : %.3f s  (%.4f us/descent)\n", tm, 1e6*tm/N);
    printf ("  mpn      : %.3f s  (%.4f us/descent)   %.2fx\n", tn, 1e6*tn/N, tm/tn);
    printf ("  mpn+fastq: %.3f s  (%.4f us/descent)   %.2fx\n", tf, 1e6*tf/N, tm/tf);
    printf ("  mismatch : %d\n", bad);
    return 0;
}
