#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>
#include "cornacchia.h"
#include "fproot.h"
#include "cminv.h"
#include "invj.h"

/*
    cm_method D p : given a discriminant D < -4 and a prime p for which the norm
    equation 4p = t^2 - v^2 D = t^2 + |D| v^2 is solvable, output the j-invariant
    of an elliptic curve E/F_p with CM by the order of discriminant D and Frobenius
    trace +/-t.

      cm_method D=<disc<-4> (p=<dec> | pbits=<n> [seed=<s>]) [ells=<l1,l2,...>] [v]
                [jobs=<n>] [tower=0|1] [tn=<n>]

    It picks the best class invariant internally (classpoly inv=-1), computes the
    small H_D^inv mod p, finds a root f0 with fproot, converts f0 to a j-invariant
    with the big-F_p invariant->j map (cm_j_from_inv), and verifies the resulting
    curve has trace +/-t.  Needs the classpoly environment (source setenv.sh).

    By default the class polynomial is not computed as such: classpoly is asked
    (CLASSPOLY_TOWER=1) for its decomposition through a subgroup of the class
    group, as in Algorithm 1 of arXiv:1009.1082 -- a polynomial V of degree m (the
    number of cosets) and interpolants W_0..W_{n-2} with m coefficients each, n the
    order of the subgroup.  A root y of V gives the factor

        U(X) = X^n + y X^{n-1} + sum_{k<=n-2} (W_k(y)/V'(y)) X^k

    of H_D^inv mod p, so f0 costs two root-findings of degrees m and n instead of
    one of degree mn, and the CRT needs fewer primes.  tower=0 turns this off;
    tn=<n> forces the subgroup order.  The height bound behind the decomposition is
    heuristic for invariants other than j and V can have a repeated root, so every
    failure on that path (V or U not splitting, V'(y) = 0, no candidate j with the
    right trace) falls back to the plain class polynomial.

    ells=<l1,l2,...> (primes <= 97, each dividing v): after verifying j, descend
    the ell-volcano to its floor for each listed ell, and output the floor
    j-invariant instead.  The H_D roots have non-cyclic ell-torsion exactly for
    the primes ell | n1 = gcd(a-1,b) (pi = a+b*omega); a floor curve has cyclic
    ell-Sylow, which (a) makes it Montgomery-representable even when p = 3 mod 4
    and N = 4 mod 8 (the unique 2-torsion point sits under a point of order 4, so
    it is halvable and f'(e1) is a square), and (b) lets a certificate m use the
    full ell-power of N (m | exponent(E) = N/n1', and n1' loses its ell part).
    Isogenous curves have the same trace, so the certificate is unaffected.
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

// ---- affine Weierstrass point arithmetic over F_p:  y^2 = x^3 + A4 x + A6 ----
typedef struct { mpz_t x, y; int inf; } wpt;
static void winit (wpt *P) { mpz_init (P->x); mpz_init (P->y); P->inf = 1; }
static void wclear (wpt *P) { mpz_clear (P->x); mpz_clear (P->y); }
static void wset (wpt *R, const wpt *P) { mpz_set (R->x, P->x); mpz_set (R->y, P->y); R->inf = P->inf; }

// R = P + Q  (R may alias neither P nor Q)
static void wadd (wpt *R, const wpt *P, const wpt *Q, const mpz_t A4, const mpz_t p)
{
    if ( P->inf ) { wset (R, Q); return; }
    if ( Q->inf ) { wset (R, P); return; }
    mpz_t lam, t, u;  mpz_init (lam); mpz_init (t); mpz_init (u);
    if ( mpz_cmp (P->x, Q->x) == 0 ) {
        mpz_add (t, P->y, Q->y);  mpz_mod (t, t, p);
        if ( mpz_cmp (P->y, Q->y) != 0 || mpz_sgn (P->y) == 0 || mpz_sgn (t) == 0 ) { R->inf = 1; goto done; }
        // doubling: lam = (3x^2 + A4) / (2y)
        mpz_mul (t, P->x, P->x);  mpz_mul_ui (t, t, 3);  mpz_add (t, t, A4);  mpz_mod (t, t, p);
        mpz_mul_ui (u, P->y, 2);  mpz_invert (u, u, p);
        mpz_mul (lam, t, u);  mpz_mod (lam, lam, p);
    } else {
        mpz_sub (t, Q->y, P->y);  mpz_sub (u, Q->x, P->x);  mpz_invert (u, u, p);
        mpz_mul (lam, t, u);  mpz_mod (lam, lam, p);
    }
    mpz_mul (t, lam, lam);  mpz_sub (t, t, P->x);  mpz_sub (t, t, Q->x);  mpz_mod (t, t, p);   // x3
    mpz_sub (u, P->x, t);  mpz_mul (u, u, lam);  mpz_sub (u, u, P->y);  mpz_mod (u, u, p);      // y3
    mpz_set (R->x, t);  mpz_set (R->y, u);  R->inf = 0;
done:
    mpz_clear (lam); mpz_clear (t); mpz_clear (u);
}

// R = [k] P
static void wmul (wpt *R, const mpz_t k, const wpt *P, const mpz_t A4, const mpz_t p)
{
    wpt acc, base, tmp;  winit (&acc); winit (&base); winit (&tmp);
    wset (&base, P);
    for ( int i = 0 ; i < (int) mpz_sizeinbase (k, 2) ; i++ ) {
        if ( mpz_tstbit (k, i) ) { wadd (&tmp, &acc, &base, A4, p);  wset (&acc, &tmp); }
        wadd (&tmp, &base, &base, A4, p);  wset (&base, &tmp);
    }
    wset (R, &acc);
    wclear (&acc); wclear (&base); wclear (&tmp);
}

// Does E: y^2=x^3+A4 x+A6 over F_p have trace +/- t?  Returns +1 if #E=p+1-t,
// -1 if #E=p+1+t, 0 if neither (j is not a CM-by-D j-invariant).
static int trace_sign (const mpz_t A4, const mpz_t A6, const mpz_t p, const mpz_t t,
                       cornacchia_ctx *cc)
{
    mpz_t rhs, x, y, np1mt, np1pt;  mpz_inits (rhs, x, y, np1mt, np1pt, NULL);
    mpz_add_ui (np1mt, p, 1);  mpz_sub (np1mt, np1mt, t);      // p+1-t
    mpz_add_ui (np1pt, p, 1);  mpz_add (np1pt, np1pt, t);      // p+1+t
    int sign = 0;
    gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, 271828);
    for ( int tries = 0 ; tries < 40 && ! sign ; tries++ ) {
        mpz_urandomm (x, rs, p);
        mpz_mul (rhs, x, x);  mpz_add (rhs, rhs, A4);  mpz_mul (rhs, rhs, x);  mpz_add (rhs, rhs, A6);  mpz_mod (rhs, rhs, p);
        if ( ! cornacchia_sqrtmodp (cc, y, rhs) ) continue;    // need rhs a QR to get a point
        wpt P, Q;  winit (&P); winit (&Q);
        mpz_set (P.x, x);  mpz_set (P.y, y);  P.inf = 0;
        wmul (&Q, np1mt, &P, A4, p);
        if ( Q.inf ) sign = 1;
        else { wmul (&Q, np1pt, &P, A4, p);  if ( Q.inf ) sign = -1; }
        wclear (&P); wclear (&Q);
    }
    gmp_randclear (rs);
    mpz_clears (rhs, x, y, np1mt, np1pt, NULL);
    return sign;
}

// ---- walk-to-the-floor on the ell-isogeny volcano ----
// jm (Montgomery) is a vertex with non-cyclic ell-torsion, so Phi_ell(X, jm)
// splits completely (ell+1 rational neighbors: <= 2 along the crater, the rest
// down).  Start three non-backtracking walkers (next = root of Phi_ell(X, j_cur)
// other than j_prev) in lockstep: at least one heads down, and below the surface
// every onward edge descends, so it reaches the floor -- detected by Phi_ell
// having no rational root besides j_prev -- within the volcano depth v_ell(v).
// Crater walkers never terminate and are cut off by maxsteps.  Replaces jm by a
// floor j-invariant and returns 1; returns 0 only on broken input assumptions.
static int volcano_floor (const fp_ctx *C, const bipoly *P, mp_limb_t *jm, int maxsteps, uint64_t seed)
{
    int s = C->s, maxr = P->maxj + 2;
    mp_limb_t *r = malloc ((size_t) maxr * s * sizeof(mp_limb_t));
    mp_limb_t cur[3][64], prev[3][64], sj[2][64];
    fp_set_ui (C, sj[0], 0);  fp_set_ui (C, sj[1], 1728);      // singular j's the curve
                                                               // formulas can't take
    int nw = invj_jroots (C, P, jm, r, maxr, seed);
    if ( nw < 1 ) { free (r); return 0; }                      // not rank-2 ell-torsion
    if ( nw > 3 ) nw = 3;
    for ( int w = 0 ; w < nw ; w++ ) { mpn_copyi (prev[w], jm, s);  mpn_copyi (cur[w], r + (size_t)w*s, s); }
    for ( int step = 0 ; step <= maxsteps ; step++ ) {
        for ( int w = 0 ; w < nw ; w++ ) {
            int n = invj_jroots (C, P, cur[w], r, maxr, seed + 999*step + w), keep = -1;
            for ( int i = 0 ; i < n && keep < 0 ; i++ )
                if ( ! fp_equal (C, r + (size_t)i*s, prev[w]) ) keep = i;
            if ( keep < 0 ) {                                  // no onward edge: the floor
                if ( fp_equal (C, cur[w], sj[0]) || fp_equal (C, cur[w], sj[1]) ) continue;
                mpn_copyi (jm, cur[w], s);  free (r);  return 1;
            }
            mpn_copyi (prev[w], cur[w], s);  mpn_copyi (cur[w], r + (size_t)keep*s, s);
        }
    }
    free (r);  return 0;
}

// Runs classpoly (inv=-1: best invariant) for D mod p into outf, with `jobs` ECRT
// workers plus a merge when jobs > 1.  tower: request the decomposition
// (CLASSPOLY_TOWER=2: sign-free where that saves computing the trace); tn > 0
// forces its subgroup order; slack > 0 adds that many percent to the height bound
// of the decomposition (CLASSPOLY_TOWER_SLACK).  Returns 0 on success.
static int run_classpoly (long D, const char *pdec, const char *outf, int jobs, int tower, long tn, int slack)
{
    static char cmd[24576];
    char tenv[192], senv[64] = "";
    if ( slack > 0 ) snprintf (senv, sizeof senv, "export CLASSPOLY_TOWER_SLACK=%d; ", slack);
    if ( tn > 0 ) snprintf (tenv, sizeof tenv, "%sexport CLASSPOLY_TOWER=2 CLASSPOLY_TOWER_N=%ld; ", senv, tn);
    else snprintf (tenv, sizeof tenv, "%sexport CLASSPOLY_TOWER=%d; unset CLASSPOLY_TOWER_N; ", senv, tower ? 2 : 0);
    remove (outf);
    if ( jobs > 1 ) {
        // parallel H_D: run `jobs` classpoly ECRT workers (each covers prime indices
        // = jobid-1 mod jobs and dumps partial sums), then a merge pass writes outf.
        char edir[] = "/tmp/cmecrtXXXXXX";
        if ( ! mkdtemp (edir) ) { fprintf (stderr, "mkdtemp failed\n"); return 1; }
        snprintf (cmd, sizeof cmd,
            "%sexport CLASSPOLY_JOBS=%d CLASSPOLY_ECRT_DIR=%s; "
            "for j in $(seq 1 %d); do CLASSPOLY_JOBID=$j classpoly %ld -1 %s %s.w$j -1 & done; wait; "
            "CLASSPOLY_JOBID=0 classpoly %ld -1 %s %s -1",
            tenv, jobs, edir, jobs, D, pdec, outf, D, pdec, outf);
        int rc = system (cmd);
        snprintf (cmd, sizeof cmd, "rm -rf %s", edir);
        if ( system (cmd) != 0 ) fprintf (stderr, "warning: cleanup of %s failed\n", edir);
        if ( rc != 0 ) { fprintf (stderr, "parallel classpoly failed\n"); return 1; }
    } else {
        snprintf (cmd, sizeof cmd, "%sclasspoly %ld -1 %s %s -1", tenv, D, pdec, outf);
        if ( system (cmd) != 0 ) { fprintf (stderr, "classpoly failed\n"); return 1; }
    }
    return 0;
}

// Next coefficient line of a classpoly file ("<integer>*<monomial>"): reads the
// integer into c.  Returns 0 at end of file or on a malformed line.
static int next_coeff (FILE *f, char *ln, int lnsz, mpz_t c)
{
    while ( fgets (ln, lnsz, f) ) {
        char *x = strchr (ln, '*');
        if ( ! x ) continue;
        *x = 0;
        return mpz_set_str (c, ln, 10) == 0;
    }
    return 0;
}

#define ROOT_MAXFAIL 30     // bits of evidence (see fp_find_root_bounded) before a polynomial is declared not to split
#define RETRY_SLACK  20     // percent added to the height bound of a decomposition when the first attempt fails (default: 5)

/*
    Root f0 of the class polynomial in the file outf, over F_p (context C).

    The file is either a plain class polynomial ("c*X^i" lines) or a decomposition
    (header lines F=m, G=n; then V_0..V_{m-1}, 1, then W_0, ..., W_{n-2} with m
    coefficients each, lowest first).  For a decomposition the file is read in one
    pass: V, then a root y of V and the powers of y, then each W_k is evaluated at
    y as its coefficients go by, so nothing of size mn is ever stored.

    A sign-free decomposition (header line S=<T^2>, see classpoly_crt.c) stores
    the coefficients that change sign under f -> -f multiplied by the trace T; they
    are divided by a square root of S here, and the root returned is f0 or -f0,
    which have the same j-invariants for the invariants concerned.

    Returns 1 on success.  *pm, *pn: the degrees used (m = 0 for a plain
    polynomial, whose degree is *pn).  tr[0], tr[1]: seconds in the two
    root-findings.  A return of 0 with *pm > 0 means the decomposition failed (the
    caller falls back); with *pm = 0 it is a parse error.

    fault (testing only) damages a decomposition as a height bound that is too
    small would: 1 adds one to a coefficient of V, 2 to a coefficient of W_0.
    (On the command line fault=1,2 damage the first attempt only, fault=3,4 do
    the same to the retry as well, so that the plain class polynomial is used.)
*/
static int root_from_file (const fp_ctx *C, const char *outf, mpz_t f0, int *pinv, int *pm, int *pn, double tr[2], uint64_t seed, int fault)
{
    static char ln[8192];
    FILE *f = fopen (outf, "r");
    if ( ! f ) { fprintf (stderr, "no %s\n", outf); return 0; }
    int inv = -2, m = 0, n = 0, s = C->s, ok = 0, signfree = 0;
    long pos = 0;
    mpz_t Sz;  mpz_init (Sz);
    for (;;) {                                              // header lines start with a letter
        pos = ftell (f);
        if ( ! fgets (ln, sizeof ln, f) ) break;
        if ( ln[0] == 'I' && ln[1] == '=' ) inv = atoi (ln+2);
        else if ( ln[0] == 'F' && ln[1] == '=' ) m = atoi (ln+2);
        else if ( ln[0] == 'G' && ln[1] == '=' ) n = atoi (ln+2);
        else if ( ln[0] == 'S' && ln[1] == '=' ) { signfree = 1;  mpz_set_str (Sz, ln+2, 10); }
        else if ( ! ( (ln[0] >= 'A' && ln[0] <= 'Z') || (ln[0] >= 'a' && ln[0] <= 'z') ) ) break;
    }
    fseek (f, pos, SEEK_SET);
    *pinv = inv;  *pm = 0;  *pn = 0;  tr[0] = tr[1] = 0.0;
    if ( inv == -2 ) { fclose (f);  mpz_clear (Sz);  fprintf (stderr, "could not parse H_D^inv\n");  return 0; }
    mpz_t c;  mpz_init (c);
    mp_limb_t rm[64];
    fp_poly h;

    if ( m < 2 || n < 2 ) {
        // plain class polynomial: coefficients in order of increasing degree, monic
        int cap = 1024, deg = -1;
        fpoly_init (C, &h, cap);
        while ( next_coeff (f, ln, sizeof ln, c) ) {
            if ( ++deg >= cap ) { cap *= 2;  h.c = realloc (h.c, (size_t)cap * s * sizeof(mp_limb_t));  h.cap = cap; }
            mpz_mod (c, c, C->pz);  fp_set_mpz (C, h.c + (size_t)deg*s, c);
        }
        h.deg = deg;  *pn = deg;
        if ( deg < 1 ) { fprintf (stderr, "could not parse H_D^inv\n"); goto done; }
        if ( h.cap < 2*deg + 2 ) { h.c = realloc (h.c, (size_t)(2*deg+2) * s * sizeof(mp_limb_t));  h.cap = 2*deg+2; }
        double t0 = wall ();
        fp_find_root (C, &h, rm, seed);
        tr[0] = wall () - t0;
        fp_get_mpz (C, f0, rm);
        ok = 1;
        goto done;
    }

    // decomposition.  V, monic of degree m
    *pm = m;  *pn = n;
    fpoly_init (C, &h, 2*(m > n ? m : n) + 4);
    mp_limb_t *V = malloc ((size_t)(m+1) * s * sizeof(mp_limb_t));
    mp_limb_t *yp = malloc ((size_t)m * s * sizeof(mp_limb_t));      // powers of y
    mp_limb_t y[64], d[64], t[64], acc[64], acc1[64], tinv[64];
    if ( signfree ) {
        // tinv = 1/tau with tau^2 = S: a root of X^2 - S
        mpz_mod (Sz, Sz, C->pz);
        if ( ! mpz_sgn (Sz) ) { fprintf (stderr, "sign-free decomposition with S = 0 mod p\n"); goto done2; }
        mpz_sub (c, C->pz, Sz);  fp_set_mpz (C, h.c, c);
        memset (h.c + s, 0, (size_t)s * sizeof(mp_limb_t));  mpn_copyi (h.c + 2*(size_t)s, C->R1, s);  h.deg = 2;
        if ( ! fp_find_root_bounded (C, &h, tinv, seed + 2, ROOT_MAXFAIL) ) { fprintf (stderr, "S is not a square mod p\n"); goto done2; }
        fp_inv (C, tinv, tinv);
    }
    for ( int i = 0 ; i <= m ; i++ ) {
        if ( ! next_coeff (f, ln, sizeof ln, c) ) { fprintf (stderr, "decomposition file truncated in V\n"); goto done2; }
        if ( fault == 1 && i == 0 ) mpz_add_ui (c, c, 1);
        mpz_mod (c, c, C->pz);  fp_set_mpz (C, V + (size_t)i*s, c);
        if ( signfree && ((m-i)&1) ) fp_mul (C, V + (size_t)i*s, V + (size_t)i*s, tinv);       // odd coefficient of V
    }
    if ( ! fp_equal (C, V + (size_t)m*s, C->R1) ) { fprintf (stderr, "V is not monic\n"); goto done2; }
    memcpy (h.c, V, (size_t)(m+1) * s * sizeof(mp_limb_t));  h.deg = m;
    double t0 = wall ();
    int found = fp_find_root_bounded (C, &h, y, seed, ROOT_MAXFAIL);
    tr[0] = wall () - t0;
    if ( ! found ) { fprintf (stderr, "V (degree %d) does not split over F_p\n", m); goto done2; }
    // d = V'(y) by Horner, and the powers 1, y, ..., y^{m-1}
    fp_set_ui (C, d, (unsigned long) m);                                // leading coefficient of V' (V monic)
    for ( int i = m-1 ; i >= 1 ; i-- ) {
        fp_mul (C, d, d, y);
        fp_set_ui (C, t, (unsigned long) i);  fp_mul (C, t, t, V + (size_t)i*s);
        fp_add (C, d, d, t);
    }
    if ( fp_is_zero (C, d) ) { fprintf (stderr, "the root of V found is a root of V'\n"); goto done2; }
    fp_inv (C, d, d);
    mpn_copyi (yp, C->R1, s);
    for ( int i = 1 ; i < m ; i++ ) fp_mul (C, yp + (size_t)i*s, yp + (size_t)(i-1)*s, y);
    // U = X^n + y X^{n-1} + sum_k (W_k(y)/V'(y)) X^k
    for ( int k = 0 ; k < n-1 ; k++ ) {
        memset (acc, 0, (size_t)s * sizeof(mp_limb_t));  memset (acc1, 0, (size_t)s * sizeof(mp_limb_t));
        for ( int i = 0 ; i < m ; i++ ) {
            if ( ! next_coeff (f, ln, sizeof ln, c) ) { fprintf (stderr, "decomposition file truncated in W_%d\n", k); goto done2; }
            if ( fault == 2 && k == 0 && i == 0 ) mpz_add_ui (c, c, 1);
            mpz_mod (c, c, C->pz);  fp_set_mpz (C, t, c);
            fp_mul (C, t, t, yp + (size_t)i*s);
            if ( signfree && ((n-k+m-1-i)&1) ) fp_add (C, acc1, acc1, t);          // odd coefficients of W_k, still times the trace
            else fp_add (C, acc, acc, t);
        }
        if ( signfree ) { fp_mul (C, acc1, acc1, tinv);  fp_add (C, acc, acc, acc1); }
        fp_mul (C, h.c + (size_t)k*s, acc, d);
    }
    mpn_copyi (h.c + (size_t)(n-1)*s, y, s);
    mpn_copyi (h.c + (size_t)n*s, C->R1, s);
    h.deg = n;
    t0 = wall ();
    found = fp_find_root_bounded (C, &h, rm, seed + 1, ROOT_MAXFAIL);
    tr[1] = wall () - t0;
    if ( ! found ) { fprintf (stderr, "U (degree %d) does not split over F_p\n", n); goto done2; }
    fp_get_mpz (C, f0, rm);
    ok = 1;
done2:
    free (V);  free (yp);
done:
    fpoly_clear (&h);  mpz_clear (c);  mpz_clear (Sz);  fclose (f);
    return ok;
}

int main (int argc, char **argv)
{
    long D = 0;  unsigned long pbits = 0, seed = 1;  int want_v = 0, jobs = 1;
    int tower = 1, fault = 0;  long tn = 0;
    int ells[16], nells = 0;
    mpz_t p;  mpz_init (p);
    for ( int i = 1 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "D=", 2) ) D = strtol (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i]+2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i]+6, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i]+5, 0, 10);
        else if ( ! strncmp (argv[i], "jobs=", 5) ) jobs = atoi (argv[i]+5);
        else if ( ! strncmp (argv[i], "tower=", 6) ) tower = atoi (argv[i]+6);
        else if ( ! strncmp (argv[i], "tn=", 3) ) tn = atol (argv[i]+3);
        else if ( ! strncmp (argv[i], "fault=", 6) ) fault = atoi (argv[i]+6);     // testing only, see root_from_file
        else if ( ! strncmp (argv[i], "ells=", 5) ) {
            for ( char *q = argv[i]+5 ; *q && nells < 16 ; ) {
                int l = (int) strtol (q, &q, 10);
                if ( l >= 2 && l <= 97 ) ells[nells++] = l;
                if ( *q == ',' ) q++;
            }
        }
        else if ( ! strcmp (argv[i], "v") ) want_v = 1;
    }
    if ( D >= -4 ) { fprintf (stderr, "need D < -4\n"); return 2; }
    if ( ! mpz_sgn (p) ) {
        if ( ! pbits ) pbits = 256;
        gmp_randstate_t rs; gmp_randinit_default (rs); gmp_randseed_ui (rs, seed);
        mpz_urandomb (p, rs, pbits); mpz_setbit (p, pbits-1); mpz_nextprime (p, p); gmp_randclear (rs);
    }
    if ( ! mpz_probab_prime_p (p, 25) ) { fprintf (stderr, "p not prime\n"); return 1; }
    const char *phidir = getenv ("CLASSPOLY_PHI_DIR");
    if ( ! phidir ) { fprintf (stderr, "set CLASSPOLY_PHI_DIR (source setenv.sh)\n"); return 1; }

    // 1. Cornacchia: 4p = t^2 + |D| v^2
    cornacchia_ctx cc;  cornacchia_init (&cc, p);
    unsigned long d = (unsigned long)(-D);
    mpz_t t, v;  mpz_init (t); mpz_init (v);
    if ( ! cornacchia_solve (&cc, d, t, v) ) { fprintf (stderr, "norm equation 4p=t^2+|D|v^2 not solvable for D=%ld\n", D); return 1; }

    // 2-4. class polynomial (or its decomposition) mod p with the best invariant,
    // a root f0 of it, the candidate j-invariants for f0, and the one with CM by D
    // (trace +/- t).  Pass 0 uses the decomposition; if anything on that path fails,
    // pass 1 repeats it with RETRY_SLACK percent more CRT primes (the usual cause is
    // the heuristic height bound), and pass 2 uses the plain class polynomial.
    char *pdec = mpz_get_str (NULL, 10, p);
    char outf[256];   snprintf (outf, sizeof outf, "/tmp/cm_%ld_%d.txt", -D, (int) getpid ());
    // Parallel classpoly only pays off when H_D is large: each worker repeats the
    // setup (class group, phi loads, ECRT precomputation), so for small h(D) the
    // single-job path is faster.  h(D) ~ sqrt(|D|), so gate on |D|.
    if ( jobs > 1 && d < 20000000UL ) jobs = 1;
    fp_ctx C;  fp_init (&C, p);
    mpz_t f0;  mpz_init (f0);
    mpz_t J[8];  for ( int i = 0 ; i < 8 ; i++ ) mpz_init (J[i]);
    mpz_t A4, A6, k, tmp;  mpz_inits (A4, A6, k, tmp, NULL);
    int found = 0, tsign = 0, inv = -2, nj = 0, dm = 0, dn = 0, fellback = 0;
    mpz_t jout;  mpz_init (jout);
    double t_hd = 0.0, tr[2] = { 0.0, 0.0 };
    if ( tn > 0 ) tower = 1;
    for ( int pass = tower ? 0 : 2 ; pass < 3 && ! found ; pass++ ) {
        double t0 = wall ();
        if ( run_classpoly (D, pdec, outf, jobs, pass < 2, pass < 2 ? tn : 0, pass == 1 ? RETRY_SLACK : 0) ) return 1;
        t_hd += wall () - t0;
        int ok = root_from_file (&C, outf, f0, &inv, &dm, &dn, tr, 0xC0FFEE + (unsigned) seed + 16*(unsigned) pass,
                                 pass == 0 ? (fault > 2 ? fault - 2 : fault) : (pass == 1 && fault > 2 ? fault - 2 : 0));
        remove (outf);
        if ( ! ok && ! dm ) return 1;                          // not a decomposition: nothing to fall back to
        if ( ok ) {
            // classpoly's invariant->j map aborts when Phi_f(f0,J) has no root, which is what a value f0 from a
            // wrong decomposition gives; on that path ask the reference implementation first, which just reports it.
            nj = ( dm ? cm_j_from_inv_ref (J, 8, f0, p, inv, phidir) : 1 );
            if ( nj > 0 ) nj = cm_j_from_inv (J, 8, f0, p, inv, phidir);
            if ( nj <= 0 && ! dm ) { fprintf (stderr, "cm_j_from_inv failed (inv=%d)\n", inv); return 1; }
            for ( int i = 0 ; i < nj && ! found ; i++ ) {
                // curve with j-invariant J[i]:  y^2 = x^3 + 3k x + 2k,  k = j/(1728-j)
                mpz_ui_sub (tmp, 1728, J[i]);  mpz_mod (tmp, tmp, p);
                if ( mpz_sgn (tmp) == 0 ) continue;                    // j=1728 (D=-4), excluded
                mpz_invert (tmp, tmp, p);  mpz_mul (k, J[i], tmp);  mpz_mod (k, k, p);
                mpz_mul_ui (A4, k, 3);  mpz_mod (A4, A4, p);
                mpz_mul_ui (A6, k, 2);  mpz_mod (A6, A6, p);
                int sg = trace_sign (A4, A6, p, t, &cc);
                if ( sg ) { found = 1;  tsign = sg;  mpz_set (jout, J[i]); }
            }
        }
        if ( found ) break;
        if ( ! dm ) { fprintf (stderr, "no CM-by-D j among %d candidates (inv=%d)\n", nj, inv); return 3; }
        fprintf (stderr, "decomposition (m=%d, n=%d, inv=%d) gave no curve with CM by D=%ld; %s\n", dm, dn, inv, D,
                 pass == 0 ? "retrying it with a larger height bound" : "falling back to the plain class polynomial");
        fellback++;
    }
    if ( ! found ) { fprintf (stderr, "no CM-by-D j (inv=%d)\n", inv); return 3; }
    free (pdec);

    // 5. optional ell-volcano descents: replace jout by a floor j-invariant
    if ( nells ) {
        const char *phid = phidir;
        mp_limb_t jm[64];
        fp_set_mpz (&C, jm, jout);
        for ( int e = 0 ; e < nells ; e++ ) {
            int ell = ells[e], depth = 0;
            mpz_set (tmp, v);
            while ( mpz_divisible_ui_p (tmp, ell) ) { mpz_divexact_ui (tmp, tmp, ell); depth++; }
            if ( ! depth ) continue;                           // ell does not divide v: already cyclic
            bipoly Phi;
            if ( ! invj_load_phi (&Phi, phid, ell) ) { fprintf (stderr, "no phi_j_%d.txt in %s\n", ell, phid); return 4; }
            int ok = volcano_floor (&C, &Phi, jm, depth + 2, 0xF10E * seed + ell);
            invj_clear (&Phi);
            if ( ! ok ) { fprintf (stderr, "ell=%d volcano descent failed (D=%ld)\n", ell, D); return 4; }
        }
        fp_get_mpz (&C, jout, jm);
        // an isogenous curve has trace +/-t; re-verify as a sanity check.  The
        // sign can flip: the k-formula model below is an arbitrary one of the two
        // quadratic twists with invariant jout, so report the sign it lands on.
        mpz_ui_sub (tmp, 1728, jout);  mpz_mod (tmp, tmp, p);
        mpz_invert (tmp, tmp, p);  mpz_mul (k, jout, tmp);  mpz_mod (k, k, p);
        mpz_mul_ui (A4, k, 3);  mpz_mod (A4, A4, p);
        mpz_mul_ui (A6, k, 2);  mpz_mod (A6, A6, p);
        tsign = trace_sign (A4, A6, p, t, &cc);
        if ( ! tsign ) { fprintf (stderr, "floor curve has wrong trace (D=%ld)\n", D); return 4; }
    }

    // 6. output
    gmp_printf ("D %ld\n", D);
    gmp_printf ("p %Zd\n", p);
    gmp_printf ("inv %d\n", inv);
    gmp_printf ("t %s%Zd\n", tsign > 0 ? "" : "-", t);         // Frobenius trace (#E = p+1-trace)
    if ( want_v ) gmp_printf ("v %Zd\n", v);
    gmp_printf ("j %Zd\n", jout);
    if ( dm ) fprintf (stderr, "[H_D^inv(deg %d = %d x %d, inv %d): %.3fs   roots: %.3fs + %.3fs   %d j-candidate(s)%s", dm*dn, dm, dn, inv, t_hd, tr[0], tr[1], nj,
                       fellback ? "  (second attempt, larger height bound)" : "");
    else fprintf (stderr, "[H_D^inv(deg %d, inv %d): %.3fs   root: %.3fs   %d j-candidate(s)%s", dn, inv, t_hd, tr[0], nj, fellback ? "  (after a failed decomposition)" : "");
    for ( int e = 0 ; e < nells ; e++ ) fprintf (stderr, "%s ell=%d floor", e ? "," : "  ", ells[e]);
    fprintf (stderr, "]\n");

    for ( int i = 0 ; i < 8 ; i++ ) mpz_clear (J[i]);
    mpz_clears (A4, A6, k, tmp, jout, t, v, f0, NULL);
    fp_clear (&C);  cornacchia_clear (&cc);  mpz_clear (p);
    return 0;
}
