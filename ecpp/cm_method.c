#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "cornacchia.h"
#include "fproot.h"
#include "cminv.h"

/*
    cm_method D p : given a discriminant D < -4 and a prime p for which the norm
    equation 4p = t^2 - v^2 D = t^2 + |D| v^2 is solvable, output the j-invariant
    of an elliptic curve E/F_p with CM by the order of discriminant D and Frobenius
    trace +/-t.

      cm_method D=<disc<-4> (p=<dec> | pbits=<n> [seed=<s>]) [v]

    It picks the best class invariant internally (classpoly inv=-1), computes the
    small H_D^inv mod p, finds a root f0 with fproot, converts f0 to a j-invariant
    with the big-F_p invariant->j map (cm_j_from_inv), and verifies the resulting
    curve has trace +/-t.  Needs the classpoly environment (source setenv.sh).
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

int main (int argc, char **argv)
{
    long D = 0;  unsigned long pbits = 0, seed = 1;  int want_v = 0, jobs = 1;
    mpz_t p;  mpz_init (p);
    for ( int i = 1 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "D=", 2) ) D = strtol (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i]+2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i]+6, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i]+5, 0, 10);
        else if ( ! strncmp (argv[i], "jobs=", 5) ) jobs = atoi (argv[i]+5);
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

    // 2. best-invariant H_D^inv mod p (classpoly inv=-1); parse I=<code> and coeffs
    char pdec[8192];  gmp_snprintf (pdec, sizeof pdec, "%Zd", p);
    char outf[256];   snprintf (outf, sizeof outf, "/tmp/cm_%ld.txt", -D);
    char cmd[16384];
    double t0 = wall (), t_hd;
    // Parallel classpoly only pays off when H_D is large: each worker repeats the
    // setup (class group, phi loads, ECRT precomputation), so for small h(D) the
    // single-job path is faster.  h(D) ~ sqrt(|D|), so gate on |D|.
    if ( jobs > 1 && d < 100000000UL ) jobs = 1;
    if ( jobs > 1 ) {
        // parallel H_D: run `jobs` classpoly ECRT workers (each covers prime indices
        // = jobid-1 mod jobs and dumps partial sums), then a merge pass writes outf.
        char edir[] = "/tmp/cmecrtXXXXXX";
        if ( ! mkdtemp (edir) ) { fprintf (stderr, "mkdtemp failed\n"); return 1; }
        snprintf (cmd, sizeof cmd,
            "export CLASSPOLY_JOBS=%d CLASSPOLY_ECRT_DIR=%s; "
            "for j in $(seq 1 %d); do CLASSPOLY_JOBID=$j classpoly %ld -1 %s %s.w$j -1 & done; wait; "
            "CLASSPOLY_JOBID=0 classpoly %ld -1 %s %s -1",
            jobs, edir, jobs, D, pdec, outf, D, pdec, outf);
        int rc = system (cmd);
        t_hd = wall () - t0;
        snprintf (cmd, sizeof cmd, "rm -rf %s", edir);  system (cmd);
        if ( rc != 0 ) { fprintf (stderr, "parallel classpoly failed\n"); return 1; }
    } else {
        snprintf (cmd, sizeof cmd, "classpoly %ld -1 %s %s -1", D, pdec, outf);
        if ( system (cmd) != 0 ) { fprintf (stderr, "classpoly failed\n"); return 1; }
        t_hd = wall () - t0;
    }
    FILE *f = fopen (outf, "r");  if ( ! f ) { fprintf (stderr, "no %s\n", outf); return 1; }
    int inv = -2, cap = 64, deg = -1;  mpz_t *hc = malloc (cap*sizeof(mpz_t));  char ln[8192];
    while ( fgets (ln, sizeof ln, f) ) {
        if ( ln[0] == 'I' && ln[1] == '=' ) { inv = atoi (ln+2); continue; }
        char *x = strstr (ln, "*X^");  if ( ! x ) continue;
        int e = atoi (x + 3);
        if ( e >= cap ) { while ( e >= cap ) cap *= 2;  hc = realloc (hc, cap*sizeof(mpz_t)); }
        while ( deg < e ) { deg++; mpz_init (hc[deg]); }
        *x = 0;  mpz_set_str (hc[e], ln, 10);
    }
    fclose (f);
    if ( inv == -2 || deg < 1 ) { fprintf (stderr, "could not parse H_D^inv\n"); return 1; }

    // 3. root f0 of H_D^inv over F_p
    fp_ctx C;  fp_init (&C, p);
    fp_poly h;  fpoly_init (&C, &h, deg + 2);
    for ( int i = 0 ; i <= deg ; i++ ) { mpz_mod (hc[i], hc[i], p);  fp_set_mpz (&C, h.c + (size_t)i*C.s, hc[i]); }
    h.deg = deg;
    t0 = wall ();
    mp_limb_t f0m[64];  fp_find_root (&C, &h, f0m, 0xC0FFEE + (unsigned) seed);
    double t_root = wall () - t0;
    mpz_t f0;  mpz_init (f0);  fp_get_mpz (&C, f0, f0m);

    // 4. f0 -> j candidates, and verify each has CM by D (trace +/- t)
    mpz_t J[8];  for ( int i = 0 ; i < 8 ; i++ ) mpz_init (J[i]);
    int nj = cm_j_from_inv (J, 8, f0, p, inv, phidir);
    if ( nj <= 0 ) { fprintf (stderr, "cm_j_from_inv failed (inv=%d)\n", inv); return 1; }

    mpz_t A4, A6, k, tmp;  mpz_inits (A4, A6, k, tmp, NULL);
    int found = 0, tsign = 0;  mpz_t jout;  mpz_init (jout);
    for ( int i = 0 ; i < nj && ! found ; i++ ) {
        // curve with j-invariant J[i]:  y^2 = x^3 + 3k x + 2k,  k = j/(1728-j)
        mpz_ui_sub (tmp, 1728, J[i]);  mpz_mod (tmp, tmp, p);
        if ( mpz_sgn (tmp) == 0 ) continue;                    // j=1728 (D=-4), excluded
        mpz_invert (tmp, tmp, p);  mpz_mul (k, J[i], tmp);  mpz_mod (k, k, p);
        mpz_mul_ui (A4, k, 3);  mpz_mod (A4, A4, p);
        mpz_mul_ui (A6, k, 2);  mpz_mod (A6, A6, p);
        int s = trace_sign (A4, A6, p, t, &cc);
        if ( s ) { found = 1;  tsign = s;  mpz_set (jout, J[i]); }
    }
    if ( ! found ) { fprintf (stderr, "no CM-by-D j among %d candidates (inv=%d)\n", nj, inv); return 3; }

    // 5. output
    gmp_printf ("D %ld\n", D);
    gmp_printf ("p %Zd\n", p);
    gmp_printf ("inv %d\n", inv);
    gmp_printf ("t %s%Zd\n", tsign > 0 ? "" : "-", t);         // Frobenius trace (#E = p+1-trace)
    if ( want_v ) gmp_printf ("v %Zd\n", v);
    gmp_printf ("j %Zd\n", jout);
    fprintf (stderr, "[H_D^inv(deg %d, inv %d): %.3fs   root: %.3fs   %d j-candidate(s)]\n", deg, inv, t_hd, t_root, nj);

    for ( int i = 0 ; i <= deg ; i++ ) mpz_clear (hc[i]);
    free (hc);
    for ( int i = 0 ; i < 8 ; i++ ) mpz_clear (J[i]);
    mpz_clears (A4, A6, k, tmp, jout, t, v, f0, NULL);
    fpoly_clear (&h);  fp_clear (&C);  cornacchia_clear (&cc);  mpz_clear (p);
    return 0;
}
