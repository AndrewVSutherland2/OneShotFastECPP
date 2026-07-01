#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "cornacchia.h"
#include "smooth.h"
#include "fproot.h"
#include "curve.h"

/*
    oneshot p : produce an n^4-smooth one-shot ECPP certificate for the prime p.

      oneshot (p=<dec> | pbits=<n> [seed=<s>]) [threads=<t>] [B=<scan bound>] [pcache=<file>]

    Pipeline (all C; shells out to the sibling tools dscan, cm_method, classpoly):
      (a) dscan finds fundamental discriminants D with 4p = t^2 + |D| v^2 solvable;
      (b) keep candidate orders N = p+1 -/+ t with N = 0 mod 4 (Montgomery-model);
      (c) batch-test the n^4-smooth part of each N; a winner has smooth part > L,
          giving a smooth m | N with L < m < L*r;
      (d) cm_method (best class invariant) gives a CM j-invariant, and mont_assemble
          builds the Montgomery (A, x0) with a point of order m.
    Prints  p A x0 m q1 ... qk  (verifiable by voneshot.py).  Needs the classpoly
    environment (source setenv.sh) and dscan/cm_method on PATH.  The prime-product
    P = prod_{q<=n^4} q is cached to pcache (default /tmp/oneshot_P_<y>.bin).
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

// cm_method D p -> j0 (best class invariant, converted to a j-invariant)
static int cm_jinvariant (long D, const char *pdec, mpz_t j0)
{
    char cmd[16384];
    snprintf (cmd, sizeof cmd, "cm_method D=%ld p=%s 2>/dev/null", D, pdec);
    FILE *pp = popen (cmd, "r");  if ( ! pp ) return 0;
    char line[8192];  int got = 0;
    while ( fgets (line, sizeof line, pp) )
        if ( line[0] == 'j' && line[1] == ' ' ) { mpz_set_str (j0, line + 2, 10);  got = 1; }
    pclose (pp);  return got;
}

int main (int argc, char **argv)
{
    mpz_t p;  mpz_init (p);
    unsigned long pbits = 0, seed = 1;  int nth = 0;  unsigned long B = 0;
    const char *pcache = NULL;
    for ( int i = 1 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i]+2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i]+6, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i]+5, 0, 10);
        else if ( ! strncmp (argv[i], "threads=", 8) ) nth = atoi (argv[i]+8);
        else if ( ! strncmp (argv[i], "B=", 2) ) B = strtoul (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "pcache=", 7) ) pcache = argv[i]+7;
    }
    if ( ! mpz_sgn (p) ) {
        if ( ! pbits ) pbits = 256;
        gmp_randstate_t rs; gmp_randinit_default (rs); gmp_randseed_ui (rs, seed);
        mpz_urandomb (p, rs, pbits); mpz_setbit (p, pbits-1); mpz_nextprime (p, p); gmp_randclear (rs);
    }
    if ( ! mpz_probab_prime_p (p, 25) ) { fprintf (stderr, "p is not a probable prime\n"); return 1; }
    if ( ! getenv ("CLASSPOLY_PHI_DIR") ) { fprintf (stderr, "set the classpoly env (source setenv.sh)\n"); return 1; }
    char pdec[8192];  gmp_snprintf (pdec, sizeof pdec, "%Zd", p);

    mpz_t L, Hass;  mpz_init (L);  mpz_init (Hass);
    unsigned long n;  uint64_t n2, n4;
    cert_bounds (p, L, Hass, &n, &n2, &n4);
    if ( ! B ) B = 30000000UL;
    fprintf (stderr, "p: %lu bits, n=%lu, n^4=%lu, B=%lu\n", (unsigned long) mpz_sizeinbase (p, 2), n, (unsigned long) n4, B);

    // Prime-product cache policy: a power-of-2 ladder, so one cached P serves a
    // whole octave of bit-lengths instead of one file per n.  Preference order:
    //   (1) an explicit pcache= file (any y);
    //   (2) an exact-y cache /tmp/oneshot_P_<n4>.bin (legacy naming);
    //   (3) the smallest CACHED power of two >= n4 (superset: primes above n4 are
    //       skipped by build_m when assembling m -- never built fresh);
    //   (4) y' = 2^floor(log2 n4) <= n4, cached or freshly built+saved -- the
    //       cheapest build (<= the exact-n4 cost).  The uncovered gap (y', n4] is
    //       recovered by near-miss top-up: bounded Pollard rho on the cofactor
    //       (a missing prime is <= n4 ~ 2^34, ~2^17 mulmods to find).
    // Cache directory: ONESHOT_PCACHE_DIR if set (setenv.sh points it at
    // <repo>/work/pcache so caches survive reboots), else /tmp.
    const char *cdir = getenv ("ONESHOT_PCACHE_DIR");  if ( ! cdir ) cdir = "/tmp";
    char cachebuf[512];
    smooth_base sb;  double t0 = wall ();  int have = 0;
    if ( pcache ) have = smooth_base_load (&sb, pcache);
    if ( ! have ) {
        snprintf (cachebuf, sizeof cachebuf, "%s/oneshot_P_%lu.bin", cdir, (unsigned long) n4);
        have = smooth_base_load (&sb, cachebuf);                     // (2) exact
    }
    int j0 = 63;  while ( ! ((n4 >> j0) & 1) ) j0--;                 // floor(log2 n4)
    for ( int j = ((1ULL << j0) == n4 ? j0 : j0 + 1) ; ! have && j <= j0 + 2 ; j++ ) {
        snprintf (cachebuf, sizeof cachebuf, "%s/oneshot_P_%llu.bin", cdir, 1ULL << j);
        have = smooth_base_load (&sb, cachebuf);                     // (3) cached >= n4
    }
    if ( ! have ) {
        uint64_t yp = 1ULL << j0;                                    // (4) round down
        snprintf (cachebuf, sizeof cachebuf, "%s/oneshot_P_%llu.bin", cdir, (unsigned long long) yp);
        if ( ! smooth_base_load (&sb, cachebuf) ) {
            fprintf (stderr, "building prime-product P (y=2^%d=%llu, ~%.0f MiB, one-time)...\n",
                     j0, (unsigned long long) yp, 1.44*yp/8/1048576.0);
            smooth_base_build (&sb, yp, nth);  smooth_base_save (&sb, cachebuf);
        }
    }
    if ( ! smooth_base_selfcheck (&sb) ) { fprintf (stderr, "P self-check failed (corrupt cache?)\n"); return 1; }
    int topup = (sb.y < n4);                                         // gap (y', n4] handled by rho
    fprintf (stderr, "P ready (y=%lu%s, %.1fs)\n", (unsigned long) sb.y,
             topup ? ", near-miss top-up on" : (sb.y > n4 ? ", oversize primes trimmed" : ""), wall () - t0);

    // (a,b) dscan -> solvable (d,t,v); keep N = p+1-/+t with N = 0 mod 4
    char cmd[16384];
    snprintf (cmd, sizeof cmd, "dscan p=%s B=%lu threads=%d dump 2>/dev/null", pdec, B, nth > 0 ? nth : 16);
    double t_scan = wall ();
    FILE *ds = popen (cmd, "r");  if ( ! ds ) { fprintf (stderr, "dscan failed (on PATH?)\n"); return 1; }
    size_t cap = 4096, nc = 0;
    unsigned long *cd = malloc (cap * sizeof(unsigned long));
    mpz_t *cN = malloc (cap * sizeof(mpz_t)), *ct = malloc (cap * sizeof(mpz_t)), *cNp = malloc (cap * sizeof(mpz_t));
    mpz_t t, v, np1, Nm;  mpz_inits (t, v, np1, Nm, NULL);  mpz_add_ui (np1, p, 1);
    unsigned long d;
    while ( gmp_fscanf (ds, "%lu %Zd %Zd", &d, t, v) == 3 ) {
        for ( int sgn = 0 ; sgn < 2 ; sgn++ ) {
            if ( sgn == 0 ) mpz_sub (Nm, np1, t);  else mpz_add (Nm, np1, t);
            if ( mpz_fdiv_ui (Nm, 4) != 0 ) continue;
            if ( nc == cap ) { cap*=2; cd=realloc(cd,cap*sizeof(unsigned long)); cN=realloc(cN,cap*sizeof(mpz_t)); ct=realloc(ct,cap*sizeof(mpz_t)); cNp=realloc(cNp,cap*sizeof(mpz_t)); }
            cd[nc] = d;  mpz_init_set (cN[nc], Nm);  mpz_init_set (ct[nc], t);  mpz_init (cNp[nc]);  nc++;
        }
    }
    pclose (ds);
    t_scan = wall () - t_scan;
    fprintf (stderr, "candidates (N=0 mod 4): %zu   [scan %.2fs]\n", nc, t_scan);

    // (c) batch y-smooth part (y = sb.y)
    mpz_t *S = malloc (nc * sizeof(mpz_t));  for ( size_t i=0;i<nc;i++ ) mpz_init (S[i]);
    double t_sm = wall ();
    smooth_parts (&sb, (const mpz_t*) cN, nc, S, nth);
    t_sm = wall () - t_sm;

    // near-miss top-up (only when sb.y < n4): a candidate whose y-smooth part
    // falls short of L by less than n4^2 could be pushed over by one or two
    // primes in (y, n4]; find them with bounded rho on the cofactor.
    int nnm = 0, nconv = 0;
    if ( topup ) {
        mpz_t nmband;  mpz_init (nmband);
        mpz_fdiv_q_ui (nmband, L, (unsigned long) n4);  mpz_fdiv_q_ui (nmband, nmband, (unsigned long) n4);
        double t_nm = wall ();
        #pragma omp parallel for num_threads(nth > 0 ? nth : 8) schedule(dynamic,8) reduction(+:nnm,nconv)
        for ( size_t i = 0 ; i < nc ; i++ ) {
            if ( mpz_cmp (S[i], L) > 0 || mpz_cmp (S[i], nmband) <= 0 ) continue;
            nnm++;
            smooth_topup (S[i], cN[i], n4);
            if ( mpz_cmp (S[i], L) > 0 ) nconv++;
        }
        mpz_clear (nmband);
        fprintf (stderr, "near-miss top-up: %d rechecked, %d converted to winners [%.2fs]\n",
                 nnm, nconv, wall () - t_nm);
    }

    // (d) winners -> cm_method -> assemble.  Process winners smallest |D| first:
    // h(D) ~ sqrt(|D|) drives both the H_D computation and the degree-h(D)
    // root-find, so a small-|D| winner is vastly cheaper (and gives a nicer cert).
    cornacchia_ctx cc;  cornacchia_init (&cc, p);
    fp_ctx C;  fp_init (&C, p);
    mpz_t A, x0, jj, mm;  mpz_inits (A, x0, jj, mm, NULL);
    uint64_t qs[64];  int nq, done = 0, nwin = 0, ntry = 0;
    size_t *order = malloc (nc * sizeof(size_t));
    size_t nord = 0;
    for ( size_t i = 0 ; i < nc ; i++ ) if ( mpz_cmp (S[i], L) > 0 ) order[nord++] = i;
    for ( size_t a = 0 ; a + 1 < nord ; a++ )                  // tiny insertion sort by |D|
        for ( size_t b = a + 1 ; b < nord ; b++ )
            if ( cd[order[b]] < cd[order[a]] ) { size_t t2 = order[a]; order[a] = order[b]; order[b] = t2; }
    for ( size_t k = 0 ; k < nord && ! done ; k++ ) {
        size_t i = order[k];
        nwin++;
        if ( ! build_m (mm, qs, &nq, S[i], L, n2, n4) ) continue;
        ntry++;
        if ( ! cm_jinvariant (-(long) cd[i], pdec, jj) ) continue;
        if ( ! mont_assemble (&C, &cc, jj, cN[i], ct[i], mm, A, x0, 0xA55E + (unsigned) cd[i]) ) continue;
        gmp_printf ("%Zd %Zd %Zd %Zd", p, A, x0, mm);
        for ( int q = 0 ; q < nq ; q++ ) printf (" %lu", (unsigned long) qs[q]);
        printf ("\n");
        fprintf (stderr, "[cert from D=-%lu after %d assembled candidate(s); smooth %.2fs, %d smoothness winners]\n",
                 cd[i], ntry, t_sm, nwin);
        done = 1;
    }
    if ( ! done ) fprintf (stderr, "no certificate (%d smoothness winners, none Montgomery-compatible); increase B\n", nwin);

    for ( size_t i=0;i<nc;i++ ) { mpz_clear (cN[i]); mpz_clear (ct[i]); mpz_clear (cNp[i]); mpz_clear (S[i]); }
    free (cd); free (cN); free (ct); free (cNp); free (S); free (order);
    mpz_clears (t, v, np1, Nm, A, x0, jj, mm, L, Hass, p, NULL);
    smooth_base_clear (&sb);  cornacchia_clear (&cc);  fp_clear (&C);
    return done ? 0 : 2;
}
