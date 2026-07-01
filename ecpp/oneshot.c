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

      oneshot (p=<dec> | pbits=<n> [seed=<s>]) [threads=<t>] [c=<ratio>]
              [B0=<initial scan bound>] [B=<max scan bound>] [pcache=<file>]

    Adaptive pipeline (all C; shells out to dscan, cm_method, classpoly): instead
    of building the full prime product P = prod_{q<=n^4} q up front (which used to
    dominate fresh runs), the smoothness bound climbs a power-of-2 ladder starting
    just above n^2, and the discriminant-scan bound B doubles independently:

      - segments  seg_j = prod of primes in (y_{j-1}, y_j]  (y_0 ~ 2^ceil(log2 n^2),
        doubling, final segment capped exactly at n^4) are built on demand and
        cached individually -- segments are bit-length independent, so one cache
        serves every prime size whose ladder passes through it;
      - every candidate order N keeps a running smooth part S = prod of its
        per-segment smooth parts; winners (S > L) can appear at any rung, and small
        rungs often suffice (certificates rarely need primes anywhere near n^4);
      - the ladder deepens only when the cumulative testing work justifies the
        next segment (work >= c * (P bits so far + next segment bits), c=1 by
        default); otherwise the candidate pool widens (incremental dscan chunk
        via Bmin=).  Both growths are geometric, so the total work is within a
        small constant of the best fixed choice in hindsight, and a cold first
        run costs about the same as a warm one.

    Prints  p A x0 m q1 ... qk  (verifiable by voneshot.py).  Needs the classpoly
    environment (source setenv.sh).  Segment caches live in $ONESHOT_PCACHE_DIR
    (default /tmp); a legacy full-product cache (oneshot_P_<y>.bin) is used as the
    ladder base when present.
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

#define MAXSEG 32

typedef struct {
    // candidate pool
    unsigned long *cd;          // |D|
    mpz_t *cN, *ct, *S;         // order, trace, running smooth part
    unsigned char *dead;        // 1 = assembly failed, skip forever
    size_t nc, cap;
    // ladder
    smooth_base seg[MAXSEG];
    int nseg;
    uint64_t ycur;
    double Pbits, Wtest;        // work accounting (bits)
} engine;

// test pool[idx[0..k)] against segment s; S *= per-segment smooth part
static void test_batch (engine *E, int s, const size_t *idx, size_t k, int nth)
{
    if ( ! k ) return;
    mpz_t *NN = malloc (k * sizeof(mpz_t)), *SS = malloc (k * sizeof(mpz_t));
    for ( size_t i = 0 ; i < k ; i++ ) { mpz_init_set (NN[i], E->cN[idx[i]]);  mpz_init (SS[i]); }
    smooth_parts (&E->seg[s], (const mpz_t*) NN, k, SS, nth);
    for ( size_t i = 0 ; i < k ; i++ ) {
        mpz_mul (E->S[idx[i]], E->S[idx[i]], SS[i]);
        mpz_clear (NN[i]);  mpz_clear (SS[i]);
    }
    free (NN);  free (SS);
    E->Wtest += mpz_sizeinbase (E->seg[s].P, 2) + (double) k * mpz_sizeinbase (E->cN[idx[0]], 2);
}

// acquire the ladder segment (E->ycur, ynext]: cache-load or build+save
static void add_segment (engine *E, uint64_t ynext, const char *cdir, int nth)
{
    smooth_base *sb = &E->seg[E->nseg];
    char path[512];
    snprintf (path, sizeof path, "%s/oneshot_Pseg_%llu_%llu.bin", cdir,
              (unsigned long long) E->ycur, (unsigned long long) ynext);
    double t0 = wall ();
    if ( ! (smooth_base_load (sb, path) && sb->lo == E->ycur && sb->y == ynext && smooth_base_selfcheck (sb)) ) {
        smooth_base_build_range (sb, E->ycur, ynext, nth);
        smooth_base_save (sb, path);
        fprintf (stderr, "segment (%llu, %llu]: built %.0f Mbit (%.1fs)\n",
                 (unsigned long long) E->ycur, (unsigned long long) ynext,
                 mpz_sizeinbase (sb->P, 2)/1e6, wall () - t0);
    }
    E->Pbits += mpz_sizeinbase (sb->P, 2);
    E->ycur = ynext;
    E->nseg++;
}

int main (int argc, char **argv)
{
    mpz_t p;  mpz_init (p);
    unsigned long pbits = 0, seed = 1;  int nth = 0;
    unsigned long B0 = 4000000UL, Bmax = 20000000000UL;
    double cfac = 1.0;
    const char *pcache = NULL;
    for ( int i = 1 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i]+2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i]+6, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i]+5, 0, 10);
        else if ( ! strncmp (argv[i], "threads=", 8) ) nth = atoi (argv[i]+8);
        else if ( ! strncmp (argv[i], "B0=", 3) ) B0 = strtoul (argv[i]+3, 0, 10);
        else if ( ! strncmp (argv[i], "B=", 2) ) Bmax = strtoul (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "c=", 2) ) cfac = strtod (argv[i]+2, 0);
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
    const char *cdir = getenv ("ONESHOT_PCACHE_DIR");  if ( ! cdir ) cdir = "/tmp";
    fprintf (stderr, "p: %lu bits, n^2=%llu, n^4=%llu, B0=%lu, Bmax=%lu, c=%.2f\n",
             (unsigned long) mpz_sizeinbase (p, 2), (unsigned long long) n2, (unsigned long long) n4, B0, Bmax, cfac);

    engine E;  memset (&E, 0, sizeof E);
    E.cap = 8192;
    E.cd = malloc (E.cap * sizeof(unsigned long));
    E.cN = malloc (E.cap * sizeof(mpz_t));  E.ct = malloc (E.cap * sizeof(mpz_t));
    E.S  = malloc (E.cap * sizeof(mpz_t));  E.dead = malloc (E.cap);

    // ladder base: explicit pcache= file, else the largest cached legacy full
    // product P(2^j) (j descending), else start empty just above n^2.
    if ( pcache && smooth_base_load (&E.seg[0], pcache) && smooth_base_selfcheck (&E.seg[0]) ) {
        E.Pbits = mpz_sizeinbase (E.seg[0].P, 2);  E.ycur = E.seg[0].y;  E.nseg = 1;
        fprintf (stderr, "base P: y=%llu from %s\n", (unsigned long long) E.ycur, pcache);
    } else {
        int jtop = 0;  while ( (1ULL << (jtop+1)) < n4 ) jtop++;      // 2^jtop <= n4 < 2^(jtop+2)
        for ( int j = jtop + 1 ; j >= 0 && ! E.nseg ; j-- ) {
            uint64_t y = 1ULL << j;
            if ( y < n2 ) break;
            char path[512];
            snprintf (path, sizeof path, "%s/oneshot_P_%llu.bin", cdir, (unsigned long long) y);
            if ( smooth_base_load (&E.seg[0], path) && E.seg[0].lo == 0 && smooth_base_selfcheck (&E.seg[0]) ) {
                E.Pbits = mpz_sizeinbase (E.seg[0].P, 2);  E.ycur = E.seg[0].y;  E.nseg = 1;
                fprintf (stderr, "base P: y=%llu (legacy cache)\n", (unsigned long long) E.ycur);
            }
        }
    }
    if ( ! E.nseg ) {                                       // first rung: (0, 2^ceil(log2 n^2)]
        uint64_t y0 = 1;  while ( y0 < n2 ) y0 <<= 1;
        if ( y0 > n4 ) y0 = n4;
        add_segment (&E, y0, cdir, nth);
    }

    cornacchia_ctx cc;  cornacchia_init (&cc, p);
    fp_ctx C;  fp_init (&C, p);
    mpz_t A, x0, jj, mm, t, v, np1, Nm;  mpz_inits (A, x0, jj, mm, t, v, np1, Nm, NULL);
    mpz_add_ui (np1, p, 1);
    uint64_t qs[64];  int nq;
    size_t *idx = NULL;  size_t idxcap = 0;
    unsigned long Bcur = 0;
    int done = 0, epoch = 0, nwin_total = 0;
    double t_start = wall ();

    while ( ! done ) {
        epoch++;
        // ---- winners so far, smallest |D| first ----
        size_t nwin = 0;
        if ( idxcap < E.nc ) { idxcap = E.nc + 16;  idx = realloc (idx, idxcap * sizeof(size_t)); }
        for ( size_t i = 0 ; i < E.nc ; i++ )
            if ( ! E.dead[i] && mpz_cmp (E.S[i], L) > 0 ) idx[nwin++] = i;
        for ( size_t a = 0 ; a + 1 < nwin ; a++ )
            for ( size_t b = a + 1 ; b < nwin ; b++ )
                if ( E.cd[idx[b]] < E.cd[idx[a]] ) { size_t w = idx[a]; idx[a] = idx[b]; idx[b] = w; }
        for ( size_t k = 0 ; k < nwin && ! done ; k++ ) {
            size_t i = idx[k];
            nwin_total++;
            if ( ! build_m (mm, qs, &nq, E.S[i], L, n2, n4) ) continue;       // may succeed at a higher rung
            if ( ! cm_jinvariant (-(long) E.cd[i], pdec, jj) ) { E.dead[i] = 1; continue; }
            if ( ! mont_assemble (&C, &cc, jj, E.cN[i], E.ct[i], mm, A, x0, 0xA55E + (unsigned) E.cd[i]) )
                { E.dead[i] = 1; continue; }
            gmp_printf ("%Zd %Zd %Zd %Zd", p, A, x0, mm);
            for ( int q = 0 ; q < nq ; q++ ) printf (" %lu", (unsigned long) qs[q]);
            printf ("\n");
            fprintf (stderr, "[cert: D=-%lu at y=%llu, B=%lu, pool=%zu, %.1fs total]\n",
                     E.cd[i], (unsigned long long) E.ycur, Bcur, E.nc, wall () - t_start);
            done = 1;
        }
        if ( done ) break;

        // ---- expand: deepen the ladder when tested work justifies it, else widen ----
        uint64_t ynext = E.ycur < n4 ? (E.ycur << 1 >= n4 ? n4 : E.ycur << 1) : 0;
        double estnext = ynext ? 1.443 * (double)(ynext - E.ycur) : 0;
        int deepen = ynext && (E.Wtest >= cfac * (E.Pbits + estnext) || Bcur >= Bmax);
        if ( ! deepen && Bcur >= Bmax && ! ynext ) {
            fprintf (stderr, "no certificate: pool %zu candidates (B=%lu), y=n^4; raise B=\n", E.nc, Bcur);
            break;
        }
        if ( deepen ) {
            add_segment (&E, ynext, cdir, nth);
            size_t k = 0;
            for ( size_t i = 0 ; i < E.nc ; i++ ) if ( ! E.dead[i] ) idx[k++] = i;
            test_batch (&E, E.nseg - 1, idx, k, nth);
            fprintf (stderr, "[deepen: y=%llu, pool=%zu]\n", (unsigned long long) E.ycur, E.nc);
        } else {
            unsigned long Bnext = Bcur ? Bcur * 2 : B0;
            if ( Bnext > Bmax ) Bnext = Bmax;
            char cmd[16384];
            snprintf (cmd, sizeof cmd, "dscan p=%s B=%lu Bmin=%lu threads=%d dump 2>/dev/null",
                      pdec, Bnext, Bcur, nth > 0 ? nth : 16);
            FILE *ds = popen (cmd, "r");
            if ( ! ds ) { fprintf (stderr, "dscan failed (on PATH?)\n"); return 1; }
            size_t first_new = E.nc;
            unsigned long d;
            while ( gmp_fscanf (ds, "%lu %Zd %Zd", &d, t, v) == 3 ) {
                for ( int sgn = 0 ; sgn < 2 ; sgn++ ) {
                    if ( sgn == 0 ) mpz_sub (Nm, np1, t);  else mpz_add (Nm, np1, t);
                    if ( mpz_fdiv_ui (Nm, 4) != 0 ) continue;          // Montgomery needs N = 0 mod 4
                    if ( E.nc == E.cap ) {
                        E.cap *= 2;
                        E.cd = realloc (E.cd, E.cap*sizeof(unsigned long));
                        E.cN = realloc (E.cN, E.cap*sizeof(mpz_t));  E.ct = realloc (E.ct, E.cap*sizeof(mpz_t));
                        E.S  = realloc (E.S,  E.cap*sizeof(mpz_t));  E.dead = realloc (E.dead, E.cap);
                    }
                    E.cd[E.nc] = d;  mpz_init_set (E.cN[E.nc], Nm);  mpz_init_set (E.ct[E.nc], t);
                    mpz_init_set_ui (E.S[E.nc], 1);  E.dead[E.nc] = 0;  E.nc++;
                }
            }
            pclose (ds);
            Bcur = Bnext;
            // test the new candidates against every rung of the ladder
            size_t knew = E.nc - first_new;
            if ( idxcap < E.nc ) { idxcap = E.nc + 16;  idx = realloc (idx, idxcap * sizeof(size_t)); }
            for ( size_t i = 0 ; i < knew ; i++ ) idx[i] = first_new + i;
            for ( int s = 0 ; s < E.nseg ; s++ ) test_batch (&E, s, idx, knew, nth);
            fprintf (stderr, "[widen: B=%lu, +%zu candidates (pool %zu), y=%llu]\n",
                     Bcur, knew, E.nc, (unsigned long long) E.ycur);
        }
    }

    for ( size_t i = 0 ; i < E.nc ; i++ ) { mpz_clear (E.cN[i]); mpz_clear (E.ct[i]); mpz_clear (E.S[i]); }
    free (E.cd); free (E.cN); free (E.ct); free (E.S); free (E.dead); free (idx);
    for ( int s = 0 ; s < E.nseg ; s++ ) smooth_base_clear (&E.seg[s]);
    mpz_clears (A, x0, jj, mm, t, v, np1, Nm, L, Hass, p, NULL);
    cornacchia_clear (&cc);  fp_clear (&C);
    return done ? 0 : 2;
}
