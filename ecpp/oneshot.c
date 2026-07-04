#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
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
    (default /tmp); pcache=<file> seeds the ladder with a prebuilt full product
    (smoothtest pbuild save=).
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

// cm_method D p -> j0 (best class invariant, converted to a j-invariant).  ells
// (may be "") lists the primes ell | n1 to clear by ell-volcano descent: the
// floor curve cm_method then returns has cyclic ell-Sylow for each one, so it is
// Montgomery-representable and a point of order m exists whenever m | N/n1h.
static int cm_jinvariant (long D, const char *pdec, int jobs, const char *ells, mpz_t j0)
{
    char cmd[16384];
    snprintf (cmd, sizeof cmd, "cm_method D=%ld p=%s jobs=%d%s%s 2>/dev/null", D, pdec,
              jobs > 0 ? jobs : 1, *ells ? " ells=" : "", ells);
    FILE *pp = popen (cmd, "r");  if ( ! pp ) return 0;
    char line[8192];  int got = 0;
    while ( fgets (line, sizeof line, pp) )
        if ( line[0] == 'j' && line[1] == ' ' ) { mpz_set_str (j0, line + 2, 10);  got = 1; }
    pclose (pp);  return got;
}

// the primes with a classical Phi_ell in the phi-file bundle (all ell <= 97):
// cm_method can descend these to the volcano floor, so only larger primes of n1
// constrain the certificate m.
static const unsigned long ELLS[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
#define NELLS (sizeof(ELLS)/sizeof(ELLS[0]))

#define MAXSEG 32

/*
    Asynchronous discriminant scan: dscan runs as a forked child while the
    smoothness engine works on the current pool, and a drain thread parses its
    dump output and precomputes the per-candidate group data (N, n1, n1h) off
    the critical path.  The next window (Bcur, 2*Bcur] is launched speculatively
    right after each ingest -- the window is fixed by the doubling schedule, so
    the result is needed at the next widen regardless of intervening deepens
    (it is only wasted if a certificate lands first, in which case the child is
    killed).
*/
typedef struct { unsigned long d; mpz_t N, t, n1, n1h; } scand;

typedef struct {
    FILE *f;  pid_t pid;
    unsigned long Bmin, B;      // scan window (Bmin, B]
    pthread_t th;
    scand *v;  size_t n, cap;   // parsed candidates (both signs, N=0 mod 4 only)
    mpz_t p;
} scanjob;

static void *scan_drain (void *arg)
{
    scanjob *J = arg;
    mpz_t t, v, np1, Nm, a2;
    mpz_inits (t, v, np1, Nm, a2, NULL);
    mpz_add_ui (np1, J->p, 1);
    unsigned long d;
    while ( gmp_fscanf (J->f, "%lu %Zd %Zd", &d, t, v) == 3 ) {
        for ( int sgn = 0 ; sgn < 2 ; sgn++ ) {
            if ( sgn == 0 ) mpz_sub (Nm, np1, t);  else mpz_add (Nm, np1, t);
            if ( mpz_fdiv_ui (Nm, 4) != 0 ) continue;          // Montgomery needs N = 0 mod 4
            // (No further representability filter: cm_method descends the
            // 2-volcano to a floor curve with cyclic 2-Sylow, so the
            // p = 3 mod 4, N = 4 mod 8 obstruction does not apply.)
            if ( J->n == J->cap ) { J->cap = J->cap ? 2*J->cap : 8192;  J->v = realloc (J->v, J->cap * sizeof(scand)); }
            scand *c = &J->v[J->n++];
            c->d = d;
            mpz_init_set (c->N, Nm);  mpz_init_set (c->t, t);
            // group structure Z/n1 x Z/(N/n1):  pi = (t + v sqrt(D))/2, n1 = gcd(a-1, b)
            // over O = Z[omega].  For trace +t (N = p+1-t):
            //   D = 0 mod 4:  n1 = gcd(t/2 - 1, v);   D = 1 mod 4:  n1 = gcd((t-v)/2 - 1, v)
            // and t -> -t for the other sign.  (Validated vs PARI ellgroup, D < -4.)
            if ( (d & 3) == 0 ) {                              // D = -d = 0 mod 4
                mpz_fdiv_q_2exp (a2, t, 1);                    // t/2
                if ( sgn == 0 ) mpz_sub_ui (a2, a2, 1); else { mpz_neg (a2, a2); mpz_sub_ui (a2, a2, 1); }
            } else {                                           // D = 1 mod 4
                if ( sgn == 0 ) mpz_sub (a2, t, v); else { mpz_neg (a2, t); mpz_sub (a2, a2, v); }
                mpz_fdiv_q_2exp (a2, a2, 1);
                mpz_sub_ui (a2, a2, 1);
            }
            mpz_init (c->n1);  mpz_gcd (c->n1, a2, v);
            // n1h = n1 stripped of the volcano-descendable primes (<= 97)
            mpz_init_set (c->n1h, c->n1);
            for ( size_t e = 0 ; e < NELLS ; e++ )
                while ( mpz_divisible_ui_p (c->n1h, ELLS[e]) )
                    mpz_divexact_ui (c->n1h, c->n1h, ELLS[e]);
        }
    }
    mpz_clears (t, v, np1, Nm, a2, NULL);
    return NULL;
}

static int scan_spawn (scanjob *J, const mpz_t p, const char *pdec, unsigned long B, unsigned long Bmin, int nth)
{
    memset (J, 0, sizeof *J);
    J->B = B;  J->Bmin = Bmin;
    int fd[2];
    if ( pipe (fd) ) return 0;
    J->pid = fork ();
    if ( J->pid < 0 ) { close (fd[0]); close (fd[1]); return 0; }
    if ( J->pid == 0 ) {
        char ap[8200], ab[32], abm[32], ath[16];
        snprintf (ap, sizeof ap, "p=%s", pdec);
        snprintf (ab, sizeof ab, "B=%lu", B);
        snprintf (abm, sizeof abm, "Bmin=%lu", Bmin);
        snprintf (ath, sizeof ath, "threads=%d", nth);
        dup2 (fd[1], 1);  close (fd[0]);  close (fd[1]);
        execlp ("dscan", "dscan", ap, ab, abm, ath, "dump", (char*) 0);
        _exit (127);
    }
    close (fd[1]);
    J->f = fdopen (fd[0], "r");
    mpz_init_set (J->p, p);
    pthread_create (&J->th, NULL, scan_drain, J);
    return 1;
}

// wait for scan + parser to finish; J->v/J->n are then owned by the caller
static int scan_join (scanjob *J)
{
    pthread_join (J->th, NULL);
    fclose (J->f);
    int st = 0;  waitpid (J->pid, &st, 0);
    mpz_clear (J->p);
    return WIFEXITED (st) && WEXITSTATUS (st) == 0;
}

static void scan_abort (scanjob *J)
{
    kill (J->pid, SIGTERM);
    pthread_join (J->th, NULL);
    fclose (J->f);
    int st = 0;  waitpid (J->pid, &st, 0);
    mpz_clear (J->p);
    for ( size_t i = 0 ; i < J->n ; i++ ) mpz_clears (J->v[i].N, J->v[i].t, J->v[i].n1, J->v[i].n1h, NULL);
    free (J->v);
}

typedef struct {
    // candidate pool
    unsigned long *cd;          // |D|
    mpz_t *cN, *ct, *S, *cj;    // order, trace, running smooth part, cached j (0 = none)
    mpz_t *n1;                  // E(F_p) = Z/n1 x Z/(N/n1): n1 = gcd(a-1, b), pi-1 = a+b*omega
    mpz_t *n1h;                 // the part of n1 at primes > 97: not clearable by volcano
                                // descent, so m must avoid it (m | N/n1h); usually 1
    unsigned char *dead;        // 1 = cm_method failed for this D, skip forever
    unsigned long *Sfail;       // bits of S at the last failed assembly (retry once S grows)
    size_t nc, cap;
    // ladder
    smooth_base seg[MAXSEG];
    int nseg;
    uint64_t ycur;
    double Pbits, Wtest;        // work accounting (bits)
} engine;

// test pool[idx[0..k)] against segments [s0,s1) in one fused pass;
// S *= product of the per-segment smooth parts
static void test_batch_range (engine *E, int s0, int s1, const size_t *idx, size_t k, int nth)
{
    if ( ! k || s1 <= s0 ) return;
    double t0 = wall ();
    mpz_t *NN = malloc (k * sizeof(mpz_t)), *SS = malloc (k * sizeof(mpz_t));
    for ( size_t i = 0 ; i < k ; i++ ) { mpz_init_set (NN[i], E->cN[idx[i]]);  mpz_init (SS[i]); }
    smooth_parts_multi (&E->seg[s0], s1 - s0, (const mpz_t*) NN, k, SS, nth);
    for ( size_t i = 0 ; i < k ; i++ ) {
        mpz_mul (E->S[idx[i]], E->S[idx[i]], SS[i]);
        mpz_clear (NN[i]);  mpz_clear (SS[i]);
    }
    free (NN);  free (SS);
    double bits = 0;
    for ( int s = s0 ; s < s1 ; s++ ) bits += mpz_sizeinbase (E->seg[s].P, 2);
    E->Wtest += bits + (double) k * (s1 - s0) * mpz_sizeinbase (E->cN[idx[0]], 2);
    fprintf (stderr, "[smooth segs %d..%d (%.0f Mbit) x %zu cand: %.1fs]\n",
             s0, s1 - 1, bits/1e6, k, wall () - t0);
}

static void test_batch (engine *E, int s, const size_t *idx, size_t k, int nth)
{ test_batch_range (E, s, s + 1, idx, k, nth); }

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


// Seff = S with the n1h-supported part removed: a point of order m needs
// m | exponent(E) = N/n1, and v_q(S) = v_q(N) for q <= y, so removing q^v_q(n1h)
// for the primes of n1h makes any m | Seff automatically divide the exponent.
// Only the >97 part n1h of n1 is removed: cm_method clears the small primes of
// n1 by descending their isogeny volcanoes to the floor (walk-to-the-floor with
// the classical Phi_ell), where the ell-Sylow subgroup is cyclic.
static void exponent_part (mpz_t Seff, mpz_t g, mpz_t r, const mpz_t S, const mpz_t n1)
{
    mpz_set (Seff, S);  mpz_set (r, n1);
    for (;;) {
        mpz_gcd (g, r, Seff);
        if ( mpz_cmp_ui (g, 1) == 0 ) break;
        mpz_divexact (Seff, Seff, g);
        mpz_divexact (r, r, g);
    }
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
    // default threads to all cores: the dscan/cm_method subprocesses need the
    // count passed explicitly, so resolve it here rather than leaving nth=0
    if ( nth <= 0 ) nth = omp_get_max_threads ();
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
    E.S  = malloc (E.cap * sizeof(mpz_t));  E.cj = malloc (E.cap * sizeof(mpz_t));
    E.n1 = malloc (E.cap * sizeof(mpz_t));  E.n1h = malloc (E.cap * sizeof(mpz_t));
    E.dead = malloc (E.cap);  E.Sfail = malloc (E.cap * sizeof(unsigned long));

    // ladder base: explicit pcache= file (a full product from smoothtest pbuild),
    // else start empty just above n^2.  (No automatic pickup of full-product
    // caches: a big base rung front-loads the whole testing cost per candidate
    // and defeats the ladder's amortization -- measured 3x slower at 299 bits.)
    if ( pcache && smooth_base_load (&E.seg[0], pcache) && smooth_base_selfcheck (&E.seg[0]) ) {
        E.Pbits = mpz_sizeinbase (E.seg[0].P, 2);  E.ycur = E.seg[0].y;  E.nseg = 1;
        fprintf (stderr, "base P: y=%llu from %s\n", (unsigned long long) E.ycur, pcache);
    }
    if ( ! E.nseg ) {                                       // first rung: (0, 2^ceil(log2 n^2)]
        uint64_t y0 = 1;  while ( y0 < n2 ) y0 <<= 1;
        if ( y0 > n4 ) y0 = n4;
        add_segment (&E, y0, cdir, nth);
    }

    cornacchia_ctx cc;  cornacchia_init (&cc, p);
    fp_ctx C;  fp_init (&C, p);
    mpz_t A, x0, jj, mm, Seff, gtmp, rtmp;
    mpz_inits (A, x0, jj, mm, Seff, gtmp, rtmp, NULL);
    uint64_t qs[64];  int nq;
    size_t *idx = NULL;  size_t idxcap = 0;
    unsigned long Bcur = 0;
    int done = 0, epoch = 0, nwin_total = 0;
    double t_start = wall ();

    // speculative first scan window (0, B0]
    scanjob pend;  int have_pend = 0;
    if ( B0 <= Bmax ) have_pend = scan_spawn (&pend, p, pdec, B0 > Bmax ? Bmax : B0, 0, nth);

    while ( ! done ) {
        epoch++;
        // ---- winners so far, smallest |D| first ----
        size_t nwin = 0;
        if ( idxcap < E.nc ) { idxcap = E.nc + 16;  idx = realloc (idx, idxcap * sizeof(size_t)); }
        for ( size_t i = 0 ; i < E.nc ; i++ ) {
            if ( E.dead[i] || mpz_cmp (E.S[i], L) <= 0 ) continue;
            if ( E.Sfail[i] && mpz_sizeinbase (E.S[i], 2) <= E.Sfail[i] ) continue;
            exponent_part (Seff, gtmp, rtmp, E.S[i], E.n1h[i]);
            if ( mpz_cmp (Seff, L) > 0 ) idx[nwin++] = i;
        }
        for ( size_t a = 0 ; a + 1 < nwin ; a++ )
            for ( size_t b = a + 1 ; b < nwin ; b++ )
                if ( E.cd[idx[b]] < E.cd[idx[a]] ) { size_t w = idx[a]; idx[a] = idx[b]; idx[b] = w; }

        // ---- winner polish: consuming a winner costs ~cm_method(h), h ~ sqrt(|D|).
        // While the best winner is expensive and one more rung is cheap, deepen
        // first -- a smaller-|D| winner may gate at the next rung and save far
        // more than the rung costs (measured: a 200-bit run once committed to an
        // h=6122 winner when h=162 was one rung away in the same pool).
        while ( nwin && E.ycur < n4 ) {
            double h_est = 0.35 * sqrt ((double) E.cd[idx[0]]);
            double t_cm  = pow (h_est / 2000.0, 1.6) * ((double) n / 256.0);      // ~seconds
            uint64_t yn2 = (E.ycur << 1 >= n4) ? n4 : E.ycur << 1;
            double t_rng = 1.443 * (double)(yn2 - E.ycur) / 1e8;                  // build+test ~seconds
            if ( t_cm <= 2 * t_rng ) break;
            fprintf (stderr, "[polish: best D=-%lu (h~%.0f) costs ~%.1fs, next rung ~%.1fs -> deepen]\n",
                     E.cd[idx[0]], h_est, t_cm, t_rng);
            add_segment (&E, yn2, cdir, nth);
            size_t kk = 0;
            for ( size_t i = 0 ; i < E.nc ; i++ ) if ( ! E.dead[i] ) idx[kk++] = i;
            test_batch (&E, E.nseg - 1, idx, kk, nth);
            nwin = 0;
            for ( size_t i = 0 ; i < E.nc ; i++ ) {
                if ( E.dead[i] || mpz_cmp (E.S[i], L) <= 0 ) continue;
                if ( E.Sfail[i] && mpz_sizeinbase (E.S[i], 2) <= E.Sfail[i] ) continue;
                exponent_part (Seff, gtmp, rtmp, E.S[i], E.n1h[i]);
                if ( mpz_cmp (Seff, L) > 0 ) idx[nwin++] = i;
            }
            for ( size_t a = 0 ; a + 1 < nwin ; a++ )
                for ( size_t b = a + 1 ; b < nwin ; b++ )
                    if ( E.cd[idx[b]] < E.cd[idx[a]] ) { size_t w = idx[a]; idx[a] = idx[b]; idx[b] = w; }
        }

        for ( size_t k = 0 ; k < nwin && ! done ; k++ ) {
            size_t i = idx[k];
            nwin_total++;
            exponent_part (Seff, gtmp, rtmp, E.S[i], E.n1h[i]);
            if ( ! build_m (mm, qs, &nq, Seff, L, n2, n4) ) continue;         // may succeed at a higher rung
            if ( ! mpz_sgn (E.cj[i]) ) {                    // j is candidate-specific: compute once, cache
                char ellstr[160];  size_t el = 0;  ellstr[0] = 0;
                for ( size_t e = 0 ; e < NELLS ; e++ )
                    if ( mpz_divisible_ui_p (E.n1[i], ELLS[e]) )
                        el += snprintf (ellstr + el, sizeof ellstr - el, "%s%lu", el ? "," : "", ELLS[e]);
                double t_cm0 = wall ();
                if ( ! cm_jinvariant (-(long) E.cd[i], pdec, nth, ellstr, jj) ) { E.dead[i] = 1; continue; }
                fprintf (stderr, "[cm_method D=-%lu: %.1fs]\n", E.cd[i], wall () - t_cm0);
                mpz_set (E.cj[i], jj);
            } else mpz_set (jj, E.cj[i]);
            if ( ! mont_assemble (&C, &cc, jj, E.cN[i], E.ct[i], mm, A, x0, 0xA55E + (unsigned) E.cd[i]) ) {
                // should be rare now that m | exponent(E) holds by construction;
                // retry if S ever grows (a different m becomes available)
                E.Sfail[i] = mpz_sizeinbase (E.S[i], 2);  continue;
            }
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
            double t_ds = wall ();
            if ( ! have_pend ) have_pend = scan_spawn (&pend, p, pdec, Bnext, Bcur, nth);
            if ( ! have_pend || pend.Bmin != Bcur || pend.B != Bnext ) { fprintf (stderr, "scan window mismatch\n"); return 1; }
            if ( ! scan_join (&pend) ) { fprintf (stderr, "dscan failed (on PATH?)\n"); return 1; }
            have_pend = 0;
            // bulk-append the parsed candidates (mpz ownership moves from pend.v)
            size_t first_new = E.nc;
            for ( size_t ci = 0 ; ci < pend.n ; ci++ ) {
                scand *c = &pend.v[ci];
                if ( E.nc == E.cap ) {
                    E.cap *= 2;
                    E.cd = realloc (E.cd, E.cap*sizeof(unsigned long));
                    E.cN = realloc (E.cN, E.cap*sizeof(mpz_t));  E.ct = realloc (E.ct, E.cap*sizeof(mpz_t));
                    E.S  = realloc (E.S,  E.cap*sizeof(mpz_t));  E.cj = realloc (E.cj, E.cap*sizeof(mpz_t));
                    E.n1 = realloc (E.n1, E.cap*sizeof(mpz_t));  E.n1h = realloc (E.n1h, E.cap*sizeof(mpz_t));
                    E.dead = realloc (E.dead, E.cap);  E.Sfail = realloc (E.Sfail, E.cap*sizeof(unsigned long));
                }
                E.cd[E.nc] = c->d;
                E.cN[E.nc][0] = c->N[0];  E.ct[E.nc][0] = c->t[0];
                E.n1[E.nc][0] = c->n1[0];  E.n1h[E.nc][0] = c->n1h[0];
                mpz_init_set_ui (E.S[E.nc], 1);  mpz_init (E.cj[E.nc]);
                E.dead[E.nc] = 0;  E.Sfail[E.nc] = 0;  E.nc++;
            }
            free (pend.v);
            Bcur = Bnext;
            fprintf (stderr, "[scan join B=%lu: %.1fs stall]\n", Bcur, wall () - t_ds);
            // launch the next window before testing: it scans while we test
            if ( Bcur < Bmax )
                have_pend = scan_spawn (&pend, p, pdec, Bcur * 2 > Bmax ? Bmax : Bcur * 2, Bcur, nth);
            // test the new candidates against every rung of the ladder (one fused pass)
            size_t knew = E.nc - first_new;
            if ( idxcap < E.nc ) { idxcap = E.nc + 16;  idx = realloc (idx, idxcap * sizeof(size_t)); }
            for ( size_t i = 0 ; i < knew ; i++ ) idx[i] = first_new + i;
            test_batch_range (&E, 0, E.nseg, idx, knew, nth);
            fprintf (stderr, "[widen: B=%lu, +%zu candidates (pool %zu), y=%llu]\n",
                     Bcur, knew, E.nc, (unsigned long long) E.ycur);
        }
    }
    if ( have_pend ) scan_abort (&pend);                    // don't leave a scan running

    for ( size_t i = 0 ; i < E.nc ; i++ ) { mpz_clear (E.cN[i]); mpz_clear (E.ct[i]); mpz_clear (E.S[i]); mpz_clear (E.cj[i]); mpz_clear (E.n1[i]); mpz_clear (E.n1h[i]); }
    free (E.cd); free (E.cN); free (E.ct); free (E.S); free (E.cj); free (E.n1); free (E.n1h); free (E.dead); free (E.Sfail); free (idx);
    for ( int s = 0 ; s < E.nseg ; s++ ) smooth_base_clear (&E.seg[s]);
    mpz_clears (A, x0, jj, mm, Seff, gtmp, rtmp, L, Hass, p, NULL);
    cornacchia_clear (&cc);  fp_clear (&C);
    return done ? 0 : 2;
}
