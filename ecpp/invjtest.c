#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "fproot.h"
#include "invj.h"

/*
    invjtest p D : classpoly picks the best invariant for D and computes
    H_D^inv mod p; we find a root f0 (fproot) and map it to j-invariants with the
    big-F_p invariant->j map (invj).  Prints "INV <str> <code>  F0 <f0>  J <j..>".
    A companion Python check cross-verifies against classpoly's word-size invtoj
    and against the Hilbert class polynomial.  Needs CLASSPOLY_PHI_DIR set.
*/

int main (int argc, char **argv)
{
    if ( argc < 3 ) { fprintf (stderr, "usage: invjtest p D\n"); return 2; }
    mpz_t p;  mpz_init_set_str (p, argv[1], 10);
    long D = strtol (argv[2], 0, 10);
    const char *phidir = getenv ("CLASSPOLY_PHI_DIR");
    if ( ! phidir ) { fprintf (stderr, "set CLASSPOLY_PHI_DIR\n"); return 1; }

    // 1. best-invariant H_D^inv mod p; capture the invariant string + code
    char pdec[4096];  gmp_snprintf (pdec, sizeof pdec, "%Zd", p);
    char outf[256];   snprintf (outf, sizeof outf, "/tmp/invjt_%ld.txt", -D);
    char cmd[8192];
    snprintf (cmd, sizeof cmd, "classpoly %ld -1 %s %s 1 2>&1 | grep -oE 'Using invariant [^ ]+ \\([0-9]+\\)'", D, pdec, outf);
    FILE *pp = popen (cmd, "r");  char line[256], invstr[64] = "";  int invcode = 0;
    if ( pp && fgets (line, sizeof line, pp) ) sscanf (line, "Using invariant %63s (%d)", invstr, &invcode);
    if ( pp ) pclose (pp);
    // strip trailing "(code)" token off invstr if grabbed
    char *paren = strchr (invstr, '(');  if ( paren ) *paren = 0;
    if ( ! invstr[0] ) { fprintf (stderr, "could not determine invariant\n"); return 1; }

    // 2. parse H_D^inv coefficients (low-degree first, monic)
    FILE *f = fopen (outf, "r");  if ( ! f ) { fprintf (stderr, "no %s\n", outf); return 1; }
    int cap = 64, deg = -1;  mpz_t *hc = malloc (cap*sizeof(mpz_t));  char ln[8192];
    while ( fgets (ln, sizeof ln, f) ) {
        char *x = strstr (ln, "*X^");  if ( ! x ) continue;
        int e = atoi (x + 3);
        if ( e >= cap ) { while ( e >= cap ) cap *= 2;  hc = realloc (hc, cap*sizeof(mpz_t)); }
        while ( deg < e ) { deg++; mpz_init (hc[deg]); }
        *x = 0;  mpz_set_str (hc[e], ln, 10);
    }
    fclose (f);
    if ( deg < 1 ) { fprintf (stderr, "degenerate H_D^inv\n"); return 1; }

    // 3. root f0 of H_D^inv
    fp_ctx C;  fp_init (&C, p);
    fp_poly h;  fpoly_init (&C, &h, deg + 2);
    for ( int i = 0 ; i <= deg ; i++ ) { mpz_mod (hc[i], hc[i], p);  fp_set_mpz (&C, h.c + (size_t)i*C.s, hc[i]); }
    h.deg = deg;
    mp_limb_t f0[64];  fp_find_root (&C, &h, f0, 0xABCDEF + (unsigned) D);
    mpz_t f0z;  mpz_init (f0z);  fp_get_mpz (&C, f0z, f0);

    // 4. map f0 -> j
    bipoly Phi;
    if ( ! invj_load (&Phi, phidir, invstr) ) { fprintf (stderr, "cannot load phi_%s_j.txt\n", invstr); return 1; }
    mp_limb_t jr[8*64];  int nj = invj_jroots (&C, &Phi, f0, jr, 8, 0x13579 + (unsigned) D);

    printf ("INV %s %d  DEG %d  F0 ", invstr, invcode, deg);
    gmp_printf ("%Zd  J", f0z);
    mpz_t jz;  mpz_init (jz);
    for ( int i = 0 ; i < nj ; i++ ) { fp_get_mpz (&C, jz, jr + (size_t)i*C.s);  gmp_printf (" %Zd", jz); }
    printf ("\n");

    mpz_clear (jz);  mpz_clear (f0z);  invj_clear (&Phi);  fpoly_clear (&h);
    for ( int i = 0 ; i <= deg ; i++ ) mpz_clear (hc[i]);
    free (hc);
    fp_clear (&C);  mpz_clear (p);
    return 0;
}
