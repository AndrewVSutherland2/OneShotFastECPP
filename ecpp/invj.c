#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <gmp.h>
#include "invj.h"

static void bp_push (bipoly *P, int fdeg, int jdeg, const mpz_t c)
{
    if ( P->n == P->cap ) { P->cap = P->cap ? 2*P->cap : 256; P->t = realloc (P->t, P->cap*sizeof(bterm)); }
    bterm *b = &P->t[P->n++];
    b->fdeg = fdeg;  b->jdeg = jdeg;  mpz_init_set (b->c, c);
    if ( fdeg > P->maxf ) P->maxf = fdeg;
    if ( jdeg > P->maxj ) P->maxj = jdeg;
}

// Parse a bivariate polynomial in X (invariant) and Y (=j).  Terms are separated
// by top-level +/-; each term is [coeff][*X[^a]][*Y[^b]] (coeff/exponents optional).
int invj_load (bipoly *P, const char *phidir, const char *invstr)
{
    char path[4096];
    snprintf (path, sizeof path, "%s/phi_%s_j.txt", phidir, invstr);
    FILE *f = fopen (path, "r");
    if ( ! f ) return 0;
    // slurp + strip whitespace
    fseek (f, 0, SEEK_END);  long sz = ftell (f);  fseek (f, 0, SEEK_SET);
    char *raw = malloc (sz + 1);  size_t got = fread (raw, 1, sz, f);  raw[got] = 0;  fclose (f);
    char *b = malloc (got + 1);  size_t bl = 0;
    for ( size_t i = 0 ; i < got ; i++ ) if ( ! isspace ((unsigned char) raw[i]) ) b[bl++] = raw[i];
    b[bl] = 0;  free (raw);

    memset (P, 0, sizeof(*P));
    mpz_t c;  mpz_init (c);
    size_t i = 0;
    while ( i < bl ) {
        int sign = 1;
        if ( b[i] == '+' ) i++;  else if ( b[i] == '-' ) { sign = -1; i++; }
        int fdeg = 0, jdeg = 0, havec = 0;
        mpz_set_ui (c, 1);
        while ( i < bl && b[i] != '+' && b[i] != '-' ) {
            if ( b[i] == '*' ) { i++; continue; }
            if ( b[i] == 'X' || b[i] == 'Y' ) {
                char var = b[i++];  int e = 1;
                if ( i < bl && b[i] == '^' ) { i++;  e = 0;  while ( i < bl && isdigit ((unsigned char) b[i]) ) e = 10*e + (b[i++]-'0'); }
                if ( var == 'X' ) fdeg = e;  else jdeg = e;
            } else if ( isdigit ((unsigned char) b[i]) ) {
                size_t j = i;  while ( j < bl && isdigit ((unsigned char) b[j]) ) j++;
                char save = b[j];  b[j] = 0;  mpz_set_str (c, b + i, 10);  b[j] = save;  i = j;  havec = 1;
            } else i++;                                         // skip anything unexpected
        }
        if ( ! havec ) mpz_set_ui (c, 1);
        if ( sign < 0 ) mpz_neg (c, c);
        bp_push (P, fdeg, jdeg, c);
    }
    mpz_clear (c);  free (b);
    return P->n > 0;
}

void invj_clear (bipoly *P)
{ for ( int i = 0 ; i < P->n ; i++ ) mpz_clear (P->t[i].c);  free (P->t);  memset (P, 0, sizeof(*P)); }

int invj_load_phi (bipoly *P, const char *phidir, int ell)
{
    char path[4096];
    snprintf (path, sizeof path, "%s/phi_j_%d.txt", phidir, ell);
    FILE *f = fopen (path, "r");
    if ( ! f ) return 0;
    memset (P, 0, sizeof(*P));
    mpz_t c;  mpz_init (c);
    char *ln = NULL;  size_t lncap = 0;                        // Phi_97 coefficients run to ~750 digits
    while ( getline (&ln, &lncap, f) > 0 ) {
        int a, b;
        if ( sscanf (ln, "[%d,%d]", &a, &b) != 2 ) continue;
        char *q = strchr (ln, ']');  q++;
        while ( *q == ' ' || *q == '\t' ) q++;
        char *e = q + strlen (q);
        while ( e > q && isspace ((unsigned char) e[-1]) ) *--e = 0;
        if ( mpz_set_str (c, q, 10) != 0 ) continue;
        bp_push (P, a, b, c);
        if ( a != b ) bp_push (P, b, a, c);                    // symmetric storage: a >= b only
    }
    free (ln);  fclose (f);  mpz_clear (c);
    return P->n > 0;
}

int invj_jroots (const fp_ctx *C, const bipoly *P, const mp_limb_t *f0,
                 mp_limb_t *jroots, int maxj, uint64_t seed)
{
    int s = C->s;
    // powers of f0 (Montgomery): pw[k] = f0^k
    mp_limb_t *pw = malloc ((size_t)(P->maxf + 1) * s * sizeof(mp_limb_t));
    mpn_copyi (pw, C->R1, s);                                  // f0^0 = 1
    for ( int k = 1 ; k <= P->maxf ; k++ ) fp_mul (C, pw + (size_t)k*s, pw + (size_t)(k-1)*s, f0);

    // g(Y) = sum_terms  (c mod p) * f0^fdeg * Y^jdeg
    fp_poly g;  fpoly_init (C, &g, P->maxj + 2);
    for ( int k = 0 ; k <= P->maxj ; k++ ) memset (g.c + (size_t)k*s, 0, (size_t)s*sizeof(mp_limb_t));
    mpz_t cr;  mpz_init (cr);
    mp_limb_t cm[64], tm[64];
    for ( int i = 0 ; i < P->n ; i++ ) {
        mpz_mod (cr, P->t[i].c, C->pz);                        // coeff in [0,p)
        fp_set_mpz (C, cm, cr);                                // -> Montgomery
        fp_mul (C, tm, cm, pw + (size_t) P->t[i].fdeg * s);    // c * f0^fdeg
        fp_add (C, g.c + (size_t) P->t[i].jdeg * s, g.c + (size_t) P->t[i].jdeg * s, tm);
    }
    g.deg = P->maxj;
    while ( g.deg >= 0 && fp_is_zero (C, g.c + (size_t)g.deg*s) ) g.deg--;
    int nr = fp_find_all_roots (C, &g, jroots, maxj, seed);
    mpz_clear (cr);  fpoly_clear (&g);  free (pw);
    return nr;
}
