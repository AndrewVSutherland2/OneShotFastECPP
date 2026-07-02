#ifndef _INVJ_H_
#define _INVJ_H_

#include <gmp.h>
#include "fproot.h"

/*
    Invariant -> j over a big prime field F_p (the big-F_p analog of classpoly's
    word-size invtoj / ff_all_j_from_inv).  classpoly picks the best class
    invariant for D (inv=-1) and computes the small H_D^inv mod p; we find a root
    f0 (the invariant value) with fproot, then map it to a j-invariant by finding
    a root of the modular polynomial Phi_inv(f0, J) over F_p.  Phi_inv is bivariate
    of small degree in J (1 or 2 for the invariants classpoly uses), read from
    phi_<invstr>_j.txt in the phi-file directory.
*/

typedef struct { int fdeg, jdeg; mpz_t c; } bterm;
typedef struct { bterm *t; int n, cap, maxf, maxj; } bipoly;

// load Phi_inv(X=invariant, Y=j) from <phidir>/phi_<invstr>_j.txt.  Returns 1 ok.
int  invj_load  (bipoly *P, const char *phidir, const char *invstr);
void invj_clear (bipoly *P);

// load the classical modular polynomial Phi_ell(X,Y) from <phidir>/phi_j_<ell>.txt
// (lines "[a,b] c" with a >= b; off-diagonal terms stand for X^a Y^b + X^b Y^a).
// The bundle has every prime ell <= 97.  Returns 1 ok.
int  invj_load_phi (bipoly *P, const char *phidir, int ell);

// given f0 (a root of H_D^inv, Montgomery form), write every j in F_p with
// Phi_inv(f0,j)=0 into jroots[i*s] (Montgomery); returns the count (<= maxj).
int  invj_jroots (const fp_ctx *C, const bipoly *P, const mp_limb_t *f0,
                  mp_limb_t *jroots, int maxj, uint64_t seed);

#endif
