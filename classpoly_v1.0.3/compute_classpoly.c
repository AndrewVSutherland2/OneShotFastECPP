#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "ff_poly.h"
#include "ff_poly/ffpolymultipoint.h"
#include "phi_poly.h"
#include "class_inv.h"
#include "classpoly.h"
#include "qform.h"
#include "pickprimes.h"
#include "classpoly_crt.h"
#include "classpoly_inv.h"
#include "evec.h"
#include "mpzutil.h"
#include "cstd.h"

/*
    Copyright 2012 Andrew V. Sutherland

    This file is part of classpoly.

    classpoly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    classpoly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with classpoly.  If not, see <http://www.gnu.org/licenses/>.
*/

#define DELTA_PERCENT		0.05
#define HEIGHT_MARGIN		256
#define MAX_SKIP_COUNT	(HEIGHT_MARGIN/64)

/*
	Decomposition of the class polynomial through a subgroup of the class group: Algorithm 1 of "Accelerating the CM method"
	(Sutherland, LMS J. Comput. Math. 15 (2012), arXiv:1009.1082), whose notation we follow.

	Let G be a subgroup of cl(D) of order n with m = h/n cosets (h/2n when we compute the square root of the class polynomial).
	The roots of the class polynomial H fall into m G-orbits, and with P_i(X) = prod_{x in orbit i}(X-x) = sum_k theta_{ik}X^k
	and y_i = theta_{i,n-1} (minus the sum of the orbit) we compute, in place of H,

		V(Y)   = prod_i (Y-y_i)					degree m
		W_k(Y) = sum_i theta_{ik} V(Y)/(Y-y_i)		degree < m, for 0 <= k <= n-2

	which have integer coefficients.  For a root y of V mod P that is not a root of V', the polynomial

		U(X) = X^n + y X^{n-1} + sum_{k<=n-2} (W_k(y)/V'(y)) X^k

	is one of the P_i mod P, so a root of H mod P is obtained from two root-findings of degrees m and n rather than one of degree mn.
	The coefficients of V and the W_k are also smaller than those of H (section 5 of the paper), so fewer CRT primes are needed.

	Which subgroups.  phi_enum_roots lists the roots in the order of the exponent vectors (e_0,...,e_{k-1}), 0 <= e_i < n[i], of the
	polycyclic presentation a_0,...,a_{k-1}, but the root in position e is [prod a_i^{s_i e_i}]x_0 for signs s_i = +/-1 that it does
	not control (one per generator: the gcd steps keep each generator's direction consistent across the whole enumeration, which is
	also what torsor_setup relies on).  So we need subgroups G whose cosets can be read off from positions whatever the signs are.
	Two families of such subgroups are used; write q_i for the image of a_i in Q = cl(D)/G.

	(1) Boxes.  G = < a_i^{c_i} : c_i < n[i] > with c_i | n[i], provided it has prod c_i cosets, which holds if and only if for
	every i with c_i < n[i] the power relation a_i^{n[i]} = prod_{j<i} a_j^{r_ij} has c_j | r_ij for all j < i.  Then c_i is the
	relative order of q_i in Q, and c_i q_i = 0 whenever a digit can exceed it, so the image of position e is
	sum s_i (e_i mod c_i) q_i; and for any signs b -> sum s_i b_i q_i is a bijection from the residue vectors to Q (peel off the last
	generator: Q modulo the earlier ones is cyclic of order c_k, on which negation is a bijection).  So the cosets are the sets of
	positions with fixed residues e_i mod c_i.  The subgroups <a_0,...,a_{d-1},a_d^e> of the paper's implementation are of this kind.

	(2) Lattices.  Let L_G be the lattice of all integer vectors x with prod a_i^{x_i} in G; it contains the relation lattice L_0 of
	the presentation, spanned by the rows n[i]e_i - sum_{j<i} r_ij e_j.  Positions e and e' are in the same coset if and only if
	(s_i(e_i-e'_i))_i lies in L_G, which certainly does not depend on the signs when L_G is mapped to itself by every change of sign
	of a coordinate.  Such lattices are easy to describe: with d_i the gcd of the i-th coordinates of L_G, the vectors 2d_i e_i lie
	in it (x minus x with its i-th sign changed), so

		L_G = { x : x_i = d_i y_i,  y mod 2 in S }     for a subspace S of F_2^k,

	and conversely every pair (d,S) for which this lattice contains L_0 is the L_G of a subgroup G with 2^(k-dim S) prod d_i cosets;
	the coset of a position is given by its residues e_i mod d_i and the class of ((e_i div d_i) mod 2)_i modulo S.  Containing L_0
	means: d_j divides n[j] and every r_ij, and S contains the rows of L_0 divided by d and reduced mod 2 (call their span S_0(d)).
	We enumerate all such d and, for each, the subspaces S_0(d) + <e_i : i in I> over subsets I of the generators.

	Neither family contains the other.  With relative orders (2,2,1061,2) and a_2^1061 = a_1 the subgroup of order 1061 generated
	by a_2^2 is the lattice d = (1,1,1,2), S = S_0(d), and no box; with relative orders (2,128,3,2) and a_2^3 = a_0 a_1^53 the box
	c = (1,4,3,2) is not a lattice of the second kind, because 4 does not divide 53 (a relation of a generator whose digit never
	exceeds its relative order in Q is not needed for (1)).  Both are stored as (d, parity checks), a box being d = c with no checks.

	We fix the symmetric function s = -e_1 (y_i is the coefficient of X^{n-1}).  That choice can fail to generate the fixed field
	of G, in which case V has a repeated root; this is detected below and the computation reverts to the plain class polynomial.
	Whoever consumes the output must check that the root of V they use is not a root of V'.

	The decomposition is requested through the environment, so that existing callers are unaffected:
		CLASSPOLY_TOWER=1		decompose, choosing G to minimize the estimated time for CRT plus root-finding mod P
							(the height bound alone when computing over Z)
		CLASSPOLY_TOWER=2		the same, and for an invariant whose sign is only fixed by the trace of the class polynomial
							write a sign-free decomposition (see classpoly_crt.c), which skips the computation of the trace
		CLASSPOLY_TOWER_N=n		decompose using a subgroup of order n (the best one of the form above with that order)
		CLASSPOLY_TOWER_C=c0,c1,..	decompose using exactly the box c (one entry per enumerated generator; for testing)
		CLASSPOLY_TOWER_C=d0,d1,../I	decompose using exactly the lattice (d, S_0(d) + <e_i : bit i of I>), as the candidate list prints it
		CLASSPOLY_TOWER_SLACK=x	add x percent to the height bound of the decomposition (default TOWER_SLACK_PERCENT, 0 for j)

	The height bound.  For j it is (5.10) of the paper and rigorous.  For other invariants the bounds on |j| are divided by the
	height factor, which is a heuristic: it is right on average over the class group, not coset by coset, and the true height has
	been seen above it.  A decomposition computed with too few primes is simply wrong, and nothing here can tell; the consumer
	must check its result (V and U split completely, the curve has the right trace) and be ready to ask again with more slack.
	The output is in the F/G format described in classpoly_crt.c (F = V, g_k = W_k).
*/
#define TOWER_MIN_DEGREE		24			// below this degree decomposing is not worth it (and tiny discriminants are where s=e_1 is known to fail)
#define TOWER_DUP_PRIMES		3			// declare V non-squarefree if the y_i collide modulo this many consecutive CRT primes
#define TOWER_SLACK_PERCENT		5.0			// added to the height bound of a decomposition for invariants other than j, where the bound is
								// heuristic (CLASSPOLY_TOWER_SLACK overrides).  Scaling the bounds on |j| by the height factor
								// is right on average over the class group but not coset by coset: measured heights run from
								// well below the unscaled bound to 3.0% above it (D=-18960475695, w_{5,7} ramified, 16 cosets).
#define TOWER_TRACE_FACTOR		1.6			// sign-free decompositions: the largest root is at most this power of what the height factor suggests
								// (1.5 for the Weber functions, whose pole at infinity has order 1/48 against a height factor of 72;
								// at most 1.44 for the eta quotients among the invariants concerned)
#define TOWER_MAX_CANDIDATES	16384			// subgroups considered
#define TOWER_MAX_EXT_BITS		12			// subsets I are taken from the first this many generators
#define TOWER_BOX				0xFFFFFFFFU		// in place of I: the candidate is a box
#define TOWER_BOUND_WORK		200000000L		// forms visited while computing height bounds; beyond this only the most promising subgroups are evaluated
// Cost model used to choose the subgroup, calibrated in November 2026 on a Zen 5 core (see reports/shortecpp-improve):
#define TOWER_CRT_NSECS			50.0			// CRT nanoseconds per bit of height bound per root of the class polynomial (35-75 measured)
#define TOWER_ROOT_NSECS		14.0			// root-finding nanoseconds per unit of degree per (bit of P)^2, on one core
#define TOWER_ROOT_PAR_DEGREE	4096			// degree above which fproot uses several cores, gaining a factor 1 + TOWER_ROOT_PAR_GAIN*lg(degree/4096)
#define TOWER_ROOT_PAR_GAIN		0.45			// in wall time (measured on 192 cores: 1.1 at degree 4682, 1.6 at 11364, 2.5 at 37456, 3.1 at 149392)
#define TOWER_ROOT_MULTI_DEGREE	48			// between this degree and TOWER_ROOT_PAR_DEGREE fproot takes its first step with several shifts
#define TOWER_ROOT_MULTI_GAIN	1.7			// in parallel, which gains this factor in wall time (measured: 12.6 -> 7.3 s, 34.6 -> 20.2 s)

struct tower_struct {
	long d[IQ_MAX_GENS];					// the cosets of G: residues e_i mod d[i], d0 <= i < k (indices as in pres->n, so the ell0 digit is index 0
	uint32_t hm[IQ_MAX_GENS];				// when present), and these nh parity checks on (e_i div d[i]) mod 2 (bit i-d0 is generator i)
	int nh;
	int d0, k;
	long n, m;							// order of G and number of cosets
	double bits;							// height bound for V and the W_k in bits, without the safety margin
	int *perm;							// perm[coset*n+j] = position in the enumeration of the j-th root of the coset
	ffpi_ctx_t ctx;						// scratch for the linear combinations
	ff_t *theta, *y, *col, *orbit;			// theta holds m rows of n+1 coefficients
};

struct tower_candidate { long d[IQ_MAX_GENS];  uint32_t hm[IQ_MAX_GENS];  uint32_t ext;  int nh;  long n, m;  double rank; };

// Reduced row echelon form over F_2, in place, pivots at the lowest set bits and rows sorted by pivot (so it is canonical for the
// span); returns the dimension.
static int tower_f2_rref (uint32_t *v, int nv)
{
	register uint32_t x, piv;
	register int i, j, r;

	for ( r = i = 0 ; i < nv ; i++ ) {
		x = v[i];
		for ( j = 0 ; j < r ; j++ ) if ( x & (v[j] & -v[j]) ) x ^= v[j];
		if ( ! x ) continue;
		piv = x & -x;
		for ( j = 0 ; j < r ; j++ ) if ( v[j] & piv ) v[j] ^= x;
		v[r++] = x;
	}
	for ( i = 1 ; i < r ; i++ ) { x = v[i];  for ( j = i ; j > 0 && (v[j-1] & -v[j-1]) > (x & -x) ; j-- ) v[j] = v[j-1];  v[j] = x; }
	return r;
}

// The subspace S = S_0(d) + <e_i : bit i of ext> for a vector d whose entries divide the columns of the relation lattice: writes its
// echelon basis to basis (room for 2k entries) and parity checks for it to hm, and returns the number of checks, k' - dim S.
static int tower_lattice (uint32_t *basis, uint32_t *hm, int *pdim, long *d, uint32_t ext, classgroup_pcp_t pres, int d0, int k)
{
	register uint32_t x;
	register int i, j, nv, r, nh;

	for ( nv = 0, i = d0 ; i < k ; i++ ) {
		x = (uint32_t)((pres->n[i]/d[i])&1) << (i-d0);
		for ( j = d0 ; j < i ; j++ ) if ( (evec_ri (pres->r, i)[j]/d[j]) & 1 ) x |= (uint32_t)1 << (j-d0);
		basis[nv++] = x;
	}
	for ( i = 0 ; i < k-d0 ; i++ ) if ( ext & ((uint32_t)1<<i) ) basis[nv++] = (uint32_t)1 << i;
	r = tower_f2_rref (basis, nv);
	// one check per column f that is not a pivot: bit f, and the pivot bit of every basis vector that has bit f set
	for ( nh = 0, i = 0 ; i < k-d0 ; i++ ) {
		for ( j = 0 ; j < r ; j++ ) if ( (basis[j] & -basis[j]) == ((uint32_t)1<<i) ) break;
		if ( j < r ) continue;
		x = (uint32_t)1 << i;
		for ( j = 0 ; j < r ; j++ ) if ( basis[j] & ((uint32_t)1<<i) ) x |= basis[j] & -basis[j];
		hm[nh++] = x;
	}
	*pdim = r;
	return nh;
}

static double tower_root_secs (double deg, double Pbits, int jobs)
{
	double t = TOWER_ROOT_NSECS * 1e-9 * deg * Pbits * Pbits * (1.0 + log2 (deg+1.0)/16.0);
	if ( jobs <= 1 ) return t;
	if ( deg > TOWER_ROOT_PAR_DEGREE ) return t / (1.0 + TOWER_ROOT_PAR_GAIN * log2 (deg/TOWER_ROOT_PAR_DEGREE));
	return ( deg >= TOWER_ROOT_MULTI_DEGREE ? t / TOWER_ROOT_MULTI_GAIN : t );
}

static int _tower_cmp_candidates (const void *a, const void *b)
{
	const struct tower_candidate *x = a, *y = b;
	if ( x->rank != y->rank ) return ( x->rank < y->rank ? -1 : 1 );
	return ( x->m < y->m ? -1 : (x->m > y->m ? 1 : 0) );
}

/*
	Chooses the subgroup.  Returns 1 and fills in T (except for the scratch space) if a decomposition should be used, 0 otherwise.
	full_bits is the height bound for the plain class polynomial; hf the height factor of the invariant.  The qform table for D
	must still be loaded.  The choice depends only on D, the presentation, P and the environment, so the workers of a multi-job
	computation all make the same one.
*/
static int tower_select (struct tower_struct *T, classgroup_pcp_t pres, int H_d, double hf, double full_bits, mpz_t P)
{
	struct tower_candidate *cand, *best;
	char *s;
	uint32_t basis[2*IQ_MAX_GENS], hm[IQ_MAX_GENS], *seen, ext, next, full, proj;
	long d[IQ_MAX_GENS], g[IQ_MAX_GENS], fc[IQ_MAX_GENS], force, sub, cosets, maxeval, hp, prodd, fext, r;
	double summax, maxdelta, bits, Pbits, cost, best_cost, crt_secs;
	int d0, k, kk, i, j, jobs, auto_on, ncand, nfc, nh, dim, nseen, lat;

	s = getenv ("CLASSPOLY_TOWER_N");  force = ( s ? atol (s) : 0 );
	s = getenv ("CLASSPOLY_TOWER");  auto_on = ( s && atoi (s) > 0 );
	nfc = lat = 0;  fext = 0;
	if ( (s = getenv ("CLASSPOLY_TOWER_C")) && *s ) {
		while ( nfc < IQ_MAX_GENS ) { fc[nfc++] = strtol (s, &s, 10);  if ( *s != ',' ) break;  s++; }
		if ( *s == '/' ) { lat = 1;  fext = strtol (s+1, 0, 10); }
		force = -1;									// the order is whatever this gives
	}
	if ( ! auto_on && ! force ) return 0;
	k = pres->k;  d0 = ( pres->ell0 ? 1 : 0 );  kk = k-d0;
	if ( k <= d0 ) return 0;
	if ( nfc && nfc != kk ) { err_printf ("CLASSPOLY_TOWER_C has %d entries, the presentation for D=%ld has %d generators\n", nfc, pres->D, kk);  return 0; }
	if ( ! force && H_d < TOWER_MIN_DEGREE ) return 0;
	s = getenv ("CLASSPOLY_JOBS");  jobs = ( s && atoi (s) > 1 ? atoi (s) : 1 );
	Pbits = ( P ? (double) mpz_sizeinbase (P, 2) : 0.0 );
	full = ( kk >= 32 ? ~(uint32_t)0 : ((uint32_t)1 << kk) - 1 );
	for ( hp = 1, i = d0 ; i < k ; i++ ) hp *= pres->n[i];
	// g[j] = gcd of the j-th column of the relation lattice: what d[j] has to divide
	for ( j = d0 ; j < k ; j++ ) {
		g[j] = pres->n[j];
		for ( i = j+1 ; i < k ; i++ ) { r = evec_ri (pres->r, i)[j];  if ( r < 0 ) r = -r;  if ( r ) g[j] = (long) ui_gcd (g[j], r); }
	}

	cand = mem_alloc (TOWER_MAX_CANDIDATES * sizeof(*cand));
	ncand = 0;
	if ( nfc ) {
		// exactly one subgroup: a box c, or a lattice (d,I)
		for ( i = d0 ; i < k ; i++ ) {
			d[i] = fc[i-d0];
			if ( d[i] <= 0 ) { err_printf ("CLASSPOLY_TOWER_C: entries must be positive\n");  mem_free (cand);  return 0; }
			if ( ( lat ? g[i] : pres->n[i] ) % d[i] ) break;
			if ( ! lat && d[i] < pres->n[i] ) { for ( j = d0 ; j < i ; j++ ) if ( evec_ri (pres->r, i)[j] % d[j] ) break;  if ( j < i ) break; }
		}
		if ( i == k ) {
			ext = ( lat ? (uint32_t) fext & full : TOWER_BOX );
			nh = ( lat ? tower_lattice (basis, hm, &dim, d, ext, pres, d0, k) : 0 );
			for ( prodd = 1, i = d0 ; i < k ; i++ ) prodd *= d[i];
			cosets = prodd << nh;  sub = hp / cosets;
			if ( sub > 1 && cosets > 1 ) {
				for ( i = d0 ; i < k ; i++ ) cand[0].d[i] = d[i];
				for ( i = 0 ; i < nh ; i++ ) cand[0].hm[i] = hm[i];
				cand[0].nh = nh;  cand[0].ext = ext;  cand[0].n = sub;  cand[0].m = cosets;  cand[0].rank = 0.0;
				ncand = 1;
			}
		}
	} else {
		// (1) the boxes: every c with c_i | n[i], like an odometer, keeping those with prod c_i cosets
		for ( i = d0 ; i < k ; i++ ) d[i] = 1;
		for (;;) {
			for ( dim = 1, i = d0 ; dim && i < k ; i++ )
				if ( d[i] < pres->n[i] ) { for ( j = d0 ; j < i ; j++ ) if ( evec_ri (pres->r, i)[j] % d[j] ) { dim = 0; break; } }
			for ( sub = cosets = 1, i = d0 ; i < k ; i++ ) { cosets *= d[i];  sub *= pres->n[i]/d[i]; }
			if ( dim && sub > 1 && cosets > 1 && (force <= 0 || sub == force) && ncand < TOWER_MAX_CANDIDATES ) {
				for ( i = d0 ; i < k ; i++ ) cand[ncand].d[i] = d[i];
				cand[ncand].nh = 0;  cand[ncand].ext = TOWER_BOX;  cand[ncand].n = sub;  cand[ncand].m = cosets;
				// rank by what we know before computing a height bound: root-finding time mod P, or (over Z) the number of cosets
				cand[ncand].rank = ( Pbits > 0 ? tower_root_secs (cosets, Pbits, jobs) + tower_root_secs (sub, Pbits, jobs) : (double) cosets );
				ncand++;
			}
			for ( i = d0 ; i < k ; i++ ) {		// next divisor of n[i], carrying into the next digit when they are exhausted
				do { d[i]++; } while ( d[i] <= pres->n[i] && (pres->n[i] % d[i]) );
				if ( d[i] <= pres->n[i] ) break;
				d[i] = 1;
			}
			if ( i == k ) break;
		}
		// (2) the lattices that are not boxes
		// every d with d[j] | g[j], like an odometer, and for each the subspaces S_0(d) + <e_i : i in I>
		seen = mem_alloc ((((size_t)1 << ( kk < TOWER_MAX_EXT_BITS ? kk : TOWER_MAX_EXT_BITS )) + 1) * IQ_MAX_GENS * sizeof(*seen));
		next = (uint32_t)1 << ( kk < TOWER_MAX_EXT_BITS ? kk : TOWER_MAX_EXT_BITS );
		for ( i = d0 ; i < k ; i++ ) d[i] = 1;
		for (;;) {
			for ( prodd = 1, i = d0 ; i < k ; i++ ) prodd *= d[i];
			for ( nseen = 0, ext = 0 ; ext < next ; ext++ ) {
				nh = tower_lattice (basis, hm, &dim, d, ext, pres, d0, k);
				for ( proj = 0, i = 0 ; i < dim ; i++ ) proj |= basis[i];
				if ( proj != full ) continue;				// some y_i is always even: the same lattice with 2d_i, found there
				for ( i = 0 ; i < nseen ; i++ ) {			// the same subspace from another I?
					if ( seen[i*IQ_MAX_GENS] != (uint32_t) dim ) continue;
					for ( j = 0 ; j < dim ; j++ ) if ( seen[i*IQ_MAX_GENS+1+j] != basis[j] ) break;
					if ( j == dim ) break;
				}
				if ( i < nseen ) continue;
				if ( dim < IQ_MAX_GENS ) { seen[nseen*IQ_MAX_GENS] = (uint32_t) dim;  for ( j = 0 ; j < dim ; j++ ) seen[nseen*IQ_MAX_GENS+1+j] = basis[j];  nseen++; }
				cosets = prodd << nh;  sub = hp / cosets;
				if ( nh && sub > 1 && cosets > 1 && (force <= 0 || sub == force) && ncand < TOWER_MAX_CANDIDATES ) {	// no checks: the box d
					for ( i = d0 ; i < k ; i++ ) cand[ncand].d[i] = d[i];
					for ( i = 0 ; i < nh ; i++ ) cand[ncand].hm[i] = hm[i];
					cand[ncand].nh = nh;  cand[ncand].ext = ext;  cand[ncand].n = sub;  cand[ncand].m = cosets;
					// rank by what we know before computing a height bound: root-finding time mod P, or (over Z) the number of cosets
					cand[ncand].rank = ( Pbits > 0 ? tower_root_secs (cosets, Pbits, jobs) + tower_root_secs (sub, Pbits, jobs) : (double) cosets );
					ncand++;
				}
			}
			for ( i = d0 ; i < k ; i++ ) {		// next divisor of g[i], carrying into the next digit when they are exhausted
				do { d[i]++; } while ( d[i] <= g[i] && (g[i] % d[i]) );
				if ( d[i] <= g[i] ) break;
				d[i] = 1;
			}
			if ( i == k ) break;
		}
		mem_free (seen);
	}
	if ( ! ncand ) {
		if ( force ) { err_printf ("The requested subgroup (CLASSPOLY_TOWER_N or CLASSPOLY_TOWER_C) is not usable for D=%ld\n", pres->D); }
		else { info_printf ("No proper subgroup available for D=%ld, not decomposing\n", pres->D); }
		mem_free (cand);
		return 0;
	}
	qsort (cand, ncand, sizeof(*cand), _tower_cmp_candidates);
	maxeval = TOWER_BOUND_WORK / pres->h;  if ( maxeval < 8 ) maxeval = 8;
	if ( ncand > maxeval ) { info_printf ("Evaluating height bounds for the %ld most promising of %d subgroups\n", maxeval, ncand);  ncand = (int) maxeval; }

	best = 0;  best_cost = 0.0;
	info_printf ("Decomposition candidates for D=%ld (degree %d, plain height bound %.0f bits, P has %.0f bits, %d job%s):\n", pres->D, H_d, full_bits, Pbits, jobs, jobs > 1 ? "s" : "");
	for ( i = 0 ; i < ncand ; i++ ) {
		// the candidates are in the order of their root-finding time, which is part of the cost: the rest cannot win
		if ( Pbits > 0 && best && cand[i].rank >= best_cost ) { info_printf ("  (%d more with at least %.1fs of root-finding)\n", ncand-i, cand[i].rank);  break; }
		qform_sum_jbounds_lat (&summax, &maxdelta, pres->n, cand[i].d, cand[i].hm, cand[i].nh, k, d0, pres->D, hf);
		sub = cand[i].n;  cosets = cand[i].m;
		// (5.10) of the paper, with the j-invariant bounds scaled by the height factor inside qform_sum_jbounds_lat
		bits = MPZ_LOG2E * (summax + maxdelta) + log2 ((double)cosets) + cosets + sub + cosets * log2 ((double)sub);
		crt_secs = TOWER_CRT_NSECS * 1e-9 * bits * H_d / jobs;
		cost = ( Pbits > 0 ? crt_secs + cand[i].rank : bits );
		if ( dbg_level >= INFO_LEVEL ) {
			out_printf ("  n=%ld m=%ld c=(", sub, cosets);
			for ( j = d0 ; j < k ; j++ ) out_printf ("%s%ld", j > d0 ? "," : "", cand[i].d[j]);
			if ( cand[i].ext == TOWER_BOX ) { out_printf (")"); } else { out_printf (")/%u", cand[i].ext); }
			out_printf (": %.0f bits, est. %.1fs CRT + %.1fs + %.1fs roots\n", bits, crt_secs, tower_root_secs (cosets, Pbits, jobs), tower_root_secs (sub, Pbits, jobs));
		}
		if ( ! best || cost < best_cost ) { best = cand + i;  best_cost = cost;  T->bits = bits; }
	}
	if ( ! force && Pbits > 0 ) {
		// decompose only if the model says it is no slower than the plain computation
		cost = TOWER_CRT_NSECS * 1e-9 * full_bits * H_d / jobs + tower_root_secs (H_d, Pbits, jobs);
		if ( cost <= best_cost ) { info_printf ("Plain class polynomial estimated at %.1fs, decomposition at %.1fs: not decomposing\n", cost, best_cost);  mem_free (cand);  return 0; }
	}
	T->d0 = d0;  T->k = k;  T->n = best->n;  T->m = best->m;  T->nh = best->nh;
	for ( i = d0 ; i < k ; i++ ) T->d[i] = best->d[i];
	for ( i = 0 ; i < best->nh ; i++ ) T->hm[i] = best->hm[i];
	info_printf ("Decomposing with subgroup of order n=%ld, m=%ld cosets, height bound %.0f bits (plain: %.0f)\n", T->n, T->m, T->bits, full_bits);
	mem_free (cand);
	return 1;
}

// Allocates the scratch space and tabulates which positions of the enumeration hold the roots of each coset (qform_lat_coset
// numbers the cosets; within one the roots are taken in the order of the enumeration).
static void tower_alloc (struct tower_struct *T, classgroup_pcp_t pres)
{
	long e[IQ_MAX_GENS], idx, coset, *fill;
	int i;

	ffpi_alloc_ctx (T->ctx, (int)T->m);
	T->theta = mem_alloc ((T->m*(T->n+1) + 2*T->m + T->n + 4) * sizeof(ff_t));
	T->y = T->theta + T->m*(T->n+1);  T->col = T->y + T->m + 1;  T->orbit = T->col + T->m + 1;
	T->perm = mem_alloc (T->m*T->n*sizeof(int));
	fill = mem_alloc (T->m*sizeof(*fill));
	for ( coset = 0 ; coset < T->m ; coset++ ) fill[coset] = 0;
	for ( i = T->d0 ; i < T->k ; i++ ) e[i] = 0;
	for ( idx = 0 ; idx < T->m*T->n ; idx++ ) {
		coset = qform_lat_coset (e, T->d, T->hm, T->nh, T->k, T->d0);
		if ( coset >= T->m || fill[coset] >= T->n ) { err_printf ("Decomposition: the cosets of the chosen subgroup do not have %ld positions each (D=%ld)\n", T->n, pres->D);  exit (0); }
		T->perm[coset*T->n+fill[coset]++] = (int) idx;
		for ( i = T->d0 ; i < T->k && ++e[i] == pres->n[i] ; i++ ) e[i] = 0;
	}
	mem_free (fill);
}

static void tower_free (struct tower_struct *T)
	{ ffpi_free_ctx (T->ctx);  mem_free (T->theta);  mem_free (T->perm); }

static int _tower_cmp_ff (const void *a, const void *b)
	{ ff_t x = *(const ff_t *)a, y = *(const ff_t *)b;  return ( x < y ? -1 : (x > y ? 1 : 0) ); }

/*
	Given the roots of the class polynomial mod _ff_p in the order produced by phi_enum_roots, writes the m*n coefficients of the
	decomposition to out (which needs room for m*n+1 entries): V_0..V_{m-1}, then W_0, ..., W_{n-2}, each with m coefficients.
	Returns 0 if two of the y_i coincide mod _ff_p (the output is still correct in that case), 1 otherwise.
*/
static int tower_build (ff_t *out, ff_t *roots, struct tower_struct *T)
{
	register long c, i, k, n, m;
	register ff_t *th;
	register int *pm;
	int distinct;

	n = T->n;  m = T->m;
	for ( c = 0, pm = T->perm ; c < m ; c++, pm += n ) {
		for ( i = 0 ; i < n ; i++ ) _ff_set (T->orbit[i], roots[pm[i]]);
		th = T->theta + c*(n+1);
		ff_poly_from_roots (th, T->orbit, (int)n);
		_ff_set (T->y[c], th[n-1]);
	}
	memcpy (T->col, T->y, m*sizeof(ff_t));
	qsort (T->col, m, sizeof(ff_t), _tower_cmp_ff);
	for ( distinct = 1, c = 1 ; c < m ; c++ ) if ( _ff_equal (T->col[c], T->col[c-1]) ) { distinct = 0; break; }
	ff_poly_from_roots (out, T->y, (int)m);					// also writes the leading 1 at out[m], which W_0 then overwrites
	ffpi_setup_ctx (T->ctx, T->y, (int)m, 1);
	for ( k = 0 ; k < n-1 ; k++ ) {
		for ( c = 0 ; c < m ; c++ ) _ff_set (T->col[c], T->theta[c*(n+1)+k]);
		ffpi_combo (out + m + k*m, T->col, (int)m, T->ctx);
	}
	return distinct;
}

int compute_classpoly (long D, int inv, mpz_t P, char *filename)
{
	time_t begin, start, end, last;
	unsigned long cntr1,cntr2;
	unsigned long fc_cycles, er_cycles, pr_cycles, tot_fc_cycles, tot_er_cycles, tot_pr_cycles, tot_cycles;
	long cbits, outbits;
	char buf[256];
	FILE *fp;
	classpoly_t H;
	classpoly_inv_t Hinv;
	classpoly_crt_t crt;
	double bits, tbits, height_factor, lastbits, percent;
	long vfilter, maxp, maxv, last_index;
	long *crt_p;
	ff_t *work, J;
	int  crt_pcnt, p1, p2;
	int H_d, ell0, N, classno, skip_count;
	struct tower_struct tower[1];
	int tower_on, tower_seen, tower_dups, tower_signfree;
	ff_t tr, *tout;
	register int i;
	
	if ( P && ! mpz_sgn(P) ) P = 0;
	if ( ! filename) { filename = buf;  if ( inv < D) sprintf(filename, "H_%ld.txt", -D); else sprintf(filename, "H_%ld_%d.txt", -D, inv); }

	// Multi-job ECRT merge mode (CLASSPOLY_JOBS=W, CLASSPOLY_JOBID=0): combine the
	// partial coefficient sums dumped by the W workers and write the class
	// polynomial; no computation of our own is needed (n, k, P all come from the
	// dump files, whose consistency ecrt_merge verifies).
	{
		char *js = getenv("CLASSPOLY_JOBS"), *jis = getenv("CLASSPOLY_JOBID");
		int mjobs = js ? atoi(js) : 0;
		if ( mjobs > 1 && jis && atoi(jis) == 0 ) {
			ecrt_context_t ecrt;
			char prefix[160], jobfile[224], *dir;
			mpz_t X;
			if ( ! P ) { err_printf ("ECRT merge mode requires a nonzero modulus P\n"); return 0; }
			memset (ecrt, 0, sizeof(ecrt[0]));
			dir = getenv("CLASSPOLY_ECRT_DIR");
			snprintf (prefix, sizeof(prefix), "%s%sH_%ld", dir ? dir : "", dir ? "/" : "", -D);
			ecrt_merge (ecrt, mjobs, prefix);
			// the tag gives the shape of a decomposition: its number of cosets (0 if none) and whether it is sign-free
			int tm = classpoly_tower_tag_m (ecrt->tag), tsf = classpoly_tower_tag_signfree (ecrt->tag);
			if ( tm && ((ecrt->k - tsf) % tm) ) { err_printf ("ECRT dump files carry tag %lu, which does not fit the coefficient count %d\n", (unsigned long)ecrt->tag, ecrt->k); return 0; }
			fp = fopen (filename, "w");
			if ( ! fp ) { err_printf ("Error opening output file %s\n", filename); return 0; }
			classpoly_crt_print_header (fp, D, inv, ecrt->P, ecrt->k, tm, tsf);
			mpz_init (X);
			for ( i = 0 ; ecrt_next_coeff (X, ecrt) ; i++ ) classpoly_crt_print_coeff (fp, X, i, ecrt->k, tm, tsf);
			mpz_clear (X);
			fclose (fp);
			for ( i = 1 ; i <= mjobs ; i++ ) {	// job dumps are consumed; remove them
				snprintf (jobfile, sizeof(jobfile), "%s_%d_%d.ecrt", prefix, mjobs, i);
				if ( ecrt->infp[i-1] ) fclose (ecrt->infp[i-1]);
				remove (jobfile);
			}
			info_printf ("ECRT merge of %d jobs complete, wrote %d coefficients to %s\n", mjobs, i-1, filename);
			return 1;
		}
	}
	
	if ( D >= -4 ) {
		if ( inv || P ) { err_printf ("D must be less than -4 unless inv=0 and P=0\n"); return 0; }
		fp = fopen (filename, "w");
		if ( ! fp ) { err_printf("Unable to create file %s\n", filename); return 0; }
		fprintf (fp, "I=0\nD=%ld\n", D);
		fprintf (fp, "%ld*X^0\n", (D==-3?0L:-1728L));
		fprintf (fp, "1*X^1\n");
		fclose (fp);
		out_printf ("Class polynomial for D=%ld written to %s\n", D,  filename);
		return 1;
	}
	begin = clock();

	// *** Alg 2 Step 1 *** Compute a polycyclic presentation for cl(D)  ***
	phi_poly_setup ();
	classpoly_init (H, D, inv);
	classpoly_setup_find_jinv (H);

	// if f(z) is invariant under the Fricke involution and all primes dividing the degree (not always the same as the level!) are ramified, then
	// every root of the classpoly is a double root and we compute the sqrt of classpoly, see section 4.1 of "Class invariants by the CRT method"
	// in this case we want to setup the presentation so that we only need to enumerate each double-root once
	N = inv_level (inv);
	ell0 = ( (inv_atkin(inv) || inv_double_eta(inv)) && inv_ramified(D,inv) ) ? inv_degree(&p1,&p2,inv) : 0;
	if ( ell0 ) { info_printf ("Computing square-root of classpoly for D=%ld\n", D); }
	else { info_printf ("Computing classpoly for D=%ld\n", D); }

	if ( ! classpoly_setup_enum_roots (H, N, ell0, inv) ) return 0;
	if ( dbg_level >= INFO_LEVEL ) classgroup_pcp_print (H->pres);
	classno = H->pres->h;  H_d = (int) classno / (ell0?2:1);
	height_factor = inv_enum_height_factor(inv);
	if ( ! height_factor ) { fprintf (stderr, "height factor error for D=%ld\n", D); assert (height_factor); }
	bits = qform_hilbert_height_bound (D, H_d) / height_factor;
	tower_on = tower_select (tower, H->pres, H_d, height_factor, bits, P);
	tower_seen = tower_dups = tower_signfree = 0;
	if ( tower_on ) {
		char *ts = getenv ("CLASSPOLY_TOWER"), *sl = getenv ("CLASSPOLY_TOWER_SLACK");
		bits = tower->bits * (1.0 + ( sl ? atof (sl) : (inv == INV_J ? 0.0 : TOWER_SLACK_PERCENT) ) / 100.0);
		if ( ts && atoi (ts) >= 2 && classpoly_inv_needs_trace (inv, H->pres->enum_cnt) ) {
			// sign-free: the odd coefficients are multiplied by the trace, which is at most H_d times the largest root
			tower_signfree = 1;
			bits += TOWER_TRACE_FACTOR * MPZ_LOG2E * jbound (1, D) / height_factor + log2 ((double)H_d) + 32;
			info_printf ("Sign-free decomposition: no trace computation, height bound %.0f bits\n", bits);
		}
	}
	if ( inv_single_eta(inv) ) bits = 1.1*bits+200;						// heuristic height estimate for single-eta quotients is too low, so bump it up a notch
	bits += HEIGHT_MARGIN;
	info_printf ("bits = %.2f (classpoly degree %d height factor %.2f, safety margin %d bits)\n", bits, H_d, height_factor, HEIGHT_MARGIN);
	vfilter = classgroup_pcp_vfilter (H->pres);

	// *** Alg 2 Step 2 *** Pick primes (note crt_p is allocated by pick_primes and freed by classpoly_crt_finish or classpoly_crt_end)
	crt_pcnt = pick_primes (&crt_p, 0, 0, D, classno, bits, inv_pfilter(inv), vfilter);
	for ( i = 0, maxp = maxv = 0, tbits = 0.0 ; i < crt_pcnt ; i++ ) {
		if ( crt_p[i] > maxp ) maxp = crt_p[i];
		tbits += log2(crt_p[i]);
	}	
	info_printf ("Selected %d primes, maxp=%ld, totbits = %.1f\n", crt_pcnt, maxp, tbits);

	// *** Alg 2 Step 3 *** CRT setup -- no explicit CRT for the moment
	classpoly_crt_start (crt, D, inv, (unsigned long *)crt_p, crt_pcnt, H_d + tower_signfree, P);
	if ( tower_on ) { crt->tower_m = (int) tower->m;  crt->tower_signfree = tower_signfree; }
	if ( ! classpoly_inv_setup_signfree (Hinv, H, inv, crt, tower_signfree) ) { classpoly_clear (H); classpoly_crt_end (crt); return 0; }

	work = mem_alloc ((2*H_d+2)*sizeof(ff_t));
	if ( tower_on ) tower_alloc (tower, H->pres);
	
	qform_table_free();	
	end = clock();
	info_printf ("Total precomputation time %ld msecs, beginning main loop...\n", delta_msecs(begin,end));
	
	percent = DELTA_PERCENT;
	last_index = 0;  lastbits = tbits = 0.0;
	fc_cycles = er_cycles = pr_cycles = tot_fc_cycles = tot_er_cycles = tot_pr_cycles = 0;

	last = start = clock();
	// *** Alg 2 Step 4 ****
	skip_count = 0 ;
	while ( classpoly_crt_next_prime (crt) ) {
		// multi-job ECRT: this worker handles only its own residue class of prime indices
		if ( crt->jobs && (crt->index % crt->jobs) != crt->jobid - 1 ) continue;

		// *** Alg 1 Steps 1+2
		cntr1 = get_cycles();
		if ( ! classpoly_crt_cached_j (&J, crt) && ! ff_classpoly_find_jinv (&J, H) ) { err_printf ("Fatal error, ff_classpoly_find_jinv failed for p=%ld\n", _ff_p); abort(); }
		dbg_printf ("j-invariant %ld is a root of H_D mod %ld\n", _ff_get_ui(J), _ff_p);
		cntr2 = get_cycles();
		fc_cycles += (unsigned long)(cntr2-cntr1);  tot_fc_cycles +=  (unsigned long)(cntr2-cntr1); 
		
		// *** Alg 1 Step 3
		cntr1 = get_cycles();
		ff_classpoly_setup (H);
		ff_classpoly_inv_from_j (&J, &J, Hinv, H);
		dbg_printf("Enumerating %s invariants mod %ld starting at %ld\n", inv_string(inv), _ff_p, _ff_get_ui(J));
		if ( ! ff_classpoly_enum_roots (H, J, inv) ) { err_printf ("Call to ff_classpoly_enum_roots failed\n"); abort(); }
		cntr2 = get_cycles();
		er_cycles += (unsigned long)(cntr2-cntr1);  tot_er_cycles +=  (unsigned long)(cntr2-cntr1); 
		dbg_printf("Successfully enumerated %d roots of H_D mod %ld\n", H_d, _ff_p);
		if ( dbg_level >= DEBUG_LEVEL ) { dbg_printf ("enum (p=%ld):\n", _ff_p); for ( i = 0 ; i < H_d ; i++ ) dbg_printf ("    [%d]=%ld\n", i,_ff_get_ui(H->roots[i])); }

		// in cases where classpoly is ambiguous, use oriented enumeration to disambiguate as described in section 4.1 of "Class invariants by the CRT method", and apply the trace trick of section 4.2 if needed
		if ( ! ff_classpoly_inv_disambiguate (H, Hinv) ) {
			if ( P ) { err_printf ("Call to ff_classpoly_inv_disambiguate failed for p=%ld, prime skipping not currently supported with ECRT\n", _ff_p); abort(); }
			if ( skip_count > MAX_SKIP_COUNT ) { err_printf ("Call to ff_classpoly_inv_disambiguate failed for p=%ld\n, exceepd MAX_SKIP_COUNT=%d\n", _ff_p, MAX_SKIP_COUNT); abort(); }
			out_printf ("Call to ff_classpoly_inv_disambiguate failed, skipping prime %ld\n", _ff_p);
			classpoly_crt_skip_prime (crt);
			continue;
		}

		// *** Alg 1 Step 4
		cntr1 = get_cycles();
		if ( tower_on ) {
			// steps 5c and 5d of Algorithm 1 in arXiv:1009.1082.  If the y_i collide modulo each of the first few primes we see, then
			// (up to a negligible chance) V has a repeated root over Z, s=e_1 is unusable for this subgroup, and we start over without it.
			tout = work + tower_signfree;
			if ( ! tower_build (tout, H->roots, tower) && tower_seen < TOWER_DUP_PRIMES ) tower_dups++;
			if ( ++tower_seen == TOWER_DUP_PRIMES && tower_dups == TOWER_DUP_PRIMES ) {
				out_printf ("Decomposition of the class polynomial for D=%ld, inv=%d with n=%ld, m=%ld has a repeated coset value, reverting to the plain class polynomial\n", D, inv, tower->n, tower->m);
				mem_free (work);  tower_free (tower);  classpoly_clear (H);  classpoly_inv_clear (Hinv);  classpoly_crt_end (crt);
				setenv ("CLASSPOLY_TOWER", "0", 1);  unsetenv ("CLASSPOLY_TOWER_N");  unsetenv ("CLASSPOLY_TOWER_C");
				return compute_classpoly (D, inv, P, filename);
			}
			if ( tower_signfree ) {
				// the sign of the roots has not been fixed: store S = T^2 and T times each coefficient that changes sign with the roots
				_ff_set_zero (tr);  for ( i = 0 ; i < H_d ; i++ ) _ff_addto (tr, H->roots[i]);
				for ( i = 0 ; i < H_d ; i++ ) if ( classpoly_tower_odd (i, (int)tower->m, (int)tower->n) ) ff_mult (tout[i], tout[i], tr);
				_ff_square (work[0], tr);
			}
		} else {
			ff_poly_from_roots (work, H->roots, H_d);
		}
		cntr2 = get_cycles();
		pr_cycles += (unsigned long)(cntr2-cntr1);  tot_pr_cycles +=  (unsigned long)(cntr2-cntr1); 
		if ( dbg_level >= DEBUG_LEVEL ) { dbg_printf ("poly (p=%ld)  x^%d",_ff_p, H_d); for ( i = H_d-1 ; i >= 0 ; i-- ) dbg_printf ("+ %ld*x^%d ", _ff_get_ui(work[i]), i); dbg_printf("\n"); }
		
		// *** Alg 2 Step 4.b CRT update
		for ( i = 0 ; i < H_d + tower_signfree ; i++ ) work[i] = _ff_get_ui(work[i]);		// convert out of Montgomery rep
		classpoly_crt_update (crt, work);
		
		// Report progress 
		tbits += log2(_ff_p);
		if ( dbg_level >= INFO_LEVEL && tbits/bits > percent ) {
			end = clock();
			tot_cycles = fc_cycles+er_cycles+pr_cycles;
			out_printf ("%12ld  (%.3f,%.3f,%.3f)  %6.1f msec/prime    %7.1f bits/s %8ld msecs (last) %10ld msecs (tot)  [%.2f]\n",
				   _ff_p, (double)fc_cycles/tot_cycles, (double)er_cycles/tot_cycles,  (double)pr_cycles/tot_cycles, (double)delta_msecs(last,end)/(crt->index-last_index),
				   1000.0*(tbits-lastbits)/delta_msecs(last,end), delta_msecs(last,end), delta_msecs(start,end), tbits/bits);
			percent += DELTA_PERCENT;  last = end;  last_index = crt->index; lastbits = tbits;  fc_cycles= er_cycles = pr_cycles = 0;
		}
	}
	mem_free (work);
	if ( tower_on ) tower_free (tower);
	classpoly_clear (H);
	classpoly_inv_clear (Hinv);
	end = clock();
	info_printf ("Completed main loop in %ld msecs\n", delta_msecs(start,end));

	if ( crt->jobs ) {	// multi-job worker: dump partial sums; the merge job writes the file
		classpoly_crt_dump (crt);
		classpoly_crt_end (crt);
		info_printf ("ECRT worker %d of %d dumped partial sums\n", crt->jobid, crt->jobs);
		return 1;
	}

	fp = fopen (filename, "w");
	if ( ! fp ) { err_printf ("Error opening file output file %s\n", filename); abort(); }
	if ( ! classpoly_crt_finish (crt, &outbits, &cbits, fp) ) { err_printf("call to classpoly_crt_finish failed!"); fclose(fp); remove (filename); abort(); }
	fclose (fp);
	end = clock();
	if ( tower_on ) out_printf ("Decomposition: F of degree %ld (cosets), G of degree %ld (subgroup)%s\n", tower->m, tower->n, tower_signfree ? ", sign-free" : "");
	if ( P ) {
		out_printf ("Class polynomial for inv=%d, D=%ld reduced mod P (%ld bits) written to %s, degree %d (%.1fs)\n", inv, D, mpz_sizeinbase(P,2), filename, H_d, (double)delta_msecs(begin,end)/1000.0);
	} else {
		out_printf ("Class polynomial for  D=%ld written to %s, degree %d, height %ld (%.0f), size %.3f MB (%.1fs)\n", D,  filename, H_d, cbits, bits-HEIGHT_MARGIN, (double)outbits/8000000.0, (double)delta_msecs(begin,end)/1000.0);
		if ( bits+height_factor-cbits < HEIGHT_MARGIN/2 ) {
			char *badfile = mem_alloc(strlen(filename)+5);
			strcpy(badfile,filename);  strcat(badfile,".bad");
			rename (filename, badfile);
			err_printf("*** Exceeded height safety margin!  Renamed output file to %s ***\n", badfile);
			mem_free(badfile);
			return 0;
		}
	}
	info_printf ("Total time for %d primes (%.1f bits) with H_d=%d was %ld msecs\n", crt_pcnt, tbits, H_d, delta_msecs(begin,end));
	tot_cycles = tot_fc_cycles+tot_er_cycles+tot_pr_cycles;	
	info_printf ("%.3f find (%.3f)   %.3f enum (%.3f) diff%.3f poly (%.3f)\n", (double)tot_fc_cycles/1000000000, (double)tot_fc_cycles/tot_cycles, (double)tot_er_cycles/1000000000,
	                 (double)tot_er_cycles/tot_cycles, (double)tot_pr_cycles/1000000000, (double)tot_pr_cycles/tot_cycles);
	return 1;
}
