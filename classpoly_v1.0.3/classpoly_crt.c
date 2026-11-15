#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "crt.h"
#include "classpoly_crt.h"
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

void classpoly_crt_start (classpoly_crt_t crt, long D, int inv, unsigned long *primes, int pcnt, int ccnt, mpz_t P)
{
	clock_t start, end;
	char prefix[160], *s, *dir;
	int fcnt;
	
	memset (crt, 0, sizeof(crt[0]));
	// Multi-job ECRT (parallel workers, cf. ecrt_dump/ecrt_merge): CLASSPOLY_JOBS=W (>1)
	// makes this process worker CLASSPOLY_JOBID (1..W), handling only prime indices
	// index = jobid-1 mod W and dumping partial coefficient sums instead of writing
	// the class polynomial; a subsequent run with CLASSPOLY_JOBID=0 merges the dumps
	// (see compute_classpoly).  CLASSPOLY_ECRT_DIR is the shared directory for the
	// dump files (it must be common to all jobs; default is the cwd).
	if ( (s = getenv("CLASSPOLY_JOBS")) && atoi(s) > 1 ) {
		crt->jobs = atoi(s);
		crt->jobid = ( (s = getenv("CLASSPOLY_JOBID")) ? atoi(s) : 1 );
		if ( crt->jobid < 1 || crt->jobid > crt->jobs ) { err_printf ("CLASSPOLY_JOBID=%d out of range 1..%d\n", crt->jobid, crt->jobs); abort(); }
	}
	crt->D = D;
	crt->inv = inv;
	crt->p = primes;
	crt->pcnt = pcnt;
	crt->ccnt = ccnt;
	if ( P ) mpz_init_set (crt->P, P); else mpz_init (crt->P);
	if ( mpz_sgn(crt->P) ) {
		start = clock();
		dir = getenv("CLASSPOLY_ECRT_DIR");
		snprintf (prefix, sizeof(prefix), "%s%sH_%ld", dir ? dir : "", dir ? "/" : "", -D);
		ecrt_init (crt->ecrt, primes, pcnt, ccnt, crt->P, crt->jobs, crt->jobs ? crt->jobid : 0, prefix);
		end = clock();
		info_printf ("Explicit CRT precomputation completed in %ld msecs\n", delta_msecs(start,end));
	} else {
		if ( ! (fcnt=crt_create_files (&crt->fp, -D, 0, ccnt)) ) { err_printf("Call to crt_create_files failed\n"); abort(); }
		info_printf ("Created %d CRT output files\n", fcnt);				
	}
	crt->index = -1;
}

void classpoly_crt_end (classpoly_crt_t crt)
{
	if ( crt->fp ) crt_close_files (crt->fp, crt->ccnt);
	else if ( mpz_sgn(crt->P) ) ecrt_clear (crt->ecrt);		// the explicit-CRT context was allocated iff P != 0
	mpz_clear (crt->P);
	if ( crt->jcache ) mem_free (crt->jcache);
	mem_free (crt->p);
}

void classpoly_crt_dump (classpoly_crt_t crt)
{
	if ( ! crt->jobs || ! mpz_sgn(crt->P) ) { err_printf ("classpoly_crt_dump requires multi-job ECRT mode\n"); abort(); }
	crt->ecrt->tag = classpoly_tower_tag (crt->tower_m, crt->tower_signfree);		// the merge job learns the shape of the output from the dump files
	ecrt_dump (crt->ecrt);
}

/*
	Output format.  A plain class polynomial is written as before: I=, D= and (if reduced) P= lines, then one line "c*X^i + " per
	coefficient and a final "1*X^deg".

	A decomposition (tower_m = m > 0, n = ccnt/m) is written in the F/G format that classpoly_load and ff_classpoly_root already
	understand: after the I=, D=, P= lines come "F=m" and "G=n", then m*n+1 coefficient lines:
		F_0, ..., F_{m-1}, 1				F(Y) = Y^m + F_{m-1}Y^{m-1} + ... + F_0 has the m coset values y as its roots
		g_{0,0}, ..., g_{0,m-1}				g_0(Y) = g_{0,0} + g_{0,1}Y + ...
		...
		g_{n-2,0}, ..., g_{n-2,m-1}
	and for a root y of F the polynomial G(X) = X^n + y X^{n-1} + sum_{j<=n-2} (g_j(y)/F'(y)) X^j divides the class polynomial.
	Each line is "c*<monomial>"; only the part before '*' is significant to a reader.

	Sign-free decompositions.  For an invariant f whose class polynomial is only determined mod p up to f -> -f, fixing the sign takes
	the trace T of the class polynomial, which costs a CRT of its own before the main one (classpoly_inv_setup).  A consumer who only
	wants a root up to sign (the invariants concerned are square roots of invariants through which j factors) can do without: under
	f -> -f the coefficient of Y^j in F is multiplied by (-1)^(m-j) and that of Y^j in g_k by (-1)^(n-k+m-1-j), so each coefficient
	has a parity, and T itself is odd.  In a sign-free decomposition every odd coefficient is stored multiplied by T, which makes all
	stored numbers independent of the sign, and S = T^2 is stored as well, in a header line "S=<S>" after G=.  Dividing the odd
	coefficients by either square root of S gives the decomposition of H(X) or of H(-X).  (The line starts with a letter, so
	classpoly_load skips it; the file must not be used as a plain F/G file.)
*/
void classpoly_crt_print_header (FILE *fp, long D, int inv, mpz_t P, int ccnt, int tower_m, int signfree)
{
	fprintf (fp, "I=%d\n", inv);
	fprintf (fp, "D=%ld\n", D);
	if ( P && mpz_sgn(P) ) gmp_fprintf (fp, "P=%Zd\n", P);
	if ( tower_m ) fprintf (fp, "F=%d\nG=%d\n", tower_m, (ccnt-(signfree?1:0))/tower_m);
}

void classpoly_crt_print_coeff (FILE *fp, mpz_t c, int i, int ccnt, int tower_m, int signfree)
{
	if ( tower_m && signfree ) {
		if ( ! i ) { gmp_fprintf (fp, "S=%Zd\n", c); return; }
		i--;
	}
	if ( ! tower_m ) {
		gmp_fprintf (fp, "%Zd*X^%d + \n", c, i);
		if ( i == ccnt-1 ) fprintf (fp, "1*X^%d\n", ccnt);
	} else if ( i < tower_m ) {
		gmp_fprintf (fp, "%Zd*Y^%d\n", c, i);
		if ( i == tower_m-1 ) fprintf (fp, "1*Y^%d\n", tower_m);
	} else {
		gmp_fprintf (fp, "%Zd*Y^%d*X^%d\n", c, (i-tower_m)%tower_m, (i-tower_m)/tower_m);
	}
}

struct classpoly_crt_output_ctx { FILE *fp; int ccnt, tower_m, signfree; };

int classpoly_crt_output_coefficient (mpz_t c, int i, void *ctx)
{
	struct classpoly_crt_output_ctx *o = (struct classpoly_crt_output_ctx *) ctx;
	classpoly_crt_print_coeff (o->fp, c, i, o->ccnt, o->tower_m, o->signfree);
	return 1;
}

int classpoly_crt_finish (classpoly_crt_t crt, long *totbits, long *maxcbits, FILE *fp)
{
	mpz_t X;
	clock_t start, end;
	register int i;

	assert (crt->index == crt->pcnt);
	if ( totbits ) *totbits = 0;
	if ( maxcbits ) *maxcbits = 0;
	start = clock();
	if ( crt->fp ) {
		crt_close_files(crt->fp, crt->ccnt);
		info_printf ("Wrote %ld bytes of data to CRT files\n", crt->ccnt*crt->pcnt*sizeof(unsigned long));
	}
	classpoly_crt_print_header (fp, crt->D, crt->inv, crt->P, crt->ccnt, crt->tower_m, crt->tower_signfree);

	if ( mpz_sgn(crt->P) ) {		
		// all ecrt cases treated as "small" in the sense of section 6 of "Computing Hilbert class polynomials with the Chinese Remainder Theorem")
		ecrt_finalize (crt->ecrt);
		mpz_init (X);
		for ( i = 0 ; ecrt_next_coeff(X, crt->ecrt) ; i++ ) classpoly_crt_print_coeff (fp, X, i, crt->ccnt, crt->tower_m, crt->tower_signfree);
		mpz_clear (X);
		if ( i != crt->ccnt ) { err_printf ("ECRT postcomputation returned only %d of %d expected coefficients\n", i, crt->ccnt); ecrt_clear (crt->ecrt); mpz_clear (crt->P); return 0; }
		if ( totbits ) *totbits = crt->ccnt * mpz_sizeinbase(crt->P,2);
		if ( maxcbits ) *maxcbits = mpz_sizeinbase(crt->P,2);
		ecrt_clear (crt->ecrt);
		mpz_clear (crt->P);
	} else {
		struct classpoly_crt_output_ctx octx = { fp, crt->ccnt, crt->tower_m, crt->tower_signfree };
		mpz_clear (crt->P);
		if ( ! crt_process_files (totbits, maxcbits, -crt->D, 0, crt->ccnt, crt->p, crt->pcnt, classpoly_crt_output_coefficient, &octx) ) { err_printf ("Error processing CRT output files\n");  return 0; }
	}
	end = clock();
	if ( crt->jcache ) mem_free (crt->jcache);
	mem_free (crt->p);
	info_printf ("Completed CRT postcomputation in %ld msecs\n", delta_msecs(start,end));
	return 1;
}
