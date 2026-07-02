#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "ff_poly.h"
#include "phi_poly.h"
#include "class_inv.h"
#include "classpoly.h"
#include "qform.h"
#include "pickprimes.h"
#include "classpoly_crt.h"
#include "classpoly_inv.h"
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
			fp = fopen (filename, "w");
			if ( ! fp ) { err_printf ("Error opening output file %s\n", filename); return 0; }
			gmp_fprintf (fp, "I=%d\n", inv);
			gmp_fprintf (fp, "D=%ld\n", D);
			gmp_fprintf (fp, "P=%Zd\n", ecrt->P);
			mpz_init (X);
			for ( i = 0 ; ecrt_next_coeff (X, ecrt) ; i++ ) gmp_fprintf (fp, "%Zd*X^%d + \n", X, i);
			mpz_clear (X);
			gmp_fprintf (fp, "1*X^%d\n", i);
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
	classpoly_crt_start (crt, D, inv, (unsigned long *)crt_p, crt_pcnt, H_d, P);
	if ( ! classpoly_inv_setup (Hinv, H, inv, crt) ) { classpoly_clear (H); classpoly_crt_end (crt); return 0; }

	work = mem_alloc (2*H_d*sizeof(ff_t));
	
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
		ff_poly_from_roots (work, H->roots, H_d);
		cntr2 = get_cycles();
		pr_cycles += (unsigned long)(cntr2-cntr1);  tot_pr_cycles +=  (unsigned long)(cntr2-cntr1); 
		if ( dbg_level >= DEBUG_LEVEL ) { dbg_printf ("poly (p=%ld)  x^%d",_ff_p, H_d); for ( i = H_d-1 ; i >= 0 ; i-- ) dbg_printf ("+ %ld*x^%d ", _ff_get_ui(work[i]), i); dbg_printf("\n"); }
		
		// *** Alg 2 Step 4.b CRT update
		for ( i = 0 ; i < H_d ; i++ ) work[i] = _ff_get_ui(work[i]);		// convert out of Montgomery rep
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
