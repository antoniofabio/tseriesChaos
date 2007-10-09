/*Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:09:44 +0100 (ven, 02 dic 2005) $*/
#include "tseriesChaos.h"
#define MEPS 1000
#define MFRAC 10

/*Space time separation plot.
in_series: time series
in_length: time series length
in_m, in_d: embedding dimension and time delay
in_steps: total time
in_idt: number of time units in each step
in_epsmax: max length scale
out: computed iso-lines of the plot
*/
void stplot(double *in_series, int *in_length, int *in_m, int *in_d, int *in_steps, int *in_idt, double *in_epsmax, double *out) {
	double tmp, need;
	int i,j, a, b, md, is, ieps, length, blength, m, d, steps, idt;
	double epsmax, *series, *hist, **stp;

/*
BIND PARAMETERS
*/
	series = in_series;
	length = *in_length;
	m = *in_m;
	d = *in_d;
	md = m*d;
	steps = *in_steps;
	idt = *in_idt;
	epsmax = sqr(*in_epsmax);
/**/
/*
INIT VARIABLES
*/
	blength = length - (m-1)*d;
	stp = (double**) R_alloc(MFRAC, sizeof(double*));
	for(i=0; i<MFRAC; i++) stp[i] = (double*) R_alloc(steps, sizeof(double));
	hist = (double*) R_alloc(MEPS, sizeof(double));
/**/

	for(i=0; i<steps; i++) { /*for each time step...*/
		for(j=0; j<MEPS; j++) hist[j] = 0.0; /*init histogram for all eps values*/
		for(j=0; j<(blength-i*idt); j++) { /*for each point...*/
			a = j; b = j+i*idt;
			DIST2(series, a, b, md, d, tmp);
			hist[MIN((long)(tmp*MEPS/epsmax), MEPS-1)]++;
		} /*end for each point*/
		for(j=0; j<MFRAC; j++) { /*update iso-lines*/
			need = (blength - i*idt)*(j+1)/(double) MFRAC;
			for(is=0, ieps=0; ieps<MEPS && is<need; ieps++)
				is +=hist[ieps];
			stp[j][i] = ieps*(epsmax/(double)MEPS);
		} /*end update iso-lines*/
	} /*end for each time step*/
	for(i=0; i<steps; i++) for(j=0; j<MFRAC; j++) /*take sqrt on all iso-lines*/
		stp[j][i] = sqrt(stp[j][i]);

	MATRIX2VEC(stp, out, MFRAC, steps); /*copy result to the output*/
}
