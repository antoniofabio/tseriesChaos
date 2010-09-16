/*Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:09:44 +0100 (ven, 02 dic 2005) $*/
#include "tseriesChaos.h"

/*
Sample correlation integral for multiple length scales and multiple embedding dimensions.
in_series: input time series
in_length: time series length
in_m, in_d, in_t: max embedding dimension, time delay and theiler window
in_neps: number of length scales to evaluate
in_epsM: max length scale
in_epsm: min length scale
out: matrix of results
*/
void d2(double *in_series, int *in_length, int *in_m, int *in_d, int *in_t, int *in_neps, double *in_epsM, double *in_epsm, double *out){
	double tmpd, **hist;
	int i,j,w;
	int length, m,d,t, neps, blength;
	double *series, epsM, epsm;
	double a, lepsM;

/*
BIND PARAMETERS
*/
	series = in_series;
	length = *in_length;
	m = *in_m;
	d = *in_d;
	t = *in_t;
	neps = *in_neps;
	epsm = sqr(*in_epsm);
	epsM = sqr(*in_epsM);
/**/
/*
INIT VARIABLES
*/
	blength = length -(m-1)*d;
	lepsM = log(epsM);
	a = log(epsm/epsM)/(double)(neps-1);
	hist = (double**) R_alloc(m, sizeof(double*));
	for(i=0; i<m; i++) {
		hist[i] = (double*) R_alloc(neps, sizeof(double));
		for(j = 0; j<neps; j++)
			hist[i][j] = 0.0;
	}
/**/

	for(i = 0; i<(blength-t); i++) { /*for each point...*/
		R_CheckUserInterrupt();
		for(j=i+t; j<blength; j++) { /*for each upper-right point...*/
			tmpd = 0.0; /*init distance to 0*/
			for(w=0; w<m; w++) { /*for each dimension...*/
				tmpd += sqr(series[i+w*d] - series[j+w*d]); /*update squared euclidean distance*/
			        int ind = (log(tmpd) - lepsM)/a; /*FIX thanks to prof. B. Ripley*/
                	        hist[w][MIN(MAX(ind, 0), neps-1)]++; /*update histogram for current dimension*/
			} /*end for each dimension*/
		} /*end for each upper-right point*/
	} /*end for each point*/

	MATRIX2VEC(hist, out, m, neps); /*copy result to output*/
}
