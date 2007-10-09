/*Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:09:44 +0100 (ven, 02 dic 2005) $*/
#include "tseriesChaos.h"

#define output(a,b) out[INDEX(a, b, ref)]

/*
Find k nearest neighbours of all points of a time series.
in_series: time series (scaled between 0 and 1)
m, d, t: embedding dimension, time delay and theiler window
length: length of the time series
eps: neighborhood size
ref: number of points to consider
in_k: max number of neighbours to search for each point
in_s: number of unit times to omit from the end
out: matrix of nearest neighbours
*/
void find_knearests(double *in_series, int *in_m, int *in_d, int *in_t,
	int *in_length, double *in_eps, int *in_ref, int *in_k, int *in_s, int *out) {
double eps, *series; 
int m,d, t, s, ref, k, length, blength;
int i,j,md;
double *dsts;
int id; int *ids;
boxSearch bs;

/*
BIND PARAMETERS
*/
	m = *in_m;
	d = *in_d;
	t = *in_t;
	s = *in_s;
	ref=*in_ref;
	k = *in_k;
	eps=*in_eps;
	series=in_series;
	length=*in_length;
/**/
/*
INIT VARIABLES
*/
	blength = length - (m-1)*d - s;
	md = m*d;
	for(i = 0; i<ref; i++) 
		for(j=0; j<k; j++) 
			output(i,j) = -1;
	dsts = (double*) R_alloc(blength, sizeof(double));
	ids = (int*) R_alloc(blength, sizeof(int));
	bs = init_boxSearch(series, m, d, blength, eps);
/**/

	for(i = 0; i<ref; i++) { /*for each reference point...*/
		find_nearests(bs, t, i, ids, dsts, &id);
		R_qsort_I(dsts, ids, 1, id);
		for(j=0; (j<k) && (j<id); j++)
			output(i, j) = ids[j]+1;
	}/*end for each reference point*/
}
