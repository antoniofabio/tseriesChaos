/*Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:09:44 +0100 (ven, 02 dic 2005) $*/
#include "tseriesChaos.h"

/*
Follows pre-computed neighbours in a time series, and averages out to obtain the 'streching factor' in time.
in_series: time series
in_m, in_d: embedding dimension and time delay
in_length: time series length
in_nref: number of reference points to follow
in_totref: total number of reference points
in_k: number of neighbours for each reference point
in_s: time steps
in_nearest: neighbours matrix
in_ref: indexes of reference points to follow
lyap: computed stretching factor in time
*/
void follow_points(double *in_series, int *in_m, int *in_d, 
	int *in_length, int *in_nref, int *in_totref, int *in_k, 
	int *in_s, int *in_nearest, int *in_ref, double *lyap){

double *series; 
int m,d, s, nref, totref, k, length, *ref;
int i,j,a,b,md, time;
double tmp, res;
int **nearest;

/*
BIND PARAMETERS
*/
	m = *in_m;
	d = *in_d;
	s = *in_s;
	nref=*in_nref;
	totref=*in_totref;
	ref = in_ref;
	k = *in_k;
	series=in_series;
	length=*in_length;
/**/
/*
INIT VARIABLES
*/
	nearest= (int**) R_alloc(totref, sizeof(int*));
	for(i = 0; i<totref; i++) {
		nearest[i] = (int*) R_alloc(k, sizeof(int));
		for (j=0; j<k; j++) 
			nearest[i][j] = in_nearest[INDEX(i, j, totref)];
	}
	for(j=0; j<s; j++) lyap[j] = 0.0;
	md = m*d;
/**/

	for(time=0; time<s; time++) { /*for each time step...*/
		for(i=0; i<nref; i++) { /*for each reference point...*/
			tmp = 0.0; /*init sum of distances*/
			for(j=0; j<k; j++) { /*for each neighbour...*/
				a = ref[i]+time-1; /*pointer to the reference point at current time*/
				b = nearest[ref[i]-1][j]+time-1; /*pointer to current neighbour at current time*/
				DIST2(series, a, b, md, d, res);
				tmp += sqrt(res); /*add distance*/
			} /*end for each neighbour*/
			lyap[time] += log(tmp/(double)k); /*add to streching factor at current time, the log-mean-"sum of distances"*/
		} /*end for each reference point*/
		lyap[time] /= (double)nref; /*divide streching factor at current time by the total number of reference points*/
	} /*end for each time step*/
}
