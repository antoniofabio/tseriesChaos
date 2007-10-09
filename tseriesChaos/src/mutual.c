/*Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:09:44 +0100 (ven, 02 dic 2005) $*/
#include "tseriesChaos.h"

#define output(i, j) out_hist[INDEX(i, j, partitions)]

/*Computes an 'auto' double histogram from a time series.
in_series: time series (scaled between 0 and 1)
in_length: time series length
in_lag: time lag
in_partitions: number of partitions to make
out_hist: matrix containing the computed double histogram
*/
void mutual(double *in_series, int *in_length, int *in_lag, 
int *in_partitions, double *out_hist) {
	int partitions, length, lag;
	int ix, iy, binx, biny, i, j;
	double *series;

	series = in_series;
	length = *in_length;
	lag =*in_lag;
	partitions=*in_partitions;

	for(i =0; i<partitions; i++) 
		for(j=0; j<partitions; j++)
			output(i, j) = 0.0;

	for(ix = 0; ix < (length-lag); ix++) {
		iy = ix + lag;
		binx = MIN((int)(series[ix]*partitions),partitions-1);
		biny = MIN((int)(series[iy]*partitions),partitions-1);
		output(binx, biny) ++;
	}
}
