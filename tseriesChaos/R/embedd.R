#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
embedd <- function(x, m, d, lags) {
	x <- as.ts(x)
	if(missing(lags)) {
		checkEmbParms(x, m, d)
		lags <- ((1:m)-1)*d
	}
	res <- lag(x, lags[1])
	for(i in 2:length(lags)) {
		res <- ts.intersect(res, lag(x, lags[i]))
	}
	res <- matrix(res, nr = nrow(res), nc = ncol(res))
	colnames(res) <- paste("lag", lags, sep="")
	res
}
