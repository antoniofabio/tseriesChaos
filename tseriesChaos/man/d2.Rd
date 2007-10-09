\name{d2}
\alias{d2}
\title{Sample correlation integral (at multiple length scales)}
\description{
Computes the sample correlation integral over a grid of \code{neps} length scales starting from \code{eps.min}, and for multiple embedding dimensions
}
\usage{
d2(series, m, d, t, eps.min, neps=100)
}
\arguments{
  \item{series}{time series}
  \item{m}{max embedding dimension}
  \item{d}{time delay}
  \item{t}{Theiler window}
  \item{eps.min}{min length scale}
  \item{neps}{number of length scales to evaluate}
}
\details{
Computes the sample correlation integral over \code{neps} length scales starting from \code{eps.min}, for embedding dimension 1,\ldots,m , considering a \code{t} time window (see references). The slope of the linear segment in the log-log plot gives an estimate of the correlation dimension (see the example).
}
\value{
Matrix.
Column 1: length scales. Column i=2, \ldots, m+1: sample correlation integral for embedding dimension i-1.
}
\examples{
d2(lorenz.ts, m=6, d=2, t=4, eps.min=2)
}
\references{
Hegger, R., Kantz, H., Schreiber, T., Practical implementation of nonlinear time series methods: The TISEAN package; CHAOS 9, 413-435 (1999)
}
\author{Antonio, Fabio Di Narzo}
\keyword{ts}
