\name{Lyapunov exponent}
\alias{lyap}
\alias{lyap_k}
\title{Tools to evaluate the maximal Lyapunov exponent of a dynamic system}
\description{
Tools to evaluate the maximal Lyapunov exponent of a dynamic system from a univariate time series
}
\usage{
lyap_k(series, m, d, t, k=1, ref, s, eps)
lyap(dsts, start, end)
}
\arguments{
  \item{series}{time series}
  \item{m}{embedding dimension}
  \item{d}{time delay}
  \item{k}{number of considered neighbours}
  \item{eps}{radius where to find nearest neighbours}
  \item{s}{iterations along which follow the neighbours of each point}
  \item{ref}{number of points to take into account}
  \item{t}{Theiler window}
  \item{dsts}{Should be the output of a call to \code{lyap_k} (see the example)}
  \item{start}{Starting time of the linear bite of \code{dsts}}
  \item{end}{Ending time of the linear bite of \code{dsts}}
}
\details{
  The function \code{lyap_k} estimates the largest Lyapunov exponent of a given scalar time series using the algorithm of Kantz.

  The function \code{lyap} computes the regression coefficients of a user specified segment of the sequence given as input.
}
\value{
\code{lyap_k} gives the  logarithm of the stretching factor in time.
 
\code{lyap} gives the regression coefficients of the specified input sequence.
}
\seealso{
  \code{\link{mutual}}, \code{\link{false.nearest}} for the choice of optimal embedding parameters.
  \code{\link{embedd}} to perform embedding.
}
\examples{
output <-lyap_k(lorenz.ts, m=3, d=2, s=200, t=40, ref=1700, k=2, eps=4)
plot(output)
lyap(output, 0.73, 2.47)
}
\references{
Hegger, R., Kantz, H., Schreiber, T., Practical implementation of nonlinear time series methods: The TISEAN package; CHAOS 9, 413-435 (1999)

M. T. Rosenstein, J. J. Collins, C. J. De Luca, A practical method for calculating largest Lyapunov exponents from small data sets, Physica D 65, 117 (1993)
}
\author{Antonio, Fabio Di Narzo}
\keyword{ts}
