\name{mutual}
\alias{mutual}
\title{Average Mutual Information}
\description{
Estimates the average mutual information index (ami) of a given time series for a specified number of lags
}
\usage{
mutual(series, partitions = 16, lag.max = 20, plot=TRUE, ...)
}

\arguments{
  \item{series}{time series}
  \item{partitions}{number of bins}
  \item{lag.max}{largest lag}
  \item{plot}{logical. If 'TRUE' (the default) the ami is plotted}
  \item{...}{further arguments to be passed to the plot method}
}
\details{
Estimates the mutual information index for a specified number of lags.
The joint probability distribution function is estimated with a simple bi-dimensional density histogram.
}
\value{
An object of class "ami", which is a vector containing the estimated mutual information index for each lag between 0 and \code{lag.max}.
}
\examples{
mutual(lorenz.ts)
}
\references{
Hegger, R., Kantz, H., Schreiber, T., Practical implementation of nonlinear time series methods: The TISEAN package; CHAOS 9, 413-435 (1999)
}
\author{ Antonio, Fabio Di Narzo }
\keyword{ts}
