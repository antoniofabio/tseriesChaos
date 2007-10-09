\name{C2}
\alias{C2}
\title{Sample correlation integral}
\description{
Sample correlation integral for the specified length scale
}
\usage{
C2(series, m, d, t, eps)
}
\arguments{
  \item{series}{time series}
  \item{m}{embedding dimension}
  \item{d}{time delay}
  \item{t}{Theiler window}
  \item{eps}{length scale}
}
\details{
Computes the sample correlation integral on the provided time series for the specified length scale, and considering a time window \code{t} (see references). It  uses a naif algorithm: simply returns the fraction of points pairs nearer than \code{eps}. Normally, you would use \code{\link{d2}}, which takes roughly the same time, but computes the correlation sum for multiple length scales and embedding dimensions at once.
}
\value{
The sample correlation integral at \code{eps} length scale.
}
\seealso{
\code{\link{d2}}
}
\references{
Hegger, R., Kantz, H., Schreiber, T., Practical implementation of nonlinear time series methods: The TISEAN package; CHAOS 9, 413-435 (1999)
}
\author{Antonio, Fabio Di Narzo}
\keyword{ts}
