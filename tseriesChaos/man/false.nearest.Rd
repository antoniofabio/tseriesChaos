\name{false.nearest}
\alias{false.nearest}
\title{Method of false nearest neighbours}
\description{
Method of false nearest neghbours to help deciding the optimal embedding dimension
}
\usage{
false.nearest(series, m, d, t, rt=10, eps=sd(series)/10)
}
\arguments{
  \item{series}{time series}
  \item{m}{maximum embedding dimension}
  \item{d}{delay parameter}
  \item{t}{Theiler window}
  \item{rt}{escape factor}
  \item{eps}{neighborhood diameter}
}
\details{
Method of false nearest neighbours to help deciding the optimal embedding dimension.
}
\value{
Fraction of false neighbors (first row) and total number of neighbors (second row) for each specified embedding dimension (columns)
}
\examples{
(fn.out <- false.nearest(rossler.ts, m=6, d=8, t=180, eps=1, rt=3))
plot(fn.out)
}
\references{
Hegger, R., Kantz, H., Schreiber, T., Practical implementation of nonlinear time series methods: The TISEAN package; CHAOS 9, 413-435 (1999)

Kennel M. B., Brown R. and Abarbanel H. D. I., Determining embedding dimension for phase-space reconstruction using a geometrical construction, Phys. Rev. A, Volume 45, 3403 (1992).
}
\author{Antonio, Fabio Di Narzo}
\keyword{ ts }
