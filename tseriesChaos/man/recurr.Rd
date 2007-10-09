\name{recurr}
\alias{recurr}
\title{Recurrence plot}
\description{
Recurrence plot
}
\usage{
recurr(series, m, d, start.time=start(series), end.time=end(series), ...)
}
\arguments{
  \item{series}{time series}
  \item{m}{embedding dimension}
  \item{d}{time delay}
  \item{start.time}{starting time window (in time units)}
  \item{end.time}{ending time window (in time units)}
  \item{...}{further parameters to be passed to \code{filled.contour}}
}
\details{
Produces the recurrence plot, as proposed by Eckmann et al. (1987). 
White is maximum distance, black is minimum.
}
\section{warning}{
Be awared that number of distances to store goes as \code{n^2}, where \code{n = length(window(series, start=start.time, end=end.time))}!
}
\examples{
recurr(lorenz.ts, m=3, d=2, start.time=15, end.time=20)
}
\references{
Eckmann J.P., Oliffson Kamphorst S. and Ruelle D., Recurrence plots of dynamical systems, Europhys. Lett., volume 4, 973 (1987)
}
\author{Antonio, Fabio Di Narzo}
\keyword{ts}
