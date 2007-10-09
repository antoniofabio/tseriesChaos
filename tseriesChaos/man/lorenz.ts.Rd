\name{lorenz.ts}
\alias{lorenz.ts}
\title{Lorenz simulated time series, without noise}
\description{
Lorenz simulated time series, without noise. Of each state of the system, we observe the euclidean norm.
}
\details{
Lorenz simulated time series, without noise, obtained with the call:
\code{lorenz.ts <- sim.cont(lorenz.syst, 0, 100, 0.05, start.x=c(5,5,5), parms=c(10, 28, -8/3), obs.fun = function(x) sqrt(sum(x^2)))}
}
\author{Antonio, Fabio Di Narzo}
\keyword{datasets}
