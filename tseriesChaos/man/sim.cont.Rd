\name{sim.cont}
\alias{sim.cont}
\title{Simulates a continuous dynamic system}
\description{
Simulates a dynamic system of provided ODEs
}
\usage{
sim.cont(syst, start.time, end.time, dt, start.x, parms=NULL, obs.fun=function(x) x[1])
}
\arguments{
  \item{syst}{ODE system}
  \item{start.time}{starting time}
  \item{end.time}{ending time}
  \item{dt}{time between observations}
  \item{start.x}{initial conditions}
  \item{parms}{parameters for the system}
  \item{obs.fun}{observed function of the state}
}
\details{
Simulates a dynamic system of provided ODEs.
Uses \code{lsoda} in \code{odesolve} for numerical integration of the system.
}
\seealso{
  \code{\link{lorenz.syst}}, \code{\link{rossler.syst}}, \code{\link{duffing.syst}}
}
\value{
The time series of the observed function of the system's state
}
\examples{
rossler.ts <- sim.cont(rossler.syst, start=0, end=650, dt=0.1, start.x=c(0,0,0), parms=c(0.15, 0.2, 10))
}

\author{Antonio, Fabio Di Narzo}
\keyword{datagen}
