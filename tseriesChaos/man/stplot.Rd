\name{stplot}
\alias{stplot}
\title{Space-time separation plot}
\description{
Space-time separation plot
}
\usage{
stplot(series, m, d, idt=1, mdt)
}
\arguments{
  \item{series}{time series}
  \item{m}{embedding dimension}
  \item{d}{time delay}
  \item{idt}{observation steps in each iteration}
  \item{mdt}{number of iterations}
}
\details{
Produces the space-time separation plot, as introduced by Provenzale et al. (1992), which can be used to decide the Theiler time window \code{t}, which is required in many other algorithms in this package.

It plots the probability that two points in the reconstructed phase-space have distance smaller than epsilon in function of epsilon and of the time t between the points, as iso-lines at levels 10\%, 20\%, \ldots, 100\%.
}
\value{
lines of costant probability at 10\%, 20\%, \ldots, 100\%.
}
\seealso{
\code{\link{false.nearest}}, \code{\link{d2}}, \code{\link{lyap_k}}
}
\examples{
stplot(rossler.ts, m=3, d=8, idt=1, mdt=250)
}
\references{
Kantz H., Schreiber T., Nonlinear time series analysis. Cambridge University Press, (1997)

Provenzale A., Smith L. A., Vio R. and Murante G., Distiguishing between low-dimensional dynamics and randomness in measured time series. Physica D., volume 58, 31 (1992)
}
\author{Antonio, Fabio Di Narzo}
\keyword{ts}
