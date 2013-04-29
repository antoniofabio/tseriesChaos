\name{plot.ami}
\alias{plot.ami}
\title{Plotting average mutual information index}
\description{
Plotting method for objects inheriting from class '"ami"'.
}
\usage{
 \method{plot}{ami}(x, main = NULL, ...)
}
\arguments{
  \item{x}{'"ami"' object}
  \item{main, ...}{additional graphical arguments}
}
\details{
Plots the ami for each lag in \code{x}.
}
\seealso{
\code{\link{mutual}}
}
\author{Antonio, Fabio Di Narzo}
\keyword{hplot}
