\name{data_y}
\alias{data_y}
\docType{data}
\title{
Simulated response variable
}
\description{
The response variable is simulated from the functional form of
\eqn{data_y = sin(2*pi*x1) + 4*(1-x2)*(1+x2) + 2*x3/(1+0.8*x3*x3) + rnorm(1,0,0.9)}, where \code{x1}, \code{x2} and \code{x3} are
simulated from a uniform distribution between 0 and 1.
}
\usage{data(data_y)}
\format{
A data matrix of 100 by 1
}
\examples{
data(data_y)
}
\keyword{datasets}

