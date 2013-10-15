\name{data_yt}
\alias{data_yt}
\docType{data}
\title{
Simulated response variable
}
\description{
The response variable is simulated from the functional form of
\eqn{data_yt = sin(2*pi*x1) + 4*(1-x2)*(1+x2) + 2*x3/(1+0.8*x3*x3) + rt(1,df=4)}, where \code{x1}, \code{x2} and \code{x3} are
simulated from a uniform distribution between 0 and 1.
}
\usage{data(data_yt)}
\format{
A data matrix of 100 by 1
}
\examples{
data(data_yt)
}
\keyword{datasets}

