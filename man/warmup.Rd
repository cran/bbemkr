\name{warmup}
\alias{warmup}
\title{
Burn-in period
}
\description{
By minimizing the cost value, the function estimates the bandwidths of the regressors and normal error variance
parameter for the burn-in period
}
\usage{
warmup(x, inicost, mutsizp, warm = 100, prob = 0.234, data_x, data_y)
}
\arguments{
  \item{x}{Log of square bandwidths}
  \item{inicost}{Cost value}
  \item{mutsizp}{Step size of random-walk Metropolis algorithm}
  \item{warm}{Number of burn-in iterations}
  \item{prob}{Optimal covergence rate of random-walk Metropolis algorithm}
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
}
\value{
    \item{x}{Log of square bandwidths}
    \item{sigma2}{Estimate of normal error variance}
    \item{cost}{Cost value}
    \item{mutsizplast}{Final step size of random-walk Metropolis algorithm}
    \item{mutsizp}{Step size of random-walk Metropolis algroithm}
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{mcmcrecord}}, \code{\link[bbemkr]{logdensity}}, \code{\link[bbemkr]{loglikelihood}}, \code{\link[bbemkr]{logpriors}}
}
\keyword{methods}
