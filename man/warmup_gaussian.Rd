\name{warmup_gaussian}
\alias{warmup_gaussian}
\title{
Burn-in period
}
\description{
By minimizing the cost value, the function estimates the bandwidths of the regressors and normal error variance
parameter for the burn-in period
}
\usage{
warmup_gaussian(x, inicost, mutsizp, warm = 100, prob = 0.234, data_x, data_y,
       prior_p = 2, prior_st = 1)
}
\arguments{
  \item{x}{Log of square bandwidths}
  \item{inicost}{Cost value}
  \item{mutsizp}{Step size of random-walk Metropolis algorithm}
  \item{warm}{Number of burn-in iterations}
  \item{prob}{Optimal covergence rate of random-walk Metropolis algorithm}
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{prior_p}{Hyperparameter of the inverse-gamma prior}
  \item{prior_st}{Hyperparameter of the inverse-gamma prior}
}
\value{
    \item{x}{Log of square bandwidths}
    \item{sigma2}{Estimate of normal error variance}
    \item{cost}{Cost value}
    \item{mutsizplast}{Final step size of random-walk Metropolis algorithm}
    \item{mutsizp}{Step size of random-walk Metropolis algorithm}
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{mcmcrecord_gaussian}}, \code{\link[bbemkr]{logdensity_gaussian}}, \code{\link[bbemkr]{loglikelihood_gaussian}}, \code{\link[bbemkr]{logpriors_gaussian}}
}
\keyword{methods}
