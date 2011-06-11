\name{bbecost}
\alias{bbecost}
\title{
Negative of log posterior associated with the bandwidths
}
\description{
Calculates the negative of log posterior, using the leave-one-out samples.
}
\usage{
bbecost(data_x, data_y, x, kerntype = c("Gaussian", "Epanechnikov", "Quartic", "Triweight", "Triangular", "Uniform"), bandx_priors, prior_p = 2, prior_st = 0.1)
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Log bandwidths of the regressors}
  \item{kerntype}{Type of kernel function. By default, Gaussian kernel is used}
  \item{bandx_priors}{Prior information of the bandwidths}
  \item{prior_p}{A tuning parameter of the prior of error variance, following inverse gamma distribution}
  \item{prior_st}{Another tuning parameter of the prior of error variance, following inverse gamma distribution}
}
\details{
The default bandx_priors is the inverse gamma distribution. However, Cauchy prior can also be used which 
achieve similar results.
}
\value{
Value of the cost function
}
\references{
X. Zhang and R.D. Brooks and M.L. King (2009), A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, \emph{Journal of Econometrics}, \bold{153}, 21-32. 
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{np_gibbs}}, \code{\link[bbemkr]{bbecost2}}
}
\keyword{methods}
