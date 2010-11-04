\name{bbecost}
\alias{bbecost}
\title{
Negative of log posterior
}
\description{
Calculates the negative of log posterior, using the leave-one-out samples.
}
\usage{
bbecost(data_x, data_y, x, bandx_priors, prior_p = 2, prior_st = 0.1)
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Bandwidths of the regressors}
  \item{bandx_priors}{Prior information}
  \item{prior_p}{Tuning parameter of the prior following inverse gamma distribution}
  \item{prior_st}{Tuning parameter of the prior following inverse gamma distribution}
}
\details{
The default bandx_priors is the inverse gamma distribution.
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
