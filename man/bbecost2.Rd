\name{bbecost2}
\alias{bbecost2}
\title{
Negative of log posterior
}
\description{
Calculates the negative of log posterior, using the leave-one-out samples.
}
\usage{
bbecost2(data_x, data_y, x, prior_st = 0.1)
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Log bandwidths of the regressors}
  \item{prior_st}{Tuning parameter of the prior following inverse gamma distribution}
}
\details{
The default bandx_priors is the inverse gamma distribution. Assuming the error density is i.i.d and follows normal distribution,
the cost2 function is designed to calculate the variance parameter of the error density.
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
\code{\link[bbemkr]{np_gibbs}}, \code{\link[bbemkr]{bbecost}}
}
\keyword{methods}
