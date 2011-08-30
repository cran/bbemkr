\name{cost}
\alias{cost}
\title{
Negative of log posterior associated with the bandwidths
}
\description{
Calculates the negative of log posterior, using the leave-one-out cross validated samples.
}
\usage{
cost(x, data_x, data_y, prior_p, prior_st)
}
\arguments{
  \item{x}{Log of square bandwidths}
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{prior_p}{A tuning parameter of the prior of error variance, following inverse gamma distribution}
  \item{prior_st}{Another tuning parameter of the prior of error variance, following inverse gamma distribution}
}
\details{
Bandwidth can be re-parameterized by a constant times optimal convergence rate, that is, \eqn{h=c*n^{rate}}. The prior of \eqn{c^2} is 
assumed to follow an inverse-gamma prior with hyperparameters \code{prior_p = 2} and \code{prior_st = 1}.
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
\code{\link[bbemkr]{np_gibbs}}, \code{\link[bbemkr]{cost2}}
}
\examples{
x = log(nrr(data_x, FALSE)^2)
inicost = cost(x, data_x, data_y, prior_p = 2, prior_st = 1)
}
\keyword{methods}

