\name{cost2}
\alias{cost2}
\title{
Negative of log posterior associated with the error variance
}
\description{
Calculates the negative of log posterior for the normal error variance, using the leave-one-out cross validated samples.
}
\usage{
cost2(x, data_x, data_y, prior_st = 1)
}
\arguments{
  \item{x}{Log of square bandwidths}
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{prior_st}{Another tuning parameter of the prior of error variance, following inverse gamma distribution}
}
\details{
The prior of normal error variance follows an inverse-gamma distribution with hyperparameter \code{prior_st = 1}.
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
\code{\link[bbemkr]{np_gibbs}}, \code{\link[bbemkr]{cost}}
}
\examples{
x = log(nrr(data_x, FALSE)^2)
cost2(x, data_x, data_y)
}
\keyword{methods}

