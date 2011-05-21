\name{bbecost2}
\alias{bbecost2}
\title{
Negative of log posterior associated with the error variance
}
\description{
Calculates the negative of log posterior, using the leave-one-out samples.
}
\usage{
bbecost2(data_x, data_y, x, kerntype = c("Gaussian", "Epanechnikov", 
         "Quartic", "Triweight", "Triangular", "Uniform"), 
         prior_st = 0.1)
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Log bandwidths of the regressors}
  \item{kerntype}{Type of kernel function. By default, Gaussian kernel is used.}
  \item{prior_st}{A tuning parameter of the prior associated with the error variance, which follows 
				  inverse gamma distribution}
}
\details{
The default bandx_priors is the inverse gamma distribution. Assuming the error density is independent 
and identically distribution and follows a normal distribution, the function is designed to calculate 
the negative of log posterior associated with the error variance component.}
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
