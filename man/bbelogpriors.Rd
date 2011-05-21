\name{bbelogpriors}
\alias{bbelogpriors}
\title{
Calculate the log prior used in the log marginal density of Chib (1995).
}
\description{
Calculate the log prior using the estimated averaged bandwidths of the regressors and the estimated 
averaged variance of the error density, obtained from the MCMC iterations
}
\usage{
bbelogpriors(x, sigma, prior_p = 2, prior_st = 0.1)
}
\arguments{
  \item{x}{Estimated averaged bandwidths of the regressors, obtained after the MCMC iterations}
  \item{sigma}{Estimated averaged variance of the error density, obtained after the MCMC iterations}
  \item{prior_p}{A tuning parameter of the prior of error variance, following inverse gamma distribution}
  \item{prior_st}{Another tuning parameter of the prior of error variance, following inverse gamma distribution}
}
\details{
According to Chib (1995), the log marginal density = loglikelihood + logprior - logdensity
}
\value{
Value of the log prior
}
\references{
S. Chib (1995) Marginal likelihood from the Gibbs output, \emph{Journal of the American Statistical Association}, \bold{90}, 432, 1313-1321.

M. A. Newton and A. E. Raftery (1994) Approximate Bayesian inference by the weighted likelihood bootstrap (with discussion), \emph{Journal of 
the Royal Statistical Society}, \bold{56}, 3-48.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{bbelogdensity}}, \code{\link[bbemkr]{bbeloglikelihood}}, \code{\link[bbemkr]{bbeMCMCrecording}}
}
\keyword{methods}
