\name{bbelogpriors}
\alias{bbelogpriors}
\title{
Calculate the log prior used in the log marginal density of Chib (1995).
}
\description{
Calculate the log prior using the estimated bandwidths of the regressors and variance of the error density obtained from the MCMC output
}
\usage{
bbelogpriors(x, sigma, prior_p = 2, prior_st = 0.1)
}
\arguments{
  \item{x}{Estimated bandwidths of the regressors, after the MCMC iterations}
  \item{sigma}{Estimated variance of the error density, after the MCMC iterations}
  \item{prior_p}{Tuning parameter in the inverse-gamma prior}
  \item{prior_st}{Tuning parameter in the inverse-gamma prior}
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
\code{\link[bbemkr]{bbelogdensity}},\code{\link[bbemkr]{bbeloglikelihood}},\code{\link[bbemkr]{bbeMCMCrecording}}
}
\keyword{methods}
