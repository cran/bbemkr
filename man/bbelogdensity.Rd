\name{bbelogdensity}
\alias{bbelogdensity}
\title{
Calculate an estimate of log posterior ordinate used in the log marginal density of Chib (1995). 
}
\description{
The calculation of log posterior ordinate estimate does not suffer from any instability problem, because it is a density value that is
averaged rather than its inverse.
}
\usage{
bbelogdensity(x, sigma, data_post)
}
\arguments{
  \item{x}{Bandwidths of the regressors}
  \item{sigma}{Variance of the error density}
  \item{data_post}{Gibbs output obtained from the MCMC iterations}
}
\details{
It should be noted that a modal value such as the posterior mode, or the maximum likelihood estimate, can be computed from the Gibbs output, 
at least approximately, if it is easy to evaluate the log-likelihood function for each draw in the simulation. 
Alternatively, one can make use of the posterior mean provided that there is no concern that it is a low density point.
}
\value{
Value of the log density.
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
\code{\link[bbemkr]{bbelogpriors}}, \code{\link[bbemkr]{bbeloglikelihood}}, \code{\link[bbemkr]{bbeMCMCrecording}}
}
\keyword{methods}

