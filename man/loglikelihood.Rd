\name{loglikelihood}
\alias{loglikelihood}
\title{
Calculate the log likelihood used in the Chib's (1995) log marginal density
}
\description{
Log marginal likelihood = Log likelihood + Log prior - Log density
}
\usage{
loglikelihood(h2, data_x)
}
\arguments{
  \item{h2}{Square of re-parameterized bandwidths and square of normal error variance}
  \item{data_x}{Regressors}
}
\details{
Calculates the log likelihood using the estimated averaged bandwidths of the regressors and estimated 
averaged variance of the error density
}
\value{
The value of log likelihood, with parameters (bandwidths + normal error variance) estimated from the MCMC iterations
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
\code{\link[bbemkr]{logpriors}}, \code{\link[bbemkr]{logdensity}}, \code{\link[bbemkr]{mcmcrecord}}
}
\keyword{methods}

