\name{logpriors}
\alias{logpriors}
\title{
Calculate the log prior used in the log marginal density of Chib (1995).
}
\description{
Log marginal likelihood = Log likelihood + Log prior - Log density
}
\usage{
logpriors(h2, data_x, prior_p, prior_st)
}
\arguments{
  \item{h2}{Square of re-parameterized bandwidths and square of normal error variance}
  \item{data_x}{Regressors}
  \item{prior_p}{Hyperparameter used in the inverse-gamma prior}
  \item{prior_st}{Hyperparameter used in the inverse-gamma prior}
}
\details{
Calculate the log prior using the estimated averaged bandwidths of the regressors and the estimated 
averaged variance of the error density, obtained from the MCMC iterations
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
\code{\link[bbemkr]{logdensity}}, \code{\link[bbemkr]{loglikelihood}}, \code{\link[bbemkr]{mcmcrecord}}
}
\keyword{methods}

