\name{bbeloglikelihood}
\alias{bbeloglikelihood}
\title{
Calculate the log likelihood used in the Chib's (1995) log marginal density
}
\description{
Calculates the log likelihood using the estimated averaged bandwidths of the regressors and estimated 
averaged variance of the error density}
\usage{
bbeloglikelihood(data_x, data_y, x, sigma, kerntype = c("Gaussian", 
    "Epanechnikov", "Quartic", "Triweight", "Triangular", "Uniform"))
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Estimated averaged bandwidths of the regressors, obtained after the MCMC iterations}
  \item{sigma}{Estimated averaged variance of the error density, obtained after the MCMC iterations}
  \item{kerntype}{Type of the kernel function. By default, Gaussian kernel is used}
}
\details{
According to Chib (1995), the log marginal density = loglikelihood + logprior - logdensity
}
\value{
The value of log likelihood, with parameters estimated from the MCMC iterations
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
\code{\link[bbemkr]{bbelogpriors}}, \code{\link[bbemkr]{bbelogdensity}}, \code{\link[bbemkr]{bbeMCMCrecording}}
}
\keyword{methods}
