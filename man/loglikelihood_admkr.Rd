\name{loglikelihood_admkr}
\alias{loglikelihood_admkr}
\title{
Calculate the log likelihood used in the Chib's (1995) log marginal density
}
\description{
Log marginal likelihood = Log likelihood + Log prior - Log density
}
\usage{
loglikelihood_admkr(h2, data_x, data_y)
}
\arguments{
	\item{h2}{Square of re-parameterized bandwidths}
	\item{data_x}{Regressors}
	\item{data_y}{Response variable}
}
\details{
Calculates the log likelihood using the estimated averaged bandwidths of the regressors and estimated 
averaged variance of the error density
}
\value{
The value of log likelihood, with parameters (bandwidths) estimated from the MCMC iterations
}
\references{
S. Chib and I. Jeliazkov (2001) Marginal likelihood from the Metropolis-Hastings output, \emph{Journal of the American Statistical Association}, \bold{96}, 453, 270-281.

S. Chib (1995) Marginal likelihood from the Gibbs output, \emph{Journal of the American Statistical Association}, \bold{90}, 432, 1313-1321.

M. A. Newton and A. E. Raftery (1994) Approximate Bayesian inference by the weighted likelihood bootstrap (with discussion), \emph{Journal of 
the Royal Statistical Society}, \bold{56}, 3-48.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{logpriors_admkr}}, \code{\link[bbemkr]{logdensity_admkr}}, \code{\link[bbemkr]{mcmcrecord_admkr}}
}
\keyword{methods}


	