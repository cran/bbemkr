\name{logdensity_admkr}
\alias{logdensity_admkr}
\title{
Calculate an estimate of log posterior ordinate used in the log marginal density of Chib (1995).
}
\description{
Log marginal likelihood = Log likelihood + Log prior - Log density
}
\usage{
logdensity_admkr(tau2, cpost)
}
\arguments{
	\item{tau2}{Square of re-parameterized bandwidths and square of normal error variance}
	\item{cpost}{Simulation output of \code{tau2} obtained from the MCMC iterations}
}
\details{
It should be noted that the posterior mode or maximum likelihood estimate can be computed from the 
simulation output at least approximately, if it is easy to evaluate the log-likelihood function for each draw
in the simulation. Alternatively, one can make use of the posterior mean provided that there is no
concern that it is a low density point.
}
\value{
Value of the log density
}
\references{
S. Chib and I. Jeliazkov (2001) Marginal likelihood from the Metropolis-Hastings output, \emph{Journal of the American Statistical Association}, \bold{96}(453), 270-281.

S. Chib (1995) Marginal likelihood from the Gibbs output, \emph{Journal of the American Statistical Association}, \bold{90}(432), 1313-1321.

M. A. Newton and A. E. Raftery (1994) Approximate Bayesian inference by the weighted likelihood bootstrap (with discussion), \emph{Journal of the Royal Statistical Society. Series B}, \bold{56}(1), 3-48.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{logpriors_admkr}}, \code{\link[bbemkr]{loglikelihood_admkr}}, \code{\link[bbemkr]{mcmcrecord_admkr}}
}
\keyword{methods}
	
