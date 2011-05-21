\name{bbemkr-package}
\alias{bbemkr-package}
\alias{bbemkr}
\docType{package}
\title{
Bayesian bandwidth estimation for multivariate kernel regression with Gaussian error assumption
}
\description{
Bayesian bandwidth estimation for Nadaraya-Watson type multivariate kernel regression with the Gaussian
assumption of the error density
}
\details{
This package designs for estimating bandwidths used in the Nadaraya-Watson kernel regression estimator.
Assuming i.i.d Gaussian error density that are uncorrelated to the regressors, the bandwidths are estimated
using Markov chain Monte Carlo (MCMC) method, in particular by the random-walk Metropolis algorithm.
}
\author{
Han Lin Shang and Xibin Zhang

Maintainer: Han Lin Shang <HanLin.Shang@monash.edu.au>
}
\references{
X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation,
\emph{Journal of Econometrics}, \bold{153}, 21-32.

X. Zhang and M. L. King and R. J. Hyndman (2006) A Bayesian approach to bandwidth selection for multivariate kernel density estimation,
\emph{Computational Statistics and Data Analysis}, \bold{50}, 3009-3031.
}
\keyword{package}
\seealso{
\code{\link[bbemkr]{np_gibbs}}
}
