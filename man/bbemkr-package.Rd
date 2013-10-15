\name{bbemkr-package}
\alias{bbemkr-package}
\alias{bbemkr}
\docType{package}
\title{
Bayesian bandwidth estimation for multivariate kernel regression
}
\description{
Bayesian bandwidth estimation for Nadaraya-Watson type multivariate kernel regression with the Gaussian 
assumption of the error density and kernel-form error density
}
\details{
This package designs for estimating bandwidths used in the Nadaraya-Watson kernel regression estimator.
Assuming iid Gaussian error density that are uncorrelated to the regressors, the bandwidths are estimated
using Markov chain Monte Carlo (MCMC) method, in particular by the random-walk Metropolis algorithm.
}
\author{
Han Lin Shang and Xibin Zhang

Maintainer: Han Lin Shang <H.Shang@soton.ac.uk>
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, \emph{Computational Statistics}, in press.

H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.

X. Zhang and M. L. King and H. L. Shang (2013) A sampling algorithm for bandwidth estimation in a nonparametric regression model with a flexible error density, http://www.buseco.monash.edu.au/ebs/pubs/wpapers/2011/wp10-11.pdf.

X. Zhang and M. L. King and H. L. Shang (2013) Bayesian bandwidth selection for a nonparametric regression model with mixed types of regressors, working paper, http://www.buseco.monash.edu.au/ebs/pubs/wpapers/2013/wp13-13.pdf.

X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation,
\emph{Journal of Econometrics}, \bold{153}, 21-32.

X. Zhang and M. L. King and R. J. Hyndman (2006) A Bayesian approach to bandwidth selection for multivariate kernel density estimation, \emph{Computational Statistics and Data Analysis}, \bold{50}, 3009-3031.
}
\keyword{package}
\seealso{
\code{\link[bbemkr]{np_gibbs}}, \code{\link[bbemkr]{mcmcrecord_gaussian}}, \code{\link[bbemkr]{warmup_gaussian}},
\code{\link[bbemkr]{gibbs_admkr_nw}}, \code{\link[bbemkr]{gibbs_admkr_erro}}, \code{\link[bbemkr]{warmup_admkr}}, \code{\link[bbemkr]{mcmcrecord_admkr}}
}
