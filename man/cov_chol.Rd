\name{cov_chol}
\alias{cov_chol}
\title{
Calculate log marginal likeliood from MCMC output
}
\description{
It is a type of candidate estimator for calculating log marginal likelihood, where the MCMC outputs are used for estimating posterior density.
}
\usage{
cov_chol(xpost, data_x, data_y, alpha, prior_p, prior_st)
}
\arguments{
  \item{xpost}{MCMC output}
  \item{data_x}{Regressor}
  \item{data_y}{Response}
  \item{alpha}{Quantile of the critical value in calculating Geweke's log marginal likelihood}
  \item{prior_p}{Hyperparameter of the inverse-gamma prior}
  \item{prior_st}{Hyperparameter of inverse-gamma prior}
}
\value{
Log marginal likelihood
}
\references{
J. Geweke (1998) Using simulation methods for Bayesian econometric models: inference, development, and communication, \emph{Econometric Reviews}, \bold{18}(1), 1-73.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{LaplaceMetropolis_gaussian}}, \code{\link[bbemkr]{logdensity_gaussian}}, \code{\link[bbemkr]{logpriors_gaussian}}, \code{\link[bbemkr]{loglikelihood_gaussian}}, \code{\link[bbemkr]{mcmcrecord_gaussian}}
}
\keyword{methods}

