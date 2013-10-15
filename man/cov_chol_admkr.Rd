\name{cov_chol_admkr}
\alias{cov_chol_admkr}
\title{
Calculate log marginal likelihood from MCMC output
}
\description{
It is a type of candidate estimator for calculating log marginal likelihood, where the MCMC outputs are used for estimating posterior density.
}
\usage{
cov_chol_admkr(xpost, alpha, data_x, data_y)
}
\arguments{
	\item{xpost}{MCMC output}
	\item{alpha}{Quantile of the critical value in calculating Geweke's log marginal likelihood}
	\item{data_x}{Regressors}
	\item{data_y}{Response variable}	
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
\code{\link[bbemkr]{LaplaceMetropolis_admkr}}, \code{\link[bbemkr]{logdensity_admkr}}, \code{\link[bbemkr]{logpriors_admkr}}, \code{\link[bbemkr]{loglikelihood_admkr}}, \code{\link[bbemkr]{mcmcrecord_admkr}}
}
\keyword{methods}
	
