\name{np_gibbs}
\alias{np_gibbs}
\title{
Estimating bandwidths of the regressors 
}
\description{
Implements the random-walk Metropolis algorithm to estimate the bandwidths of the regressors
}
\usage{
np_gibbs(xh, inicost, k, mutsizp, prob, data_x, data_y, prior_p, prior_st)
}
\arguments{
  \item{xh}{Log of square bandwidths}
  \item{inicost}{Cost value}
  \item{k}{Iteration number}
  \item{mutsizp}{Step size of random-walk Metropolis algorithm}
  \item{prob}{Optimal covergence rate}
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{prior_p}{Hyperparameter used in the inverse-gamma prior}
  \item{prior_st}{Hyperparameter used in the inverse-gamma prior}
}
\details{
1) The log bandwidths of the regressors are initialized using the normal reference rule of Silverman (1986).

2) Conditioning on the variance parameter of the error density, we implement random-walk Metropolis 
   algorithm to update the bandwidths, in order to achieve the minimum cost value.	

3) The variance of the error density can be directly sampled.

4) Iterate steps 2) and 3) until the cost value is minimized.

5) Check the convergence of the parameters by examining the simulation inefficient factor (sif) value.
   The smaller the sif value is, the better convergence of the parameters is.
}
\value{
  \item{x}{}
  \item{sigma2}{Estimated variance of the normal error density}
  \item{cost}{Cost value}
  \item{accept_h}{Accept or reject. \code{accept_h=1} indicates acceptance, while \code{accept_h=0} indicates rejection.}
  \item{mutsizp}{Step size of random-walk Metropolis}
}
\references{
X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, 
\emph{Journal of Econometrics}, \bold{153}, 21-32. 
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{mcmcrecord}}, \code{\link[bbemkr]{logdensity}}, \code{\link[bbemkr]{loglikelihood}}, \code{\link[bbemkr]{logpriors}}
}
\keyword{methods}
