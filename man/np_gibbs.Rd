\name{np_gibbs}
\alias{np_gibbs}
\title{
Estimating bandwidths of the regressors and the variance of error density
}
\description{
Implements the random-walk Metropolis algorithm to estimate the bandwidths of the regressors and the
variance of the error density for finite samples
}
\usage{
np_gibbs(data_x, data_y, xh, inicost, kerntype = c("Gaussian", "Epanechnikov", "Quartic", "Triweight", "Triangular", "Uniform"), prior_p = 2, sizep, ...)
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{xh}{Log bandwidths of the regressors}
  \item{inicost}{Initial cost value}
  \item{kerntype}{Type of kernel function. By default, Gaussian kernel is used}
  \item{prior_p}{A tuning parameter in the prior of the error density}
  \item{sizep}{A tuning parameter in the random-walk Metropolis algorithm. A large value of \code{sizep}
			   decreases the acceptance rate, whereas a small value of \code{sizep} increases the 
			   acceptance rate}
  \item{\dots}{Other arguments}
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
  \item{xh}{Log bandwidths of the regressors}
  \item{sigma}{Variance of the error density}
  \item{inicost}{Initial cost value}
  \item{accept_h}{Acceptance rate of the random-walk Metropolis algorithm}
}
\references{
X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, 
\emph{Journal of Econometrics}, \bold{153}, 21-32. 
}
\author{
Han Lin Shang
}

\seealso{
\code{\link[bbemkr]{bbeMCMCrecording}}, \code{\link[bbemkr]{bbelogdensity}}
}
\examples{
inicost = bbecost(data_x, data_y, nrr(data_x))
np_gibbs(data_x, data_y, nrr(data_x), inicost, prior_p = 2, sizep = 1.2)
}
\keyword{methods}
