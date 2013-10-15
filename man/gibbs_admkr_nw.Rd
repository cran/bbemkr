\name{gibbs_admkr_nw}
\alias{gibbs_admkr_nw}
\title{
Estimating bandwidths of the regressors
}
\description{
Implements the random-walk Metropolis algorithm to estimate the bandwidths of the regressors
}
\usage{
gibbs_admkr_nw(xh, inicost, k, mutsizp, prob, data_x, data_y)
}
\arguments{
	\item{xh}{Log of square bandwidths in the regression function}
	\item{inicost}{Cost value}
	\item{k}{Iteration number}
	\item{mutsizp}{Step size of random-walk Metropolis algorithm}
	\item{prob}{Optimal convergence rate for drawing single or multiple parameters}
	\item{data_x}{Regressors}
	 \item{data_y}{Response variable}
}
\details{
1) The log bandwidths of the regressors are initialized using the normal reference rule of Silverman (1986).

2) Conditioning on the variance parameter of the error density, we implement random-walk Metropolis 
   algorithm to update the bandwidths, in order to achieve the minimum cost value.	

3) The bandwidth of the kernel-form error density can be directly sampled.

4) Iterate steps 2) and 3) until the cost value is minimized.

5) Check the convergence of the parameters by examining the simulation inefficient factor (sif) value.
   The smaller the sif value is, the better convergence of the parameters is.
}
\value{
	\item{x}{Estimated bandwidths of the regression function}
	\item{cost}{Cost value, that is negative of log posterior}
	\item{accept_h}{Accept or reject. \code{accept_h=1} indicates acceptance, while \code{accept_h=0} indicates rejection.}
  	\item{mutsizp}{Step size of the random-walk Metropolis algorithm}
}
\references{
X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, \emph{Journal of Econometrics}, \bold{153}, 21-32. 

B. W. Silverman (1986) Density Estimation for Statistics and Data Analysis. Chapman and Hall, New York.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{mcmcrecord_admkr}}, \code{\link[bbemkr]{logdensity_admkr}}, \code{\link[bbemkr]{loglikelihood_admkr}}, \code{\link[bbemkr]{logpriors_admkr}}, \code{\link[bbemkr]{gibbs_admkr_erro}}
}
\keyword{methods}

