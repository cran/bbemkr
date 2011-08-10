\name{mcmcrecord}
\alias{mcmcrecord}
\title{
MCMC iterations
}
\description{
Estimated averaged bandwidths of the regressors and averaged variance parameter of the normal error density
}
\usage{
mcmcrecord(x, inicost, mutsizp, warm = 100, M = 100, prob = 0.234, 
           num_batch = 10, step = 10, data_x, data_y, xm, 
           prior_p = 2, prior_st = 1)
}
\arguments{
  \item{x}{Log of square bandwidth}
  \item{inicost}{Initial cost value}
  \item{mutsizp}{Step size of random-walk Metropolis algorithm. At each iteration, the value of \code{mutsizp} will alter depending on acceprance or rejection. As the number of iteration increases, the final acceptance probability will converge to the optimal rate, which is 0.234 for multiple parameters}
  \item{warm}{Burn-in period}
  \item{M}{Number of MCMC iteration}
  \item{prob}{Optimal acceptance rate of random-walk Metropolis algorithm}
  \item{num_batch}{Number of batch samples}
  \item{step}{Recording value at a specific step, in order to achieve iid samples and eliminate correlation}
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{xm}{Values of true regression function}
  \item{prior_p}{Hyperparameter of inverse-gamma prior}
  \item{prior_st}{Hyperparameter of inverse-gamma prior}
}
\details{
Akin to the burn-in period, it determines the retained bandwidths for the regressors and the variance of 
the error density for finite samples. It also calculates the simulation inefficient factor (SIF) value,
R square, mean square error, and log marginal density by Newton and Raftery (1994) and Chib (1995).
}
\value{
    \item{sum_h}{Estimated parameters in an order of the bandwidths of the regressors, the variance
				 parameter of the error density and cost value}
    \item{h2}{Estimated parameters in an order of the square bandwidths of the regressors, the square variance parameter of the error density}
    \item{sif}{Simulation inefficient factor. The small it is, the better the method is in general}
    \item{mutsizp}{Step size of random-walk Metropolis algroithm for each iteration of \code{MCMCrecord}}
    \item{cpost}{Simulation output of square bandwidths and square normal error variance obtained from MCMC}
    \item{accept}{Acceptance rate of random-walk Metropolis algorithm}
    \item{Chib}{Log marginal likelihood of Chib (1995)}
    \item{R2}{R square}
    \item{MSE}{Mean square error}
}
\note{
Time-consuming for large iterations. 
}
\references{
X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, 
\emph{Journal of Econometrics}, \bold{153}, 21-32. 

S. Chib (1995) Marginal likelihood from the Gibbs output, \emph{Journal of the American Statistical Association}, \bold{90}, 432, 1313-1321.

M. A. Newton and A. E. Raftery (1994) Approximate Bayesian inference by the weighted likelihood bootstrap (with discussion), \emph{Journal of 
the Royal Statistical Society}, \bold{56}, 3-48.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{logdensity}}, \code{\link[bbemkr]{logpriors}}, \code{\link[bbemkr]{loglikelihood}}
}
\examples{
htm = proc.time()
x = log(nrr(data_x, FALSE)^2)
inicost = cost(x, data_x, data_y)
warmres = warmup(x, inicost, mutsizp=1.00, warm = 100, prob=0.234,
                 data_x=data_x, data_y=data_y)
mcmcres = mcmcrecord(x = warmres$x, inicost = warmres$cost, 
                     mutsizp = warmres$mutsizplast, M = 100, 
                     data_x = data_x, data_y = data_y, xm = xm)
proc.time()-htm	
}
\keyword{methods}
