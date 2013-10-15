\name{mcmcrecord_admkr}
\alias{mcmcrecord_admkr}
\title{
MCMC iterations
}
\description{
Estimated averaged bandwidths of the regressors of the kernel-form error density
}
\usage{
mcmcrecord_admkr (x, inicost, mutsizp, errorsizp, warm = 100, M = 100, prob = 0.234, 
    errorprob = 0.44, num_batch = 10, step = 10, data_x, data_y, xm, alpha = 0.05, 
      mlike = c("Chib", "Geweke", "LaplaceMetropolis", "all"))  
}
\arguments{
  \item{x}{Log of square bandwidth}
  \item{inicost}{Initial cost value}
  \item{mutsizp}{Step size of random-walk Metropolis algorithm. At each iteration, the value of \code{mutsizp} will alter depending on acceprance or rejection. As the number of iteration increases, the final acceptance probability will converge to the optimal rate, which is 0.234 for multiple parameters}
    \item{errorsizp}{Step size of random-walk Metropolis algorithm. At each iteration, the value of \code{errorsizp} will alter depending on acceprance or rejection. As the number of iteration increases, the final acceptance probability will converge to the optimal rate, which is 0.44 for single parameter}
  \item{warm}{Burn-in period}
  \item{M}{Number of MCMC iteration}
  \item{prob}{Optimal acceptance rate of random-walk Metropolis algorithm for the regression function}
  \item{errorprob}{Optimal acceptance rate of random-walk Metropolis algorithm for the error density}  
  \item{num_batch}{Number of batch samples}
  \item{step}{Recording value at a specific step, in order to achieve iid samples and eliminate correlation}
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{xm}{Values of true regression function}
  \item{alpha}{Quantile of the critical value in calculating Geweke's log marginal likelihood}
  \item{mlike}{Method for calculating log marginal likelihood}
}
\details{
Akin to the burn-in period, it determines the retained bandwidths for the regressors and the variance of 
the error density for finite samples. It also calculates the simulation inefficient factor (SIF) value,
R square, mean square error, and log marginal density by Chib (1995), Geweke (1999) and 
the Laplace Metropolis method describe in Raftery (1996).
}
\value{
    \item{sum_h}{Estimated parameters in an order of the bandwidths of the regressors, the variance
				 parameter of the error density and cost value}
    \item{h2}{Estimated parameters in an order of the square bandwidths of the regressors, the square variance parameter of the error density}
    \item{sif}{Simulation inefficient factor. The small it is, the better the method is in general}
    \item{mutsizp}{Step size of random-walk Metropolis algroithm for each iteration of \code{MCMCrecord}}
    \item{cpost}{Simulation output of square bandwidths obtained from MCMC}
    \item{ghost}{Simulation output of square bandwidths obtained from MCMC}
    \item{accept_nw}{Acceptance rate of random-walk Metropolis algorithm for the regression function}
    \item{accept_erro}{Acceptance rate of random-walk Metropolis algorithm for the kernel-form error density}    
    \item{marginalike}{Log marginal likelihood}
    \item{R2}{R square}
    \item{MSE}{Mean square error}
}
\note{
Time-consuming for large iterations. 
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density, \emph{Computational Statistics}, in press.

H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.

X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, \emph{Journal of Econometrics}, \bold{153}, 21-32.

S. Chib and I. Jeliazkov (2001) Marginal likelihood from the Metropolis-Hastings output, \emph{Journal of the American Statistical Association}, \bold{96}, 453, 270-281.

S. Chib (1995) Marginal likelihood from the Gibbs output, \emph{Journal of the American Statistical Association}, \bold{90}, 432, 1313-1321.

M. A. Newton and A. E. Raftery (1994) Approximate Bayesian inference by the weighted likelihood bootstrap (with discussion), \emph{Journal of 
the Royal Statistical Society}, \bold{56}, 3-48.

J. Geweke (1998) Using simulation methods for Bayesian econometric models: inference, development, and communication, \emph{Econometric Reviews}, \bold{18}(1), 1-73.

A. E. Raftery (1996) Hypothesis testing and model selection, in Markov Chain Monte Carlo In Practice by W. R. Gilks, S. Richardson and D. J. Spiegelhalter, Chapman and Hall, London.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{logdensity_admkr}}, \code{\link[bbemkr]{logpriors_admkr}}, \code{\link[bbemkr]{loglikelihood_admkr}}, \code{\link[bbemkr]{warmup_admkr}}
}
\keyword{methods}
