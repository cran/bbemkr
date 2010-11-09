\name{bbeMCMCrecording}
\alias{bbeMCMCrecording}
\title{
MCMC iterations
}
\description{
Select the final bandwidths of the regressors and
the variance parameter of the error density}
\usage{
bbeMCMCrecording(data_x, data_y, x, costpara, num_batch = 50, 
                 M = 10000, step = 20, sizep = 1.2)
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Retained log bandwidths of the regressors, after the warmup or burnin period}
  \item{costpara}{Retained cost value, after the warmup period}
  \item{num_batch}{Number of batch samples}
  \item{M}{Number of iterations}
  \item{step}{Recording value at a specific step, in order to achieve i.i.d. samples and and eliminate or reduce correlation}
  \item{sizep}{Tuning parameter of the bandwidths}
}
\details{
Similar to the warmup period, it determines the optimal
bandwidths for the regressors and the variance of the error density
for finite sample size. It also calculates the SIF value,
R square and log marginal density by Newton and Raftery (1994) and Chib (1995).
}
\value{
 \item{accept_raterecording}{Acceptance rate of the random-walk Metropolis algorithm}
 \item{sum_h}{Selected parameters in an order of bandwidths of the regressors, variance parameter of the error density, likelihood and cost values}
 \item{std_h}{Standard deviation of the selected parameters}
 \item{batch_h}{Standard deviation of the selected parameters from different draws (equal to num_batch)}
 \item{total_sd}{Total standard deviation of the selected parameters}
 \item{sif}{Simulation inefficient factor. The small it is, the better the method is in general}
 \item{R2}{R square value for determining the goodness of fit}
 \item{data_post}{Gibbs output useful for calculating the Chib's (1995) log marginal density}
 \item{logmarginalNR}{Newton-Raftery log marginal density}
 \item{loglikelihoods}{Log likelihood for the Chib's (1995) log marginal density}
 \item{logprior}{Log prior for the Chib's (1995) log marginal density}
 \item{logdensity}{Log posterior density calculated from the Gibbs output}
 \item{logmarginalChib}{Chib's (1995) log marginal density}
}
\note{
Time-consuming for large sample size
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
\code{\link[bbemkr]{bbelogdensity}}, \code{\link[bbemkr]{bbelogpriors}}, \code{\link[bbemkr]{bbeloglikelihood}}
}
\examples{
dummy = bbewarmup(nrr(data_x), bbecost(data_x, data_y, nrr(data_x)), warm = 2)
bbeMCMCrecording(data_x, data_y, dummy$xh, dummy$cost, num_batch = 2, M = 4, step = 2)
}
\keyword{methods}
