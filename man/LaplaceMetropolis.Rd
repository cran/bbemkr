\name{LaplaceMetropolis}
\alias{LaplaceMetropolis}
\title{
Laplace-Metropolis estimator of log marginal likelihood
}
\description{
As pointed out by Raftery (1996), the Laplace-Metropolis estimator performs well in calculating log marginal likelihood among other methods considered.  
}
\usage{
LaplaceMetropolis(theta, data = NULL, prior_p, prior_st, 
        method = c("likelihood","L1center","median"))
}
\arguments{
  \item{theta}{MCMC output}
  \item{data}{Regressors}
  \item{prior_p}{Hyperparameter of the inverse-gamma prior}
  \item{prior_st}{Hyperparameter of the inverse-gamma prior}
  \item{method}{Computing method. \code{L1center} and \code{median} are computationally fast}
}
\details{
The idea of the Laplace-Metropolis estimator is to avoid the limitations of the Laplace method by using 
posterior simulation to estimate the quantities it needs. The Laplace method for integrals is based on
a Taylor series expansion of the real-valued function \eqn{f(u)} of the \eqn{d}-dimensional vector
\eqn{u}, and yields the approximation \eqn{P(D)\approx (2*pi)^(d/2)|A|^(1/2)P(D|\theta)P(\theta)},
where \eqn{\theta} is the posterior mode of \eqn{h(\theta)=log(P(D|\theta)P(\theta))}, \eqn{A} is minus 
the inverse Hessian of \eqn{h(\theta)} evaluated at \eqn{theta}, and \eqn{d} is the dimension of \eqn{\theta}.

The simplest way to estimate \eqn{\theta} from posterior simulation output, and probably the most accurate,
is to compute \eqn{h(\theta^(t))} for each \eqn{t=1,\dots,T} and take the value for which it is largest.
}
\value{
Log marginal likelihood
}
\references{
A. E. Raftery (1996) Hypothesis testing and model selection, in Markov Chain Monte Carlo In Practice by W. R. Gilks, S. Richardson and D. J. Spiegelhalter, Chapman and Hall, London.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{logdensity}}, \code{\link[bbemkr]{logpriors}}, \code{\link[bbemkr]{loglikelihood}}
}
\keyword{methods}

