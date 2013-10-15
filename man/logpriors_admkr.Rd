\name{logpriors_admkr}
\alias{logpriors_admkr}
\title{
Calculate the log prior used in the log marginal density of Chib (1995).
}
\description{
Log marginal likelihood = Log likelihood + Log prior - Log density
}
\usage{
logpriors_admkr(h2, data_x)
}
\arguments{
  \item{h2}{Square of re-parameterized bandwidths}
  \item{data_x}{Regressors}
}
\details{
Calculate the log prior using the estimated averaged bandwidths of the regressors, obtained from the MCMC iterations
}
\value{
Value of the log prior
}
\references{
S. Chib and I. Jeliazkov (2001) Marginal likelihood from the Metropolis-Hastings output, \emph{Journal of the American Statistical Association}, \bold{96}(453), 270-281.

S. Chib (1995) Marginal likelihood from the Gibbs output, \emph{Journal of the American Statistical Association}, \bold{90}(432), 1313-1321.

M. A. Newton and A. E. Raftery (1994) Approximate Bayesian inference by the weighted likelihood bootstrap (with discussion), \emph{Journal of the Royal Statistical Society. Series B}, \bold{56}(1), 3-48.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{logdensity_admkr}}, \code{\link[bbemkr]{loglikelihood_admkr}}, \code{\link[bbemkr]{mcmcrecord_admkr}}
}
\keyword{methods}
  