\name{bbelike}
\alias{bbelike}
\title{
Log likelihood 
}
\description{
Calculate the log likelihood, using the leave-one-out samples.
}
\usage{
bbelike(data_x, data_y, x, sigma)
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Log bandwidth of the regressors}
  \item{sigma}{Variance of the error density}
}
\value{
Log likelihood value
}
\references{
X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, 
\emph{Journal of Econometrics}, \bold{153}, 21-32. 
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{bbeMCMCrecording}}
}
\keyword{methods}
