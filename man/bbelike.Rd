\name{bbelike}
\alias{bbelike}
\title{
Log likelihood constructed via Gaussian error assumption
}
\description{
Calculate the log likelihood, using the leave-one-out samples.
}
\usage{
bbelike(data_x, data_y, x, sigma, kerntype = c("Gaussian", "Epanechnikov", 
        "Quartic", "Triweight", "Triangular", "Uniform"))
}
\arguments{
  \item{data_x}{Regressors}
  \item{data_y}{Response variable}
  \item{x}{Log bandwidth of the regressors}
  \item{sigma}{Variance of the error density}
  \item{kerntype}{Type of the kernel function. By default, Gaussian kernel is used.}
}
\details{
With the assumption of i.i.d Gaussian error distribution, the likelihood can be constructed and 
maximized. 
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
