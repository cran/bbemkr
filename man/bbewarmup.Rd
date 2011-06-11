\name{bbewarmup}
\alias{bbewarmup}
\title{
Burn-in period
}
\description{
By minimizing the cost value, the function estimates the bandwidths of the regressors and error variance
parameter of the error density for the warmup or burn-in period
}
\usage{
bbewarmup(x, costpara, kerntype = c("Gaussian", "Epanechnikov", "Quartic", "Triweight", "Triangular", "Uniform"), warm = 4, sizep = 2)
}
\arguments{
  \item{x}{Log bandwidths of the regressors}
  \item{costpara}{Initial cost value}
  \item{kerntype}{Type of kernel function. By default, Gaussian kernel is used}
  \item{warm}{Number of iterations in the warmup or burn-in period}
  \item{sizep}{A tuning parameter used in the random-walk Metropolis algorithm}
}
\value{
  \item{xh}{Log bandwidths of the regressors, obtained after the warmup period}
  \item{xhtest}{Log bandwidths of the regressors, obtained during the warmup period}
  \item{sigma}{Variance of the error density, obtained after the warmup period}
  \item{sigmatest}{Variance of the error density, obtained during the warmup period}
  \item{cost}{Cost value, obtained after the warmup period}
  \item{costtest}{Cost value, obtained during the warmup period}
  \item{accept_rate}{Acceptance rate of the random-walk Metropolis algorithm}
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{bbeMCMCrecording}}, \code{\link[bbemkr]{bbelogdensity}}
}
\keyword{methods}
