\name{bbewarmup}
\alias{bbewarmup}
\title{
Burn-in period
}
\description{
To minimize the cost value, it estimates the bandwidths of the regressors and variance parameter of the error density for the warmup or burn-in period.
}
\usage{
bbewarmup(x, costpara, warm = 4, sizep = 2)
}
\arguments{
  \item{x}{Log bandwidths of the regressors}
  \item{costpara}{Initial cost value}
  \item{warm}{Number of iterations in the warmup or burn-in period}
  \item{sizep}{Tuning parameter used in the random-walk Metropolis algorithm}
}
\value{
\item{xh}{Log bandwidths of the regressors, after the warmup period}
\item{xhtest}{Log bandwidths of the regressors, during the warmup period}
\item{sigma}{Variance of the error density, after the warmup period}
\item{sigmatest}{Variance of the error density, during the warmup period}
\item{cost}{Cost value, after the warmup period}
\item{costtest}{Cost value, during the warmup period}
\item{accept_rate}{Acceptance rate of the random-walk Metropolis algorithm}
}
\author{
Han Lin Shang
}

\seealso{
\code{\link[bbemkr]{bbeMCMCrecording}}, \code{\link[bbemkr]{bbelogdensity}}
}
\keyword{methods}
