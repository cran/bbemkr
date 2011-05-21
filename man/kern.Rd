\name{kern}
\alias{kern}
\title{
Calculate the R square value as a measure of goodness of fit
}
\description{
To determine the goodness of fit, we calculate the R square value for leave-one-out regressors.
}
\usage{
kern(data_x, x, kerntype = c("Gaussian", "Epanechnikov", "Quartic", 
     "Triweight", "Triangular", "Uniform"))
}
\arguments{
  \item{data_x}{Regressors}
  \item{x}{Bandwidths of the regressors, after the MCMC iterations}
  \item{kerntype}{Type of kernel function. By default, Gaussian kernel is used}
}
\value{
Value of R square
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{nrr}}
}
\keyword{methods}

