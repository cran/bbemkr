\name{kern}
\alias{kern}
\title{
Calculate R square for goodness of fit
}
\description{
To determine the goodness of fit, we calculate the R square for leave-one-out regressors.
}
\usage{
kern(data_x, x)
}
\arguments{
  \item{data_x}{Regressors}
  \item{x}{Bandwidths of the regressors, after the MCMC iterations}
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

