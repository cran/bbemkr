\name{kern}
\alias{kern}
\title{
Calculate the R square value and mean square error as measures of goodness of fit
}
\description{
To determine the goodness of fit, we calculate the R square value for leave-one-out cross validated regressors.
To determine the googness of approximating regression function, we calculate the mean square error for leave-one-out cross validated regressors.
}
\usage{
kern(h, data_x, data_y, xm)
}
\arguments{
  \item{h}{Bandwidth chosen after MCMC}
  \item{data_x}{Regressors}
  \item{data_y}{Response}
  \item{xm}{Values of true regression function}
}
\value{
Values of R square and mean square error
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{nrr}}
}
\keyword{methods}

