\name{nrr}
\alias{nrr}
\title{
Normal reference rule
}
\description{
Rule of thumb for estimating bandwidths.
}
\usage{
nrr(data_x, logband = TRUE)
}
\arguments{
  \item{data_x}{Regressors}
  \item{logband}{When \code{logband=TRUE}, log bandwidths are given. When \code{logband=FALSE},
				 bandwidths are given.}
}
\value{
Bandwidths.
}
\references{
X. Zhang and R. D. Brooks and M. L. King (2009) A Bayesian approach to bandwidth selection for multivariate kernel regression with an application to state-price density estimation, 
\emph{Journal of Econometrics}, \bold{153}, 21-32. 

D. W. Scott (1992) Multivariate density estimation: theory, practice, and visualization, John Wiley and Sons, New York.

A. W. Bowman and A. Azzalini (1997) Applied smoothing techniques for data analysis, Oxford University Press, London.
}
\author{
Han Lin Shang
}
\examples{
nrr(data_x)
}
\keyword{methods}

