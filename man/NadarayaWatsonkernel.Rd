\name{NadarayaWatsonkernel}
\alias{NadarayaWatsonkernel}
\title{
Nadaraya-Watson kernel estimator
}
\description{
Nadaraya (1964) and Watson (1964) proposed to estimate \code{m} as a locally weighted average, using a kernel as a weighting function. 
}
\usage{
NadarayaWatsonkernel(x, y, h, gridpoint)
}
\arguments{
  \item{x}{A set of \code{x} observations.}
  \item{y}{A set of \code{y} observations.}
  \item{h}{Optimal bandwidth chosen by the user.}
  \item{gridpoint}{A set of gridpoints.}
}
\details{
\eqn{\frac{\sum^n_{i=1}K_h(x-x_i)y_i}{\sum^n_{j=1}K_h(x-x_j)}},
where \eqn{K} is a kernel function with a bandwidth \code{h}.
}
\value{
 \item{gridpoint}{A set of gridpoints.}
 \item{mh}{Density values corresponding to the set of gridpoints.}
}
\references{
M. Rosenblatt (1956) Remarks on some nonparametric estimates of a density function, \emph{The Annals of Mathematical Statistics}, \bold{27}(3), 832-837.

E. Parzen (1962) On estimation of a probability density function and mode, \emph{The Annals of Mathematical Statistics}, \bold{33}(3), 1065-1076.

E. A. Nadaraya (1964) On estimating regression, \emph{Theory of probability and its applications}, \bold{9}(1), 141-142.

G. S. Watson (1964) Smooth regression analysis, \emph{Sankhya: The Indian Journal of Statistics (Series A)}, \bold{26}(4), 359-372.
}
\author{
Han Lin Shang
}
\examples{
x = rnorm(100)
y = rnorm(100)
NadarayaWatsonkernel(x, y, h = 2, gridpoint = seq(-3, 3, length.out = 100))
}
\keyword{method}
