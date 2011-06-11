\name{PriestleyChaokernel}
\alias{PriestleyChaokernel}
\title{
Priestley-Chao kernel estimator
}
\description{
Computes Priestley-Chao kernel estimator
}
\usage{
PriestleyChaokernel(x, y, h, gridpoint)
}
\arguments{
  \item{x}{A set of \code{x} observations.}
  \item{y}{A set of \code{y} observations.}
  \item{h}{Optimal bandwidth chosen by the user.}
  \item{gridpoint}{A set of gridpoints.}
}
\details{
\eqn{\frac{1}{h}\sum^n_{i=1}(x_i-x_{i-1})K(\frac{x-x_i}{h})y_i}
}
\value{
  \item{gridpoint}{A set of gridpoints.}
  \item{mh}{Density values corresponding to the set of gridpoints.}
}
\references{
M. B. Priestley and M. T. Chao (1972) Nonparametric function fitting, Journal of Royal Statistical Society (Series B), 34(3), 385-392.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{GasserMullerkernel}}, \code{\link[bbemkr]{NadarayaWatsonkernel}}
}
\examples{
x = rnorm(100)
y = rnorm(100)
PriestleyChaokernel(x, y, h = 2, gridpoint = seq(-3,3,length.out=100))
}
\keyword{method}
