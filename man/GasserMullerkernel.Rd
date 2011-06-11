\name{GasserMullerkernel}
\alias{GasserMullerkernel}
\title{
Gasser-Muller kernel estimator
}
\description{
Computes the Gasser-Muller kernel estimator 
}
\usage{
GasserMullerkernel(gridpoint, x, y, h)
}
\arguments{
  \item{gridpoint}{One gridpoint within the range of \code{x} variables. A set of gridpoints can be achieved by looping this function.}
  \item{x}{A set of \code{x} observations.}
  \item{y}{A set of \code{y} observations.}
  \item{h}{Optimal bandwidth chosen by the user.}
}
\details{
\eqn{\frac{1}{h}\sum^n_{i=1}[\int^{s_i}_{s_{i-1}}K(\frac{x-u}{h})du]y_i},
where
\eqn{s_i=\frac{x_{i-1}+x_i}{2}}.
}
\value{
Density value
}
\references{
T. Gasser, H-G. Muller and V. Mammitzsch (1985), Kernels for nonparametric curve estimation, Journal of the Royal Statistical Society (Series B), 47(2), 238-252.

Q. Li and J. S. Racine (2007), Nonparametric econometrics: Theory and Practice. Princeton University Press.

J. S. Simonoff (1996), Smoothing methods in Statistics, Springer.
}
\author{
Han Lin Shang
}
\note{
Density can be negative sometimes, although it happens sparsely on the boundary. 
}
\seealso{
\code{\link[bbemkr]{NadarayaWatsonkernel}}, \code{\link[bbemkr]{PriestleyChaokernel}}
}
\examples{
x = sort(rnorm(100))
y = seq(1, 100, length.out = 100)
gridpoint = 1
GasserMullerkernel(gridpoint, x, y, h = 0.3)
}
\keyword{method}

