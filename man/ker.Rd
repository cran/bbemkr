\name{ker}
\alias{ker}
\title{
Type of kernel function
}
\description{
For data that have infinite support, Gaussian kernel is suggested. For data that have [-1, 1] support, other types of kernel can be used.
}
\usage{
ker(u, kerntype = c("Gaussian", "Epanechnikov", "Quartic", 
                    "Triweight", "Triangular", "Uniform"))
}
\arguments{
  \item{u}{A numeric object}
  \item{kerntype}{Type of kernel function}
}
\details{
Oftentimes, we deal with numeric values of infinite support, Gaussian kernel is commonly used. 
However, Epanechnikov kernel is the optimal kernel as measured by Mean Integrated Square Error.
The difference among kernel functions is minor, but the influence of bandwidths is vital.
}
\value{
Kernel value
}
\references{
J. Fan and I. Gijbels (1996) Local Polynomial Modelling and Its Application. Chapman and Hall, London.

Q. Li and J. Racine (2007) Nonparametric Econometrics: Theory and Practice. Princeton University Press, New Jersey.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{np_gibbs}}, \code{\link[bbemkr]{gibbs_admkr_nw}}, \code{\link[bbemkr]{gibbs_admkr_erro}}
}
\keyword{methods}
