\name{cost_admkr}
\alias{cost_admkr}
\title{
Negative of log posterior associated with the bandwidths
}
\description{
Calculates the negative of log posterior, using the leave-one-out cross validated samples.
}
\usage{
cost_admkr(x, data_x, data_y)
}
\arguments{
	\item{x}{Log of square bandwidths}
	\item{data_x}{Regressors}
	\item{data_y}{Response variable}
}
\details{
Bandwidth can be re-parameterized by a constant time optimal convergence rate, that is, \eqn{h = c*n^{rate}}.
}
\value{
Value of the cost function
}
\references{
H. L. Shang (2013) Bayesian bandwidth estimation for a nonparametric functional regression model with unknown error density, \emph{Computational Statistics and Data Analysis}, \bold{67}, 185-198.

X. Zhang, M. L. King and H. L. Shang (2013). A sampling algorithm for bandwidth estimation in a nonparametric regression model with a flexible error density. Working paper, http://users.monash.edu.au/~xzhang/zhang.king.shang.2013.pdf

X. Zhang, M. L. King and H. L. Shang (2013). Bayesian bandwidth selection for a nonparametric regression model with mixed types of regressors. Working paper, http://www.buseco.monash.edu.au/ebs/pubs/wpapers/2013/wp13-13.pdf

X. Zhang and M. L. King (2013). Gaussian kernel GARCH models. Working paper, http://users.monash.edu.au/~xzhang/zhang.king.2013.rev.pdf
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{gibbs_admkr_nw}}, \code{\link[bbemkr]{gibbs_admkr_erro}}
}
\examples{
x = log(c(nrr(data_x, FALSE),2)^2)
inicost = cost_admkr(x, data_x = data_x, data_y = data_ynorm)
}
