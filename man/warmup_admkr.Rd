\name{warmup_admkr}
\alias{warmup_admkr}
\title{
Burn-in period
}
\description{
By minimising the cost value, the function estimates the bandwidths of the regressors and kernel-form error density for the burn-in period
}
\usage{
warmup_admkr(x, inicost, mutsizp, errorsizp, warm = 100, prob = 0.234, 
    errorprob = 0.44, data_x, data_y) 
}
\arguments{
	\item{x}{Log of square bandwidths}
	\item{inicost}{Cost value}
	\item{mutsizp}{Step size of random-walk Metropolis algorithm for the regressors}
	\item{errorsizp}{Step size of random-walk Metropolis algorithm for the kernel-form error density}
	\item{warm}{Number of burn-in iterations}
	\item{prob}{Optimal covergence rate of random-walk Metropolis algorithm for the regressors}
	\item{errorprob}{Optimal covergence rate of random-walk Metropolis algorithm for the kernel-form error density}
	\item{data_x}{Regressors}
	\item{data_y}{Response variable}
}
\value{
	\item{x}{Log of square bandwidths}
	\item{cost}{Cost value}
	\item{mutsizp}{Step size of random-walk Metropolis algorithm for the regressors}
	\item{errorsizp}{Step size of random-walk Metropolis algorithm for the kernel-form error density}
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{mcmcrecord_admkr}}, \code{\link[bbemkr]{logdensity_admkr}}, \code{\link[bbemkr]{loglikelihood_admkr}}, \code{\link[bbemkr]{logpriors_admkr}}
}
\keyword{methods}


	