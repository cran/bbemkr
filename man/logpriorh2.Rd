\name{logpriorh2}
\alias{logpriorh2}
\title{
Prior of square bandwidths
}
\description{
Prior of square bandwidths
}
\usage{
logpriorh2(h2, prior_alpha = 1, prior_beta = 0.05)
}
\arguments{
  \item{h2}{Square bandwidths}
  \item{prior_alpha}{Hyperparameter of inverse-gamma prior}
  \item{prior_beta}{Hyperparameter of inverse-gamma prior}
}
\details{
Prior choice of bandwidths
}
\value{
Prior value
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[bbemkr]{mcmcrecord_gaussian}}, \code{\link[bbemkr]{logdensity_gaussian}}, \code{\link[bbemkr]{loglikelihood_gaussian}}, \code{\link[bbemkr]{logpriors_gaussian}}
}
\keyword{methods}

