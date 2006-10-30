\name{rob.browne.test}
\alias{rob.browne.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{An overall goodness-of-fit of the model}
\description{
  Computing an overall goodness-of-fit of the model
}
\usage{
rob.browne.test(estcov, amat, w.par, smpl.frm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estcov}{The covariance matrix estimate}
  \item{amat}{The adjacent matrix corresponding to \code{estcov}}
  \item{w.par}{A non-negative real value, the robustness tuning parameter
    correspoinding to \code{estcov}}
  \item{smpl.frm}{A data frame, a sample data}
}
\value{
  \item{chi.val}{The test statistics}
  \item{df}{The degree of freedom}
  \item{p.val}{The probability}
}
\author{Masashi Miyamura}
\seealso{\code{\link{rob.fitConGraph}}}
\keyword{robust}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
