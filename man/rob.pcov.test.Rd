\name{rob.pcov.test}
\alias{rob.pcov.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Testing zero partial covariances}
\description{
  Note that issues with regard to multiple testing are not concerned.
}
\usage{
rob.pcov.test(estcov, amat, w.par, smpl.frm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estcov}{The estimated covariance matrix}
  \item{amat}{An adjacent matrix specifying one model}
  \item{w.par}{A non-negative real value, the robustness tuning parameter}
  \item{smpl.frm}{A data frame, a sample data}
}
\value{
  A list including all pairs of vertices. Each elements consists of
  \item{vs}{A character pair indicating variable names}
  \item{val}{A real value, the value of the partial covariance}
  \item{aval}{A rael value, an asymptotic value of the partial covariance}
  \item{tsta}{A real value, a test statistics}
  \item{pval}{A real value, a p-value}
}
\author{Masashi Miyamura}
\keyword{robust}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
\keyword{utilities}% __ONLY ONE__ keyword per line
