\name{rob.estimConGraph}
\alias{rob.estimConGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robustfied estimators for a given Gaussian concentration graph models}
\description{
  Computing robustified mean vector and covariance matrix.
}
\usage{
rob.estimConGraph(amat, w.par, smpl.frm, it.limit = 200, tol = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{An adjacent matrix}
  \item{w.par}{A non-negative real value, a robustness tuning parameter}
  \item{smpl.frm}{A data frame, a sample set}
  \item{it.limit}{a positive integer indicating the limit times of the iteration}
  \item{tol}{A small positive number idicating the tolerance used in the
  convergence tests}
}
\value{
  \item{mu}{A mean vector}
  \item{Sigma}{A covariance matrix}
  \item{weight}{A vector consisting of weight value for each observation}
  \item{w.par}{A non-negative real value, the specified robustness tuning parameter}
  \item{it}{A number of iterations}
  \item{sucess}{A logical value, indicating convergence of iteration}
}
\author{Masashi Miyamura}
\seealso{\code{\link{rob.fitConGraph}}}

\keyword{robust}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
