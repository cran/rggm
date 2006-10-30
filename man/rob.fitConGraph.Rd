\name{rob.fitConGraph}
\alias{rob.fitConGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting a concentration graph model based on a robustified method}
\description{
  You can obtain robustified estimates of a Gaussian concentration graph model.
}
\usage{
  rob.fitConGraph(amat, w.par, smpl.frm, it.limit = 200, tol = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{An adjacent matrix}
  \item{w.par}{A non-negative real value, a robustness tuning parameter}
  \item{smpl.frm}{A data frame, a sample data}
  \item{it.limit}{a positive integer indicating the limit times of the iteration}
  \item{tol}{A small positive number idicating the tolerance used in the
  convergence tests}
}
\details{
  \item{it.limit}{When the case the number of the iteration is over
    \code{it.limit}, the fitting procedure will be turminated.}
}
\value{
  \item{mhat}{Estimates of the mean vector}
  \item{Shat}{Estimates of the covariance matrix}
  \item{w.vec}{Estimates of the weight vector}
  \item{w.par}{The specified robustness tuning parameter}
  \item{it}{The number of the iteration}
  \item{tstat}{The test statistics for the overall goodness-of-fit}
  \item{df}{The degree of freedom for the model}
  \item{p.val}{The p-value of the model}
}
\references{
  Masashi Miyamura and Yutaka Kano, 2006, Robust Gaussian graphical
  modeling, Journal of Multivariate Analysis, 97, 7, 1525--1550
}
\author{Masashi Miyamura}
\examples{
  require(mvtnorm)
  R <- matrix(c(1.0, 0.8, 0.0,
                0.8, 1.0, 0.0,
                0.0, 0.0, 1.0),
              3, 3)
  # A sample data including five outliers
  df <- as.data.frame(rbind(rmvnorm(95, rep(0, 3), R),
                            rmvnorm(5, rep(1, 3), R)))

  ug <- UG(~a*b*c)
  names(df) <- vertices(ug)

  rob.fitConGraph(ug, 0.3, df)
}
\keyword{robust}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
