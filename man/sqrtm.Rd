\name{sqrtm}
\alias{sqrtm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Squared root matrix of a given symmetric matrix}
\description{
  Comupting the squared root matrix of a symmetric matrix based on the
  eigenvalue decomposition.
}
\usage{
sqrtm(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A positive definite matrix}
}
\value{
  A symmetrix matrix with dimesion \code{dim(m)}.
}
\author{Masashi Miyamura}
\examples{
m <- matrix(c(1.0, 0.8, 0.3,
              0.8, 1.0, 0.1,
              0.3, 0.1, 1.0),
            3, 3)
sqrtm(m) \%*\% sqrtm(m) # identical to m
}
\keyword{utilities}% at least one, from doc/KEYWORDS

