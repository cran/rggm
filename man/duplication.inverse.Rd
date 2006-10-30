\name{duplication.inverse}
\alias{duplication.inverse}
\alias{elimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The inverse of a duplication matrix}
\description{
  Return the inverse of a duplication matrix, so called an elimination matrix
}
\usage{
duplication.inverse(n)
elimination(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a positive integer, indicating a dimmension}
}
\value{
  An \code{n*(n+1)/2} rows \code{n} columns elimination matrix.
}
\references{
  Jan R. Magnus and Heinz Neudecker, (1999) Matrix Differential Calculus
  with Applications in Statistics and Econometrics, Second Edition, John Wiley & Sons\cr
}
\author{Masashi Miyamura}
\seealso{\code{\link{vech}}, \code{\link{duplication}}, \code{\link{elimination}}}
\examples{
  m <- matrix(c(1, 2, 3,
                2, 5, 6,
                3, 6, 9),
              3, 3)
  elimination(3) \%*\% as.vector(m) # identical to vech(m)

  duplication.inverse(3) \%*\% duplication(3) # An identity matrix
  elimination(3) \%*\% duplication(3) # The same thing
}
\keyword{utilities}% at least one, from doc/KEYWORDS
