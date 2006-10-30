\name{duplication}
\alias{duplication}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A duplication matrix}
\description{
  %  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
duplication(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a positive integer}
}
\value{
  A \code{n*(n+1)/2} matrix with real values.
}
\references{
  Jan R. Magnus and Heinz Neudecker, (1999) \emph{Matrix Differential Calculus
  with Applications in Statistics and Econometrics}, Second Edition, John Wiley & Sons  
}
\author{Masashi Miyamura}
\seealso{\code{\link{vech}}, \code{\link{duplication.inverse}}, \code{\link{elimination}}}
\examples{
  m <- matrix(c(1, 2, 3,
                2, 5, 6,
                3, 6, 9),
              3, 3)
  duplication(3) \%*\% vech(m) # identical to as.vector(m)
}
\keyword{utilities}% at least one, from doc/KEYWORDS

