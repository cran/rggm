\name{vech}
\alias{vech}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The half vectororized elements of a given matrix}
\description{
%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vech(m, diag = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A symmetric matrix}
  \item{diag}{a logical value. If TRUE, diag elements will be included.}
}
\value{
  A \code{p*(p+1)/2} x \code{1} matrix (a vector) obtained by stacking
  all columns of a given matrix
}
\references{
  Jan R. Magnus and Heinz Neudecker, (1999) \emph{Matrix Differential Calculus
  with Applications in Statistics and Econometrics}, Second Edition, John Wiley & Sons\cr
}
\author{Masashi Miyamura}
\seealso{\code{\link{duplication}}}
\keyword{utilities}% at least one, from doc/KEYWORDS

