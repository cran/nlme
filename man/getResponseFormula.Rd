\name{getResponseFormula}
\title{Extract Formula Specifying Response Variable}
\usage{
getResponseFormula(object)
}
\alias{getResponseFormula}
\arguments{
 \item{object}{any object from which a formula can be extracted.}
}
\description{
  The left hand side of \code{formula{object}} is returned as a
  one-sided formula. 
}
\value{
  a one-sided formula with the response variable associated with
  \code{formula{object}}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getResponse}}}

\examples{
library(nlme)
getResponseFormula(y ~ x | g)
}
\keyword{models}

