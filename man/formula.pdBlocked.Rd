% $Id: formula.pdBlocked.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{formula.pdBlocked}
\title{Extract pdBlocked Formula}
\usage{
formula(object, asList)
}
\alias{formula.pdBlocked}
\arguments{
 \item{object}{an object inheriting from class \code{pdBlocked},
   representing a positive definite block diagonal matrix.} 
 \item{asList}{an optional logical value. If \code{TRUE}, a list with
   the formulas for the individual block diagonal elements of
   \code{object} is returned; else, if \code{FALSE}, a one-sided formula
   combining all individual formulas is returned. Defaults to
   \code{FALSE}.} 
}
\description{
  The \code{formula} attributes of the \code{pdMat} elements of
  \code{object} are extracted and returned as a list, in case
  \code{asList=TRUE}, or converted to a single one-sided formula when
  \code{asList=FALSE}. If the \code{pdMat} elements do not have a
  \code{formula} attribute, a \code{NULL} value is returned.
}
\value{
  a list of one-sided formulas, or a single one-sided formula, or
  \code{NULL}. 
}

\author{Jose Pinheiro and Douglas Bates}


\seealso{\code{\link{pdBlocked}}, \code{\link{pdMat}}}

\examples{
pd1 <- pdBlocked(list(~ age, ~ Sex - 1))
formula(pd1)
formula(pd1, asList = TRUE)
}
\keyword{models}
