% $Id: formula.pdMat.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{formula.pdMat}
\title{Extract pdMat Formula}
\usage{
\method{formula}{pdMat}(object, asList)
}
\alias{formula.pdMat}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.}
 \item{asList}{logical.  Should the asList argument be applied to each
   of the components?  Never used.}
}
\description{
  This method function extracts the formula associated with a
  \code{pdMat} object, in which the column and row names are specified.
}
\value{
  if \code{object} has a \code{formula} attribute, its value is
  returned, else \code{NULL} is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\note{
  Because factors may be present in \code{formula(object)}, the
  \code{pdMat} object needs to have access to a data frame where the
  variables named in the formula can be evaluated, before it can resolve
  its row and column names from the formula.
}

\seealso{\code{\link{pdMat}}}

\examples{
pd1 <- pdSymm(~Sex*age)
formula(pd1)
}
\keyword{models}
