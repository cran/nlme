% $Id: getGroupsFormula.reStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{getGroupsFormula.reStruct}
\title{Extract reStruct Grouping Formula}
\usage{
getGroupsFormula(object, asList)
}
\alias{getGroupsFormula.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{asList}{an optional logical value. If \code{TRUE} the returned
   value with be a list of formulas; else, if \code{FALSE} the returned
   value will be a one-sided formula. Defaults to \code{FALSE}.}
}
\description{
  The names of the \code{object} components are used to construct a
  one-sided formula, or a named list of formulas, depending on the value
  of \code{asList}.  The components of the returned list are ordered
  from outermost to innermost level.
}
\value{
  a one-sided formula, or a list of one-sided formulas, with the
  grouping structure associated with \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{reStruct}}, \code{\link{getGroups}}}

\examples{
library(nlme)
rs1 <- reStruct(list(A = pdDiag(diag(2), ~age), B = ~1))
getGroupsFormula(rs1)
}
\keyword{models}
