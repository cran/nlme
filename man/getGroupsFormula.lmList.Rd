% $Id: getGroupsFormula.lmList.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
\name{getGroupsFormula.lmList}
\title{Extract lmList Object Grouping Formula}
\usage{
getGroupsFormula(object, asList)
}
\alias{getGroupsFormula.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{asList}{an optional logical value. If \code{TRUE} the returned
   value with be a list of formulas; else, if \code{FALSE} the returned
   value will be a one-sided formula. Defaults to \code{FALSE}.
 }
}
\description{
  A formula representing the grouping factor determining the
  partitioning of the observations used to produce the \code{lm}
  components of \code{object} is obtained and returned as a list with a
  single component, or as a one-sided formula.  
}
\value{
  a one-sided formula, or a list with a single one-sided formula,
  representing the grouping factor corresponding to the \code{lm}
  components of \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{ \code{\link{lmList}}, \code{\link{getGroups.lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
getGroupsFormula(fm1)
}
\keyword{models}
