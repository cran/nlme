% $Id: getGroupsFormula.lme.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{getGroupsFormula.lme}
\title{Extract lme Object Grouping Formula}
\usage{
getGroupsFormula(object, asList)
}
\alias{getGroupsFormula.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{asList}{an optional logical value. If \code{TRUE} the returned
   value with be a list of formulas; else, if \code{FALSE} the returned
   value will be a one-sided formula. Defaults to \code{FALSE}.}
}
\description{
  The grouping formula associated with the random effects structure
  (\code{reStruct}) of \code{object} is returned either as
  a named list of one-sided formulas, or a single one-sided formula,
  depending  on the value of \code{asList}. The components of the
  returned list are ordered from outermost to innermost level and are
  named after the grouping factor expression.
}
\value{
  a one-sided formula, or a list of one-sided formulas, with the
  grouping structure associated with the random effects structure of
  \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{ \code{\link{reStruct}}, \code{\link{getGroups.lme}}}

\examples{
library(nlme)
data(Orthodont)
fm1 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
getGroupsFormula(fm1)
}
\keyword{models}
