% $Id: getGroupsFormula.gls.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{getGroupsFormula.gls}
\title{Extract gls Object Grouping Formula}
\usage{
getGroupsFormula(object, asList)
}
\alias{getGroupsFormula.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{asList}{an optional logical value. If \code{TRUE} the returned
   value with be a list of formulas; else, if \code{FALSE} the returned
   value will be a one-sided formula. Defaults to \code{FALSE}.}
}
\description{
  If present, the grouping formula associated with the correlation
  structure (\code{corStruct}) of \code{object} is returned either as
  a named list with a single one-sided formula, or a single one-sided
  formula, depending on the value of \code{asList}. If \code{object}
  does not include a correlation structure, or if the correlation
  structure does not include groups, \code{NULL} is returned.
}
\value{
  if a correlation structure with groups is included in \code{object}, a
  one-sided formula, or a list with a single one-sided formula, with the
  corresponding grouping structure, else \code{NULL}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{ \code{\link{corClasses}}, \code{\link{getGroups.gls}}}

\examples{
library(nlme)
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
getGroupsFormula(fm1)
%}
}
\keyword{models}
