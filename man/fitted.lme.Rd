% $Id: fitted.lme.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{fitted.lme}
\title{Extract lme Fitted Values}
\usage{
fitted(object, level, asList)
}
\alias{fitted.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{level}{an optional integer vector giving the level(s) of grouping
   to be used in extracting the fitted values from \code{object}. Level
   values increase from outermost to innermost grouping, with
   level zero corresponding to the population fitted values. Defaults to
   the highest or innermost level of grouping.} 
 \item{asList}{an optional logical value. If \code{TRUE} and a single
   value is given in \code{level}, the returned object is a list with
   the fitted values split by groups; else the returned value is
   either a vector or a data frame, according to the length of
   \code{level}. Defaults to \code{FALSE}.}
}
\description{
  The fitted values at level \eqn{i} are obtained by adding together the
  population fitted values (based only on the fixed effects estimates)
  and the estimated contributions of the random effects to the fitted
  values at grouping levels less or equal to \eqn{i}. The resulting
  values estimate the best linear unbiased predictions (BLUPs) at level
  \eqn{i}.  
}
\value{
  If a single level of grouping is specified in \code{level}, the
  returned value is either a list with the fitted values split by groups
  (\code{asList = TRUE}) or a vector with the fitted values
  (\code{asList = FALSE}); else, when multiple grouping levels are
  specified in \code{level}, the returned object is a data frame with
  columns given by the fitted values at different levels and the
  grouping factors. 
}
\references{
  Bates, D.M. and Pinheiro, J.C. (1998) "Computational methods for
  multilevel models" available in PostScript or PDF formats at
  http://nlme.stat.wisc.edu/pub/NLME/
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}, \code{\link{residuals.lme}} }

\examples{
data(Orthodont)
fm1 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
fitted(fm1, level = 0:1)
}
\keyword{models}
