% $Id: getGroups.lme.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{getGroups.lme}
\title{Extract lme Object Groups}
\usage{
getGroups(object, form, level)
}
\alias{getGroups.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{form}{this argument is included to make the method function
   compatible with the generic and is ignored in this method.}
 \item{level}{an optional integer vector giving the level(s) of grouping
   to be extracted from \code{object}. Defaults to the highest or
   innermost level of grouping.}  
}
\description{
  The grouping factors corresponding to the linear mixed-effects model
  represented by \code{object} are extracted. If more than one level is
  indicated in \code{level}, the corresponding grouping factors are
  combined  into a data frame; else the selected grouping factor is
  returned as a vector. 
}
\value{
  either a data frame with columns given by the grouping factors
  indicated in \code{level}, or, when a single level is requested, a
  factor representing the selected grouping factor. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}}

\examples{
data(Pixel)
fm1 <- lme(pixel ~ day + day^2, Pixel,
  random = list(Dog = ~day, Side = ~1))
getGroups(fm1, level = 1:2)
}
\keyword{models}
