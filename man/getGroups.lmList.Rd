% $Id: getGroups.lmList.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{getGroups.lmList}
\title{Extract lmList Object Groups}
\usage{
getGroups(object)
}
\alias{getGroups.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  The grouping factor determining the partitioning of the observations
  used to produce the \code{lm} components of \code{object} is
  extracted.  
}
\value{
  a vector with the grouping factor corresponding to the \code{lm}
  components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
library(nlme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
getGroups(fm1)
}
\keyword{models}
