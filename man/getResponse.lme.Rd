% $Id: getResponse.lme.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{getResponse.lme}
\title{Extract lme Object Response}
\usage{
getResponse(object)
}
\alias{getResponse.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
}
\description{
  This method function extracts the response variable used in fitting
  the linear mixed-effects model corresponding to \code{object}.
}
\value{
  a vector with the response variable corresponding to the linear
  mixed-effects model represented by \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}}

\examples{
library(nlme)
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
getResponse(fm1)
}
\keyword{models}
