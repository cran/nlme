% $Id: getResponse.lme.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{getResponse.lme}
\title{Extract lme Object Response}
\usage{
\method{getResponse}{lme}(object, form)
}
\alias{getResponse.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{form}{an optional two-sided formula. Defaults to
   \code{formula(object)}.}
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
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
getResponse(fm1)
}
\keyword{models}
