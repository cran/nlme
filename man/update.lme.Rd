% $Id: update.lme.Rd,v 1.7 2001/10/30 20:51:15 bates Exp $
\name{update.lme}
\title{Update an lme Object}
\usage{
\method{update}{lme}(object, fixed, data, random, correlation, weights,
      subset, method, na.action, control, \dots)
}

\alias{update.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{other arguments}{defined as in \code{\link{lme}}. See that function's
   documentation for descriptions of and default values for these arguments.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The non-missing arguments in the call to the \code{update.lme} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{lme} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{lme} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}}

\examples{
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
fm2 <- update(fm1, distance ~ age * Sex)
summary(fm1)
summary(fm2)
}
\keyword{models}
