% $Id: update.lme.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
\name{update.lme}
\title{Update an lme Object}
\usage{
update(object, fixed, data, random, correlation, weights, subset, method,
       na.action, control)
}

\alias{update.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{other arguments}{defined as in \code{\link{lme}}. See that function's
   documentation for descriptions of and default values for these arguments.}
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
}
\keyword{models}
