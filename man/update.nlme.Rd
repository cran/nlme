% $Id: update.nlme.Rd,v 1.9 2001/06/18 21:16:50 bates Exp $
\name{update.nlme}
\title{Update an nlme Object}
\usage{
\method{update}{nlme}(object, model, data, fixed, random, groups, start,
       correlation, weights, subset, method, na.action, naPattern,
       control, verbose)
}
\alias{update.nlme}
\arguments{
 \item{object}{an object inheriting from class \code{nlme}, representing
   a fitted nonlinear mixed-effects model.}
 \item{other arguments}{defined as in \code{\link{nlme}}. See that
     function's documentation for descriptions of and default values
     for these arguments.}
}
\description{
  The non-missing arguments in the call to the \code{update.nlme} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{nlme} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{nlme} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlme}}}

\examples{
data(Loblolly)
fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
            data = Loblolly,
            fixed = Asym + R0 + lrc ~ 1,
            random = Asym ~ 1,
            start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
fm2 <- update(fm1, random = pdDiag(Asym + lrc ~ 1))
summary(fm1)
summary(fm2)
}
\keyword{models}
