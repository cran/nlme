% $Id: update.lmList.Rd,v 1.5 2001/01/10 19:21:05 bates Exp $
\name{update.lmList}
\title{Update an lmList Object}
\usage{
update(object, formula, data, level, subset, na.action, control, pool)
}

\alias{update.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} fitted objects.
 }
 \item{formula}{
   a two-sided linear formula with the common model for the individuals
   \code{lm} fits.
 }
 \item{other arguments}{defined as in \code{\link{lmList}}. See
   that function's documentation for descriptions of and default values
   for these arguments.
 }
}
\description{
  The non-missing arguments in the call to the \code{update.lmList} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{lmList} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{lmList} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(Orthodont)
fm2 <- update(fm1, distance ~ I(age - 11))
summary(fm1)
summary(fm2)
}
\keyword{models}
