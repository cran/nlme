% $Id: update.nlsList.Rd,v 1.5 2001/01/10 19:21:05 bates Exp $
\name{update.nlsList}
\title{Update an nlsList Object}
\usage{
update(object, model, data, start, control, level, subset, na.action,
       control, pool)
}

\alias{update.nlsList}
\arguments{
  \item{object}{an object inheriting from class \code{nlsList}, representing
      a list of fitted \code{nls} objects.}
  \item{other arguments}{defined as in \code{\link{nlsList}}. See that
      function's documentation for descriptions of and default values for
      these arguments.}
}
\description{
  The non-missing arguments in the call to the \code{update.nlsList} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{nlsList} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{nlsList} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlsList}}}

\examples{
data(Soybean)
fm1 <- nlsList(weight ~ SSlogis(Time, Asym, xmid, scal) | Plot, Soybean)
fm2 <- update(fm1, start = list(Asym = 23, xmid = 57, scal = 9))
}
\keyword{models}
