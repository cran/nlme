% $Id: update.gnls.Rd,v 1.5 2000/07/03 18:22:48 bates Exp $
\name{update.gnls}
\title{Update a gnls Object}
\usage{
update(object, model, data, params, start, correlation, weights, subset,
       na.action, naPattern, control, verbose) 
}

\alias{update.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
 \item{other arguments}{defined as in \code{\link{gnls}}. See
   that function's documentation for descriptions of and default values
   for these arguments.}
}
\description{
  The non-missing arguments in the call to the \code{update.gnls} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{gnls} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{gnls} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}}

\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
\dontrun{
fm2 <- update(fm1, correlation = corAR1())
}
}
\keyword{models}
