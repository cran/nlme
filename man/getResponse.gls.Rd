% $Id: getResponse.gls.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{getResponse.gls}
\title{Extract gls Object Response}
\usage{
\method{getResponse}{gls}(object, form)
}
\alias{getResponse.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{form}{an optional two-sided formula. Defaults to
   \code{formula(object)}.}
}
\description{
  This method function extracts the response variable used in fitting
  the linear model corresponding to \code{object}.
}
\value{
  a vector with the response variable corresponding to the linear
  model represented by \code{object}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gls}}}

\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
getResponse(fm1)
}
\keyword{models}
