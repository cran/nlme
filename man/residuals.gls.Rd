% $Id: residuals.gls.Rd,v 1.5 2001/01/10 19:21:05 bates Exp $
\name{residuals.gls}
\title{Extract gls Residuals}
\usage{
residuals(object, type)
}
\alias{residuals.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
  \item{type}{an optional character string specifying the type of
    residuals to be used. If \code{"response"}, the "raw" residuals
    (observed - fitted) are used; else, if \code{"pearson"}, the
    standardized residuals (raw residuals divided by the corresponding
    standard errors) are used; else, if \code{"normalized"}, the
    normalized residuals (standardized residuals pre-multiplied by the
    inverse square-root factor of the estimated error correlation
    matrix) are used. Partial matching of arguments is used, so only the
    first character needs to be provided. Defaults to \code{"pearson"}.} 
}
\description{
  The residuals for the linear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the residuals for the linear model represented by
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}, \code{\link{fitted.gls}} }

\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
residuals(fm1)
}
\keyword{models}
