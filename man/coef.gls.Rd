% $Id: coef.gls.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{coef.gls}
\title{Extract gls Coefficients}
\usage{
coef(object, allCoef)
}
\alias{coef.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{allCoef}{an optional logical value indicating whether all
   coefficients, including the non-estimable, should be
   extracted. Non-estimable coefficients are assigned
   \code{NAs}. Default is \code{FALSE}.}
}
\description{
  The estimated coefficients for the linear model represented by
  \code{object} are extracted. 
}
\value{
  a vector with the estimated coefficients for the linear model
  represented by \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}} 

\examples{
library(nlme)
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
coef(fm1)
%}
}
\keyword{models}
