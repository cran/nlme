% $Id: coef.gls.Rd,v 1.6 2001/10/30 20:51:14 bates Exp $
\name{coef.gls}
\title{Extract gls Coefficients}
\usage{
\method{coef}{gls}(object, allCoef, \dots)
}
\alias{coef.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{allCoef}{an optional logical value indicating whether all
   coefficients, including the non-estimable, should be
   extracted. Non-estimable coefficients are assigned
   \code{NAs}. Default is \code{FALSE}.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The estimated coefficients for the linear model represented by
  \code{object} are extracted. 
}
\value{
  a vector with the estimated coefficients for the linear model
  represented by \code{object}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gls}}} 

\examples{
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
coef(fm1)
%}
}
\keyword{models}
