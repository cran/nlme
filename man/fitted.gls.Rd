% $Id: fitted.gls.Rd,v 1.8 2002/03/05 14:59:39 bates Exp $
\name{fitted.gls}
\title{Extract gls Fitted Values}
\usage{
\method{fitted}{gls}(object, \dots)
}
\alias{fitted.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The fitted values for the linear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the fitted values for the linear model represented by
  \code{object}. 
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{gls}}, \code{\link{residuals.gls}} }
\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
fitted(fm1)
}
\keyword{models}
