% $Id: fitted.gls.Rd,v 1.4.2.1 2000/12/02 21:06:08 bates Exp $
\name{fitted.gls}
\title{Extract gls Fitted Values}
\usage{
fitted(object)
}
\alias{fitted.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
}
\description{
  The fitted values for the linear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the fitted values for the linear model represented by
  \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{gls}}, \code{\link{residuals.gls}} }
\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
fitted(fm1)
}
\keyword{models}
