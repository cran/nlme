% $Id: fitted.gls.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
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
library(nlme)
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
fitted(fm1)
%}
}
\keyword{models}
