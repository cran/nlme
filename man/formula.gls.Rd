% $Id: formula.gls.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{formula.gls}
\title{Extract gls Object Formula}
\usage{
formula(object)
}
\alias{formula.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
}
\description{
  This method function extracts the linear model formula
  associated with \code{object}. 
}
\value{
  a two-sided linear formula specifying the linear model used to
  obtain \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}}

\examples{
library(nlme)
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
formula(fm1)
%}
}
\keyword{models}
