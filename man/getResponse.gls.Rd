% $Id: getResponse.gls.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{getResponse.gls}
\title{Extract gls Object Response}
\usage{
getResponse(object)
}
\alias{getResponse.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
}
\description{
  This method function extracts the response variable used in fitting
  the linear model corresponding to \code{object}.
}
\value{
  a vector with the response variable corresponding to the linear
  model represented by \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}}

\examples{
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
getResponse(fm1)
%}
}
\keyword{models}
