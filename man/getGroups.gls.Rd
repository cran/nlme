% $Id: getGroups.gls.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
\name{getGroups.gls}
\title{Extract gls Object Groups}
\usage{
getGroups(object)
}
\alias{getGroups.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
}
\description{
  If present, the grouping factor associated to the correlation
  structure for the linear model represented by \code{object} is extracted. 
}
\value{
  if the linear model represented by \code{object} incorporates a
  correlation structure and the corresponding \code{corStruct} object
  has a grouping factor, a vector with the group values is returned;
  else, \code{NULL} is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}, \code{\link{corClasses}}}

\examples{
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
getGroups(fm1)
%}
}
\keyword{models}
