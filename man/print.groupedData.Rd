% $Id: print.groupedData.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{print.groupedData}
\title{Print a groupedData Object}
\usage{
print(x, ...)
}
\alias{print.groupedData}
\arguments{
 \item{x}{an object inheriting from class \code{groupedData}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  Prints the display formula and the data frame associated with
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{groupedData}}}

\examples{
library(nlme)
data(Orthodont)
print(Orthodont)
}
\keyword{models}
