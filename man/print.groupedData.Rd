% $Id: print.groupedData.Rd,v 1.2 2000/04/17 19:03:40 bates Exp $
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
data(Orthodont)
print(Orthodont)
}
\keyword{models}
