% $Id: collapse.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{collapse}
\title{Collapse According to Groups}
\usage{
collapse(object, ...)
}
\alias{collapse}
\arguments{
 \item{object}{an object to be collapsed, usually a data frame.}
 \item{\dots}{some methods for the generic may require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Currently, only a \code{groupedData}
  method is available. 
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{collapse.groupedData}}}

\examples{
## see the method function documentation
}
\keyword{models}
