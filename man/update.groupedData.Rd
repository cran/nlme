% $Id: update.groupedData.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{update.groupedData}
\title{Update a groupedData Object}
\usage{
\method{update}{groupedData}(object, formula, data, order.groups, FUN,
outer, inner, labels, units, \dots)
}
\alias{update.groupedData}
\arguments{
 \item{object}{an object inheriting from class \code{groupedData}.}
 \item{other arguments}{defined as in \code{\link{groupedData}}. See
   that function's documentation for descriptions of and default values
   for these arguments.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The non-missing arguments in the call to the \code{update.groupedData}
  method replace the corresponding arguments in the original call used to
  produce \code{object} and \code{groupedData} is used with the modified
  call to produce an updated fitted object.
}
\value{
  an updated \code{groupedData} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{groupedData}}}

\examples{
data(Orthodont)
Orthodont2 <- update(Orthodont, FUN = mean)
}
\keyword{models}
