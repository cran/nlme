% $Id: update.groupedData.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{update.groupedData}
\title{Update a groupedData Object}
\usage{
update(object, formula, data, order.groups, FUN, outer, inner, labels, units)
}
\alias{update.groupedData}
\arguments{
 \item{object}{an object inheriting from class \code{groupedData}.}
 \item{other arguments}{defined as in \code{\link{groupedData}}. See
   that function's documentation for descriptions of and default values
   for these arguments.}
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
library(nlme)
data(Orthodont)
Orthodont2 <- update(Orthodont, FUN = mean)
}
\keyword{models}
