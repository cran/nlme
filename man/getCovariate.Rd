% $Id: getCovariate.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{getCovariate}
\title{Extract Covariate from an Object}
\usage{
getCovariate(object, form, data)
}
\alias{getCovariate}
\arguments{
 \item{object}{any object with a \code{covariate} component}
 \item{form}{an optional one-sided formula specifying the covariate(s)
   to be extracted. Defaults to \code{formula(object)}.}
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{form}. }
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{corStruct}, \code{corSpatial},
  \code{data.frame}, and \code{varFunc}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{
    \code{\link{getCovariate.corStruct}},
    \code{\link{getCovariate.data.frame}},
    \code{\link{getCovariate.varFunc}},
    \code{\link{getCovariateFormula}}
}
\examples{
## see the method function documentation
}
\keyword{models}
