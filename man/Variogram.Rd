% $Id: Variogram.Rd,v 1.3 2000/07/03 18:22:46 bates Exp $
\name{Variogram}
\title{Calculate Semi-variogram}
\usage{
Variogram(object, distance, ...)
}
\alias{Variogram}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{default}, \code{gls} and \code{lme}. See
  the appropriate method documentation for a description of the
  arguments. 
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\references{
  Cressie, N.A.C. (1993), "Statistics for Spatial Data", J. Wiley & Sons.
}  

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{Variogram.default}},\code{\link{Variogram.gls}},
  \code{\link{Variogram.lme}}, \code{\link{plot.Variogram}}}
\examples{
## see the method function documentation
}
\keyword{models}
