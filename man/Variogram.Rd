% $Id: Variogram.Rd,v 1.5.2.1 2003/04/17 22:28:45 bates Exp $
\name{Variogram}
\title{Calculate Semi-variogram}
\usage{
Variogram(object, distance, \dots)
}
\alias{Variogram}
\arguments{
  \item{object}{a numeric vector with the values to be used for
    calculating the semi-variogram, usually a residual vector from a
    fitted model.}
  \item{distance}{a numeric vector with the pairwise distances
    corresponding to the elements of \code{object}. The order of the
    elements in \code{distance} must correspond to the pairs
    \code{(1,2), (1,3), \dots, (n-1,n)}, with \code{n} representing the
    length of \code{object}, and must have length \code{n(n-1)/2}.}
  \item{\dots}{some methods for this generic function require additional
    arguments.}
}
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

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{Variogram.default}},\code{\link{Variogram.gls}},
  \code{\link{Variogram.lme}}, \code{\link{plot.Variogram}}}
\examples{
## see the method function documentation
}
\keyword{models}
