% $Id: initialize.corStruct.Rd,v 1.3 2000/04/17 19:03:39 bates Exp $
\name{initialize.corStruct}
\title{Initialize corStruct Object}
\usage{
initialize(object, data, ...)
}
\alias{initialize.corStruct}
\alias{initialize.corAR1}
\alias{initialize.corARMA}
\alias{initialize.corCAR1}
\alias{initialize.corCompSymm}
\alias{initialize.corHF}
\alias{initialize.corIdent}
\alias{initialize.corLin}
\alias{initialize.corNatural}
\alias{initialize.corSpatial}
\alias{initialize.corSpher}
\alias{initialize.corSymm}
\arguments{
 \item{object}{an object inheriting from class \code{corStruct}
   representing a correlation structure.}
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{formula(object)}.}
 \item{\dots}{this argument is included to make this method compatible
   with the generic.}
}
\description{
  This method initializes \code{object} by evaluating its associated
  covariate(s) and grouping factor, if any is present, in \code{data},
  calculating various dimensions and constants used by optimization
  algorithms involving \code{corStruct} objects (see the appropriate
  \code{Dim} method documentation), and assigning initial values for
  the coefficients in \code{object}, if none were present.
}
\value{
  an initialized object with the same class as \code{object}
  representing a correlation structure.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{Dim.corStruct}}}

\examples{
data(Orthodont)
cs1 <- corAR1(form = ~ 1 | Subject)
cs1 <- initialize(cs1, data = Orthodont)
}
\keyword{models}
