% $Id: nlmeStruct.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{nlmeStruct}
\title{Nonlinear Mixed-Effects Structure}
\usage{
nlmeStruct(reStruct, corStruct, varStruct)
}
\alias{nlmeStruct}
\arguments{
 \item{reStruct}{a \code{reStruct} representing a random effects
   structure.}
 \item{corStruct}{an optional \code{corStruct} object, representing a
   correlation structure. Default is \code{NULL}.}
 \item{varStruct}{an optional \code{varFunc} object, representing a
   variance function structure. Default is \code{NULL}.}
}
\description{
  A nonlinear mixed-effects structure is a list of model components
  representing different sets of parameters in the nonlinear mixed-effects
  model. An \code{nlmeStruct} list must contain at least a
  \code{reStruct} object, but may also contain \code{corStruct} and
  \code{varFunc} objects. \code{NULL} arguments are not included in the
  \code{nlmeStruct} list. 
}
\value{
  a list of model components determining the parameters to be estimated
  for the associated nonlinear mixed-effects model.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlme}}, \code{\link{reStruct}},
  \code{\link{corClasses}}, \code{\link{varFunc}} }

\examples{
nlms1 <- nlmeStruct(reStruct(~age), corAR1(), varPower())
}
\keyword{models}
