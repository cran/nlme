% $Id: gnlsStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{gnlsStruct}
\title{Generalized Nonlinear Least Squares Structure}
\usage{
gnlsStruct(corStruct, varStruct)
}
\alias{gnlsStruct}
\arguments{
 \item{corStruct}{an optional \code{corStruct} object, representing a
   correlation structure. Default is \code{NULL}.}
 \item{varStruct}{an optional \code{varFunc} object, representing a
   variance function structure. Default is \code{NULL}.}
}
\description{
  A generalized nonlinear least squares structure is a list of model
  components representing different sets of parameters in the nonlinear 
  model. A \code{gnlsStruct}  may contain \code{corStruct} and
  \code{varFunc} objects. \code{NULL} arguments are not included in the
  \code{gnlsStruct} list. 
}
\value{
  a list of model variance-covariance components determining the
  parameters to be estimated for the associated nonlinear model.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}, \code{\link{corClasses}},
  \code{\link{varFunc}} } 

\examples{
library(nlme)
gnls1 <- gnlsStruct(corAR1(), varPower())
}
\keyword{models}
