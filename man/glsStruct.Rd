% $Id: glsStruct.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
\name{glsStruct}
\title{Generalized Least Squares Structure}
\usage{
glsStruct(corStruct, varStruct)
}
\alias{glsStruct}
\arguments{
 \item{corStruct}{an optional \code{corStruct} object, representing a
   correlation structure. Default is \code{NULL}.}
 \item{varStruct}{an optional \code{varFunc} object, representing a
   variance function structure. Default is \code{NULL}.}
}
\description{
  A generalized least squares structure is a list of model components
  representing different sets of parameters in the linear 
  model. A \code{glsStruct}  may contain \code{corStruct} and
  \code{varFunc} objects. \code{NULL} arguments are not included in the
  \code{glsStruct} list. 
}
\value{
  a list of model variance-covariance components determining the
  parameters to be estimated for the associated linear model.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gls}}, \code{\link{corClasses}},
  \code{\link{varFunc}} } 

\examples{
gls1 <- glsStruct(corAR1(), varPower())
}
\keyword{models}
