% $Id: isInitialized.reStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{isInitialized.reStruct}
\title{Check if an reStruct Object is Initialized}
\usage{
isInitialized(object)
}
\alias{isInitialized.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
}
\description{
  Checks if all \code{pdMat} components of \code{object} have been
  initialized.
}
\value{
  a logical value indicating whether all components of \code{object} have
  been initialized. 
}

\author{Jose Pinheiro and Douglas Bates }

\seealso{\code{\link{initialize}}, \code{\link{reStruct}}}

\examples{
library(nlme)
rs1 <- reStruct(~age|Subject)
isInitialized(rs1)
}
\keyword{models}
