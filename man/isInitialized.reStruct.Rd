% $Id: isInitialized.reStruct.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{isInitialized.reStruct}
\title{Check if an reStruct Object is Initialized}
\usage{
\method{isInitialized}{reStruct}(object)
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
rs1 <- reStruct(~age|Subject)
isInitialized(rs1)
}
\keyword{models}
