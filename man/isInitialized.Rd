% $Id: isInitialized.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{isInitialized}
\title{Check if Object is Initialized}
\usage{
isInitialized(object)
}
\alias{isInitialized}
\alias{isInitialized.pdMat}
\alias{isInitialized.pdBlocked}
\arguments{
 \item{object}{any object requiring initialization.}
}
\description{
  Checks if \code{object} has been initialized (generally through a call
  to \code{initialize}), by searching for components and attributes
  which are modified during initialization.
}
\value{
  a logical value indicating whether \code{object} has been
  initialized. 
}
\author{Jose Pinheiro and Douglas Bates }

\seealso{\code{\link{initialize}}}

\examples{
library(nlme)
pd1 <- pdDiag(~age)
isInitialized(pd1)
}
\keyword{models}
