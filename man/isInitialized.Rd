% $Id: isInitialized.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
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
  to \code{Initialize}), by searching for components and attributes
  which are modified during initialization.
}
\value{
  a logical value indicating whether \code{object} has been
  initialized. 
}
\author{Jose Pinheiro and Douglas Bates }

\seealso{\code{\link{Initialize}}}

\examples{
pd1 <- pdDiag(~age)
isInitialized(pd1)
}
\keyword{models}
