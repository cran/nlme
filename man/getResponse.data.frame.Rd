% $Id: getResponse.data.frame.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{getResponse.data.frame}
\title{Extract Response from a Data Frame}
\usage{
\method{getResponse}{data.frame}(object, form)
}
\alias{getResponse.data.frame}
\arguments{
 \item{object}{an object inheriting from class \code{data.frame}.}
 \item{form}{an optional formula specifying the response to be
   evaluated in \code{object}. Defaults to \code{formula(object)}.}
}
\description{
  The left hand side of \code{form} is evaluated in \code{object}.  
}
\value{
  the value of the left hand side of \code{form} evaluated in
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getResponseFormula}}}

\examples{
data(Orthodont)
getResponse(Orthodont)
}
\keyword{models}
