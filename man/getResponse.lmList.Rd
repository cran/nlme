% $Id: getResponse.lmList.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{getResponse.lmList}
\title{Extract lmList Object Response}
\usage{
getResponse(object)
}
\alias{getResponse.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  The response vectors from each of the \code{lm} components of
  \code{object} are extracted and combined into a single vector.  
}
\value{
  a vector with the response vectors corresponding to the \code{lm}
  components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
library(nlme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
getResponse(fm1)
}
\keyword{models}
