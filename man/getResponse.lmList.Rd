% $Id: getResponse.lmList.Rd,v 1.6 2002/03/05 14:59:39 bates Exp $
\name{getResponse.lmList}
\title{Extract lmList Object Response}
\usage{
\method{getResponse}{lmList}(object, form)
}
\alias{getResponse.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{form}{an optional two-sided formula. Defaults to
   \code{formula(object)}.}
}
\description{
  The response vectors from each of the \code{lm} components of
  \code{object} are extracted and combined into a single vector.  
}
\value{
  a vector with the response vectors corresponding to the \code{lm}
  components of \code{object}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
getResponse(fm1)
}
\keyword{models}
