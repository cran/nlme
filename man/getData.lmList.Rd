% $Id: getData.lmList.Rd,v 1.5 2002/03/05 14:59:39 bates Exp $
\name{getData.lmList}
\title{Extract lmList Object Data}
\usage{
getData(object)
}
\alias{getData.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  If present in the calling sequence used to produce \code{object}, the
  data frame used to fit the model is obtained. 
}
\value{
  if a \code{data} argument is present in the calling sequence that
  produced \code{object}, the corresponding data frame (with
  \code{na.action} and \code{subset} applied to it, if also present in
  the call that produced \code{object}) is returned;
  else, \code{NULL} is returned.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
getData(fm1)
}
\keyword{models}
