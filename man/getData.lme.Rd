% $Id: getData.lme.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{getData.lme}
\title{Extract lme Object Data}
\usage{
\method{getData}{lme}(object)
}
\alias{getData.lme}
\alias{getData.nlme}
\alias{getData.nls}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a linear mixed-effects fitted model.}
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

\seealso{\code{\link{lme}}}

\examples{
data(Ovary)
fm1 <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary,
           random = ~ sin(2*pi*Time))
getData(fm1)
}
\keyword{models}
