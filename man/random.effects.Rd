% $Id: random.effects.Rd,v 1.4.4.1 2003/04/17 22:28:45 bates Exp $
\name{random.effects}
\title{Extract Random Effects}
\usage{
random.effects(object, \dots)
ranef(object, \dots)
}
\alias{random.effects}
\alias{ranef}
\alias{print.ranef}
\arguments{
 \item{object}{any fitted model object from which random effects
   estimates can be extracted.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{lmList} and \code{lme}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\seealso{\code{\link{ranef.lmList}},\code{\link{ranef.lme}}} 
\examples{
## see the method function documentation
}
\keyword{models}
