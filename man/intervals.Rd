% $Id: intervals.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{intervals}
\title{Confidence Intervals on Coefficients}
\usage{
intervals(object, level, ...)
}
\alias{intervals}
\arguments{
 \item{object}{a fitted model object from which parameter estimates can
   be extracted.}
 \item{level}{an optional numeric value for the interval confidence
   level. Defaults to 0.95.}
 \item{\dots}{some methods for the generic may require additional
   arguments.}
}
\description{
  Confidence intervals on the parameters associated with the model
  represented by \code{object} are obtained. This function is generic;
  method functions can be written to handle specific classes of
  objects. Classes which already have methods for this function include:
  \code{gls}, \code{lme}, and \code{lmList}. 
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\seealso{\code{\link{intervals.gls}}, \code{\link{intervals.lme}},
    \code{\link{intervals.lmList}}
} 

\author{Jose Pinheiro and Douglas Bates}
\examples{
## see the method documentation
}
\keyword{models}
