% $Id: formula.lmList.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{formula.lmList}
\title{Extract lmList Object Formula}
\usage{
\method{formula}{lmList}(x, \dots)
}
\alias{formula.lmList}
\arguments{
 \item{x}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the common linear model formula
  associated with each \code{lm} component of \code{x}. 
}
\value{
  a two-sided linear formula specifying the linear model used to
  obtain the \code{lm} components of \code{x}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
formula(fm1)
}
\keyword{models}
