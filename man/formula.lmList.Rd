% $Id: formula.lmList.Rd,v 1.4 2000/07/03 18:22:46 bates Exp $
\name{formula.lmList}
\title{Extract lmList Object Formula}
\usage{
formula(object)
}
\alias{formula.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  This method function extracts the common linear model formula
  associated with each \code{lm} component of \code{object}. 
}
\value{
  a two-sided linear formula specifying the linear model used to
  obtain the \code{lm} components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
formula(fm1)
}
\keyword{models}
