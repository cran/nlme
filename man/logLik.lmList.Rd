% $Id: logLik.lmList.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{logLik.lmList}
\title{Log-Likelihood of an lmList Object}
\usage{
logLik(object, REML, pool)
}
\alias{logLik.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{REML}{an optional logical value. If \code{TRUE} the restricted
   log-likelihood is returned, else, if \code{FALSE}, the log-likelihood
   is returned. Defaults to \code{FALSE}.
 }
 \item{pool}{an optional logical value indicating whether all \code{lm}
   components of \code{object} may be assumed to have the same error
   variance. Default is \code{attr(object, "pool")}.
 }
}
\description{
 If \code{pool=FALSE}, the (restricted) log-likelihoods of the \code{lm}
 components of \code{object} are summed together. Else, the (restricted)
 log-likelihood of the \code{lm} fit with different coefficients for
 each level of the grouping factor associated with the partitioning of
 the \code{object} components is obtained.
}
\value{
  either the sum of the (restricted) log-likelihoods of each \code{lm}
  component in \code{object}, or the (restricted) log-likelihood for the
  \code{lm} fit with separate coefficients for each component of
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
library(nlme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
logLik(fm1)   # returns NA when it should not
}
\keyword{models}
