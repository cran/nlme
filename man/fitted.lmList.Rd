% $Id: fitted.lmList.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{fitted.lmList}
\title{Extract lmList Fitted Values}
\usage{
\method{fitted}{lmList}(object, subset, asList)
}
\alias{fitted.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{subset}{an optional character or integer vector naming the
   \code{lm} components of \code{object} from which the fitted values
   are to be extracted. Default is \code{NULL}, in which case all
   components are used.
 }
 \item{asList}{an optional logical value. If \code{TRUE}, the returned
   object is a list with the fitted values split by groups; else the
   returned value is a vector. Defaults to \code{FALSE}.
 }
}
\description{
  The fitted values are extracted from each \code{lm} component of
  \code{object} and arranged into a list with as many components as
  \code{object}, or combined into a single vector.
}
\value{
  a list with components given by the fitted values of each \code{lm}
  component of \code{object}, or a vector with the fitted values for all
  \code{lm} components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}, \code{\link{residuals.lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
fitted(fm1)
}
\keyword{models}
