% $Id: predict.lmList.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
\name{predict.lmList}
\title{Predictions from an lmList Object}
\usage{
predict(object, newdata, subset, pool, asList, se.fit)
}
\alias{predict.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{newdata}{an optional data frame to be used for obtaining the
   predictions. All variables used in the \code{object} model formula
   must be present in the data frame. If missing, the same data frame
   used to produce \code{object} is used.
 } 
 \item{subset}{an optional character or integer vector naming the
   \code{lm} components of \code{object} from which the predictions
   are to be extracted. Default is \code{NULL}, in which case all
   components are used.
 }
 \item{asList}{an optional logical value. If \code{TRUE}, the returned
   object is a list with the predictions split by groups; else the
   returned value is a vector. Defaults to \code{FALSE}.
 }
 \item{pool}{an optional logical value indicating whether a pooled
   estimate of the residual standard error should be used. Default is
   \code{attr(object, "pool")}.
 }
 \item{se.fit}{an optional logical value indicating whether pointwise
   standard errors should be computed along with the
   predictions. Default is \code{FALSE}.
 }
}
\description{
  If the grouping factor corresponding to \code{object} is included
  in \code{newdata}, the data frame is partitioned according to the
  grouping factor levels; else, \code{newdata} is repeated for all
  \code{lm} components. The predictions and, optionally, the standard
  errors for the predictions, are obtained for each \code{lm}
  component of \code{object}, using the corresponding element of the
  partitioned \code{newdata}, and arranged into a list with as many
  components as \code{object}, or combined into a single vector or data
  frame (if \code{se.fit=TRUE}).
}
\value{
  a list with components given by the predictions (and, optionally, the
  standard errors for the predictions) from each \code{lm}
  component of \code{object},  a vector with the predictions from all 
  \code{lm} components of \code{object}, or a data frame with columns
  given by the predictions and their corresponding standard errors.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}, \code{\link{predict.lm}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
predict(fm1, se.fit = TRUE)
}
\keyword{models}
