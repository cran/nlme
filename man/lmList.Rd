% $Id: lmList.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{lmList}
\title{List of lm Objects with a Common Model}
\usage{
lmList(object, data, level, subset, na.action, pool)
}
\alias{lmList}
\alias{lmList.formula}
\arguments{
  \item{object}{
    either a linear formula object of the form \code{y ~ x1+...+xn | g}
    or a \code{groupedData} object. In the formula object, \code{y}
    represents the response, \code{x1,...,xn} the covariates, and
    \code{g} the grouping factor specifying the partitioning of the data
    according to which different \code{lm} fits should be performed. The
    grouping factor \code{g} may be omitted from the formula, in which
    case the grouping structure will be obtained from \code{data}, which
    must inherit from class \code{groupedData}. The method function 
    \code{lmList.groupedData} is documented separately.
  }
  \item{data}{
    a data frame in which to interpret the variables named in
    \code{object}. 
  }
  \item{level}{
    an optional integer specifying the level of grouping to be used when 
    multiple nested levels of grouping are present.
  }
 \item{subset}{an optional expression indicating which subset of the rows of
   \code{data} should  be  used in the fit. This can be a logical
   vector, or a numeric vector indicating which observation numbers are
   to be included, or a  character  vector of the row names to be
   included.  All observations are included by default.}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{lmList} to print an error message and terminate if there are any
   incomplete observations.
 }
 \item{pool}{
   an optional logical value that is preserved as an attribute of the
   returned value.  This will be used as the default for \code{pool} in
   calculations of standard deviations or standard errors for summaries.
 }
}
\description{
  \code{Data} is partitioned according to the levels of the grouping
  factor \code{g} and individual \code{lm} fits are obtained for each
  \code{data} partition, using the model defined in \code{object}.
}

\value{
  a list of \code{lm} objects with as many components as the number of
  groups defined by the grouping factor. Generic functions such as
  \code{coef}, \code{fixed.effects}, \code{lme}, \code{pairs},
  \code{plot}, \code{predict}, \code{random.effects}, \code{summary},
  and \code{update} have methods that can be applied to an \code{lmList}
  object. 
}

\seealso{
\code{\link{lm}}, \code{\link{lme.lmList}}.
}
\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
summary(fm1)
}
\keyword{models}

