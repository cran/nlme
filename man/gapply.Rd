% $Id: gapply.Rd,v 1.1 2000/03/17 22:21:22 saikat Exp $
\name{gapply}
\title{Apply a Function by Groups}
\usage{
gapply(object, which, FUN, form, level, groups, ...)
}
\alias{gapply}
\arguments{
 \item{object}{an object to which the function will be applied - usually
   a \code{groupedData} object or a \code{data.frame}. Must inherit from
   class \code{data.frame}.
 }
 \item{which}{an optional character or positive integer vector
   specifying which columns of \code{object} should be used with
   \code{FUN}. Defaults to all columns in \code{object}.
 }
 \item{FUN}{function to apply to the distinct sets of rows
   of the data frame \code{object} defined by the values of
   \code{groups}.
 } 
 \item{form}{an optional one-sided formula that defines the groups.
   When this formula is given the right-hand side is evaluated in
   \code{object}, converted to a factor if necessary, and the unique
   levels are used to define the groups.  Defaults to
   \code{formula(object)}.
 }
 \item{level}{an optional positive integer giving the level of grouping
   to be used in an object with multiple nested grouping levels.
   Defaults to the highest or innermost level of grouping.
 }
 \item{groups}{an optional factor that will be used to split the 
   rows into groups.  Defaults to \code{getGroups(object, form, level)}.
 }
 \item{\dots}{optional additional arguments to the summary function
   \code{FUN}.  Often it is helpful to specify \code{na.rm = TRUE}.
 }
}

\description{
  Applies the function to the distinct sets of rows of the data frame
  defined by \code{groups}.
}

\value{
  Returns a data frame with as many rows as there are levels in the
  \code{groups} argument.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gsummary}}}

\examples{
\dontrun{
library( nlme )
data( Quinidine )
## Find number of non-missing "conc" observations for each Subject
gapply( Quinidine, FUN = function(x) sum(!is.na(x$conc)) )
}
}
\keyword{data}
