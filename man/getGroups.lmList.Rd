% $Id: getGroups.lmList.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{getGroups.lmList}
\title{Extract lmList Object Groups}
\usage{
\method{getGroups}{lmList}(object, form, level, data, sep)
}
\alias{getGroups.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.}
 \item{form}{an optional formula with a conditioning expression on its
   right hand side (i.e. an expression involving the \code{|}
   operator). Defaults to \code{formula(object)}.  Not used.}
 \item{level}{a positive integer vector with the level(s) of grouping to
   be used when multiple nested levels of grouping are present. This
   argument is optional for most methods of this generic function and
   defaults to all levels of nesting.  Not used.}
 \item{data}{a data frame in which to interpret the variables named in
   \code{form}. Optional for most methods.  Not used.}
 \item{sep}{character, the separator to use between group levels when
   multiple levels are collapsed.  The default is \code{'/'}.  Not used.}
}
\description{
  The grouping factor determining the partitioning of the observations
  used to produce the \code{lm} components of \code{object} is
  extracted.  
}
\value{
  a vector with the grouping factor corresponding to the \code{lm}
  components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
getGroups(fm1)
}
\keyword{models}
