% $Id: getGroups.Rd,v 1.1 2000/03/17 22:21:22 saikat Exp $
\name{getGroups}
\title{Extract Grouping Factors from an Object}
\usage{
getGroups(object, form, level, data)
}
\alias{getGroups}
\arguments{
 \item{object}{any object}
 \item{form}{an optional formula with a conditioning expression on its
   right hand side (i.e. an expression involving the \code{|}
   operator). Defaults to \code{formula(object)}.}
 \item{level}{a positive integer vector with the level(s) of grouping to
   be used when multiple nested levels of grouping are present. This
   argument is optional for most methods of this generic function and
   defaults to all levels of nesting.}
 \item{data}{a data frame in which to interpret the variables named in
   \code{form}. Optional for most methods.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{corStruct}, \code{data.frame},
  \code{gls}, \code{lme}, \code{lmList}, and \code{varFunc}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getGroupsFormula}}, \code{\link{getGroups.data.frame}},
\code{\link{getGroups.gls}}, \code{\link{getGroups.lmList}},
\code{\link{getGroups.lme}}}
\examples{
## see the method function documentation
}
\keyword{models}
