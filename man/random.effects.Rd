% File nlme/man/random.effects.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{random.effects}
\title{Extract Random Effects}
\usage{
random.effects(object, \dots)
ranef(object, \dots)
}
\alias{random.effects}
\alias{ranef}
%% \alias{print.ranef}
\arguments{
 \item{object}{any fitted model object from which random effects
   estimates can be extracted.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{lmList} and \code{lme}.

  \code{random.effects} is an alias for \code{ranef}; methods are
  implemented for the latter.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\references{
 Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models
 in S and S-PLUS", Springer, esp. pp. 100, 461.  
}
\seealso{\code{\link{ranef.lmList}}, \code{\link{ranef.lme}}}
\keyword{models}
