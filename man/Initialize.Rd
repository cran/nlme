% File nlme/man/Initialize.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{Initialize}
\title{Initialize Object}
\usage{
Initialize(object, data, \dots)
}
\alias{Initialize}
\arguments{
 \item{object}{any object requiring initialization, e.g. "plug-in"
   structures such as \code{corStruct} and \code{varFunc} objects. }
 \item{data}{a data frame to be used in the initialization procedure.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include: \code{corStruct}, \code{glsStruct}, \code{lmeStruct},
  \code{reStruct}, and \code{varFunc}.
}
\value{
  an initialized object with the same class as \code{object}. Changes
  introduced by the initialization procedure will depend on the method
  function used; see the appropriate documentation. 
}
\references{
  Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models
  in S and S-PLUS", Springer.  
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{Initialize.corStruct}},
  \code{\link{Initialize.glsStruct}},
  \code{\link{Initialize.lmeStruct}},
  \code{\link{Initialize.reStruct}},
  \code{\link{Initialize.varFunc}},
  \code{\link{isInitialized}}
}
\keyword{models}
