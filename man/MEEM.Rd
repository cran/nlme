\name{MEEM}
\alias{MEEM}
\title{EM iterations for mixed-effects}
\description{
  Perform a fixed number of EM iterations on a mixed-effects model
}
\usage{
MEEM(object, conLin, niter=0)
}
\arguments{
  \item{object}{an object that inherits from class \code{reStruct},
    representing the random-effects structure in a mixed-effects model.}
  \item{conLin}{a condensed linear model providing the data and model
    matrices corresponding to the \code{object}}
  \item{niter}{a positive integer - the number of iterations to be preformed.}
}
\value{
  the updated \code{object}
}
\author{Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{MEestimate}}
\keyword{models}
