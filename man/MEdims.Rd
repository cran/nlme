\name{MEdims}
\alias{MEdims}
\title{Calculate dimensions in mixed-effects}
\description{
  Calculate the dimensions of various arrays used to represent a
  mixed-effects model.}
}
\usage{
MEdims(groups, ncols)
}
\arguments{
  \item{groups}{a factor or list of factors describing the grouping of
    the observations}
  \item{ncols}{a vector of positive integers listing the number of
    columns in the model matrices at each level of random effects terms
    or fixed effects
}
\value{
  a list of dimensions
  \item{N}{total number of observations in the original data.}
  \item{ZXrows}{number of rows in the \code{ZXy} array.  This can be
    different from \code{N} if the condensed linear model has been
    decomposed by \code{MEdecomp}.}
  \item{ZXcols}{number of columns in the \code{ZXy} array.}
  \item{Q}{number of levels of random effects.}
  \item{StrRows}{total number of rows in the arrays to be stored if
    conditional estimates of random effects and fixed effects are to be
    calculated.}
  \item{qvec}{a vector of length \code{Q+2} giving the number of columns
    in the model matrices for each level of random effects, the fixed
    effects and the response.}
  \item{ngrps}{a vector of length \code{Q} giving the number of groups
    at each level of the random effects.}
  \item{DmOff}{a vector of length \code{Q} giving the offsets into the
    concatenated relative precision factor arrays.}
  \item{ncols}{a vector of length \code{Q} giving the number of columns
    in the model matrices at each level of random effects.}
  \item{nrot}{a vector of length \code{Q+2} giving the number of columns
    to the right of the current level in the \code{ZXy} array that must
    be rotated in the decomposition.}
  \item{ZXoff}{a list of vectors of length \code{Q+2} giving the offsets
    into the \code{ZXy} array for the beginning of the model matrix in
    each group at each level.}
  \item{ZXlen}{a list of vectors of length \code{Q+2} giving the lengths
    (number of rows) of the model matrix in each group at each level.}
  \item{SToff}{a list of vectors of length \code{Q+2} giving the offsets
    into the stored data array for the beginning of the model matrix in
    each group at each level.}
  \item{DecOff}{a list of vectors of length \code{Q+2} giving the offsets
    into the decomposed data array for the beginning of the model matrix in
    each group at each level.  These offsets differ from \code{SToff}
    when there are nested random effects.  The \code{DecOff} entry for
    an outer grouping factor points to the beginning of the block and the
    \code{SToff} entry points to the beginning of the last section of
    the block.}
  \item{DecLen}{a list of vectors of length \code{Q+2} giving the
    lengths (number of rows) in the block to be decomposed for each
    group at each level.}
}
\author{Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{MEestimate}}}
\keyword{models}

