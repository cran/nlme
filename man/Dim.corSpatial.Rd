% $Id: Dim.corSpatial.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{Dim.corSpatial}
\title{Dimensions of a corSpatial Object}
\usage{
\method{Dim}{corSpatial}(object, groups)
}
\alias{Dim.corSpatial}
\arguments{
 \item{object}{an object inheriting from class \code{corSpatial},
   representing a spatial correlation structure.}
 \item{groups}{an optional factor defining the grouping of the
   observations; observations within a group are correlated and
   observations in different groups are uncorrelated.}
}
\description{
  if \code{groups} is missing, it returns the \code{Dim} attribute of
  \code{object}; otherwise, calculates the dimensions associated with
  the grouping factor.
}
\value{
  a list with components:
  \item{N}{length of \code{groups}}
  \item{M}{number of groups}
  \item{spClass}{an integer representing the spatial correlation class;
    0 = user defined class, 1 = \code{corSpher}, 2 = \code{corExp}, 3 =
    \code{corGaus}, 4 = \code{corLin}}
  \item{sumLenSq}{sum of the squares of the number of observations per
    group} 
  \item{len}{an integer vector with the number of observations per
    group}
  \item{start}{an integer vector with the starting position for the
    distance vectors in each group, beginning from zero}
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{Dim}}, \code{\link{Dim.corStruct}}}

\examples{
data(Orthodont)
Dim(corGaus(), getGroups(Orthodont))
}
\keyword{models}
