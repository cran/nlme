% $Id: model.matrix.reStruct.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{model.matrix.reStruct}
\title{reStruct Model Matrix}
\usage{
model.matrix(object, data, contr)
}
\alias{model.matrix.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{formula(object)}.}
 \item{contr}{an optional named list specifying the contrasts to be used
   for representing the \code{factor} variables in \code{data}. The
   components names should match the names of the variables in
   \code{data} for which the contrasts are to be specified. The
   components of this list will be used as the \code{contrasts}
   attribute of the corresponding factor. If missing, the default
   contrast specification is used.}  
}
\description{
  The model matrices for each element of \code{formula(object)},
  calculated using \code{data}, are bound together column-wise. When
  multiple grouping levels are present (i.e. when \code{length(object) >
    1}), the individual model matrices are combined from innermost (at
  the leftmost position) to outermost (at the rightmost position). 
}
\value{
  a matrix obtained by binding together, column-wise, the model matrices
  for each element of \code{formula(object)}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{model.matrix}}, \code{\link{contrasts}},
  \code{\link{reStruct}}, \code{\link{formula.reStruct}}} 

\examples{
data(Pixel)
rs1 <- reStruct(list(Dog = ~day, Side = ~1), data = Pixel)
model.matrix(rs1, Pixel)
}
\keyword{models}
