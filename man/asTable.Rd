% $Id: asTable.Rd,v 1.5 2001/01/10 19:21:04 bates Exp $
\name{asTable}
\alias{asTable}
\alias{asTable.groupedData}
\title{Convert groupedData to a matrix}
\usage{
asTable(object)
}
\arguments{
 \item{object}{A balanced \code{groupedData} object}
}
\description{
  Create a tabular representation of the response in a balanced
  groupedData object.
}
\details{
  A balanced groupedData object can be represented as a matrix or table
  of response values corresponding to the values of a primary covariate
  for each level of a grouping factor.  This function creates such a
  matrix representation of the data in \code{object}.
}
\value{
  A matrix.  The data in the matrix are the values of the response.  The
  columns correspond to the distinct values of the primary covariate and
  are labelled as such.  The rows correspond to the distinct levels of
  the grouping factor and are labelled as such.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{groupedData}}, \code{\link{isBalanced}},
  \code{\link{balancedGrouped}}}
\examples{
data(Orthodont)
asTable(Orthodont)
}
\keyword{ data }
