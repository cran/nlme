% $Id: formula.groupedData.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{formula.groupedData}
\title{Extract groupedData Formula}
\usage{
formula(object)
}
\alias{formula.groupedData}
\arguments{
 \item{object}{an object inheriting from class \code{groupedData}.}
}
\description{
  This method function extracts the display formula associated with a
  \code{groupedData} object. This is a two-sided formula of the form
  \code{resp ~ cov | group}, where \code{resp} is the response,
  \code{cov} is the primary covariate, and \code{group} is the grouping
  structure. 
}
\value{
  a two-sided formula with a conditioning expression, representing the
  display formula for \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{groupedData}}}

\examples{
data(Orthodont)
formula(Orthodont)
}
\keyword{models}
