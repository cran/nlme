% $Id: formula.groupedData.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{formula.groupedData}
\title{Extract groupedData Formula}
\usage{
\method{formula}{groupedData}(x, \dots)
}
\alias{formula.groupedData}
\arguments{
 \item{x}{an object inheriting from class \code{groupedData}.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
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
  display formula for \code{x}.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{groupedData}}}

\examples{
data(Orthodont)
formula(Orthodont)
}
\keyword{models}
