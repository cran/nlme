% $Id: pruneLevels.Rd,v 1.5 2002/03/05 14:59:39 bates Exp $
\name{pruneLevels}
\title{Prune Factor Levels}
\usage{
pruneLevels(object)
}
\alias{pruneLevels}
\alias{pruneLevels.factor}
\alias{pruneLevels.ordered}
\arguments{
 \item{object}{an object inheriting from class \code{factor}.}
}
\description{
  The \code{levels} attribute of \code{object} are pruned to contain
  only the levels occurring in the factor.
}
\value{
  an object identical to \code{object}, but with the \code{levels}
  attribute containing only value occurring in the factor.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{factor}}, \code{\link{ordered}}}

\examples{
f1 <- factor(c(1,1,2,3,3,4,5))
levels(f1)
f2 <- f1[4:7]
levels(f2)
levels(pruneLevels(f2))
}
\keyword{models}
