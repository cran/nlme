% $Id: varComb.Rd,v 1.5 2002/03/05 14:59:39 bates Exp $
\name{varComb}
\title{Combination of Variance Functions}
\usage{
varComb(...)
}
\alias{varComb}
\arguments{
 \item{\dots}{objects inheriting from class \code{varFunc} representing
   variance function structures.}
}
\description{
  This function is a constructor for the \code{varComb} class,
  representing a combination of variance functions. The corresponding
  variance function is equal to the product of the variance functions of
  the \code{varFunc} objects listed in \code{\dots}.
}
\value{
  a \code{varComb} object representing a combination of variance
  functions, also inheriting from class \code{varFunc}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{varWeights.varComb}}, \code{\link{coef.varComb}}}

\examples{
vf1 <- varComb(varIdent(form = ~1|Sex), varPower())
}
\keyword{models}
