% $Id: varComb.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
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

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{varWeights.varComb}}, \code{\link{coef.varComb}}}

\examples{
library(nlme)
vf1 <- varComb(varIdent(form = ~1|Sex), varPower())
}
\keyword{models}
