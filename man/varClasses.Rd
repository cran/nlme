% $Id: varClasses.Rd,v 1.4 2000/07/03 18:22:49 bates Exp $
\name{varClasses}
\title{Variance Function Classes}
\alias{varClasses}

\description{
  Standard classes of variance function structures (\code{varFunc}) 
  available in the \code{nlme} library. Covariates included in the
  variance function, denoted by variance covariates, may involve
  functions of the fitted model object, such as the fitted values and
  the residuals. Different coefficients may be assigned to the levels of
  a classification factor.
}
\value{
  Available standard classes ():
  \item{varExp}{exponential of a variance covariate.}
  \item{varPower}{power of a variance covariate.}
  \item{varConstPower}{constant plus power of a variance covariate.}
  \item{varIdent}{constant variance(s), generally used to allow
    different variances according to the levels of a classification
    factor.} 
  \item{varFixed}{fixed weights, determined by a variance covariate.}
  \item{varComb}{combination of variance functions.}
}
\author{Jose Pinheiro and Douglas Bates}
\note{
  Users may define their own \code{varFunc} classes by specifying a
  \code{constructor} function and, at a minimum, methods for the
  functions \code{coef}, \code{coef<-}, and \code{initialize}. For
  examples of these functions, see the methods for class
  \code{varPower}. 
}
\seealso{\code{\link{varExp}}, \code{\link{varPower}},
  \code{\link{varConstPower}}, \code{\link{varIdent}},
  \code{\link{varFixed}}, \code{\link{varComb}}
}
\keyword{models}
