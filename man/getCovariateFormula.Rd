\name{getCovariateFormula}
\title{Extract Covariates Formula}
\usage{
getCovariateFormula(object)
}
\alias{getCovariateFormula}
\arguments{
 \item{object}{any object from which a formula can be extracted.}
}
\description{
  The right hand side of \code{formula(object)}, without any
  conditioning expressions (i.e. any expressions after a \code{|}
  operator) is returned as a one-sided formula.
}
\value{
  a one-sided formula describing the covariates associated with
  \code{formula(object)}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getCovariate}}}

\examples{
getCovariateFormula(y ~ x | g)
getCovariateFormula(y ~ x)
}
\keyword{models}

