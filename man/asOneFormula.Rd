% $Id: asOneFormula.Rd,v 1.5 2002/03/05 14:59:39 bates Exp $
\name{asOneFormula}
\title{Combine Formulas of a Set of Objects}
\usage{
asOneFormula(..., omit)
}
\alias{asOneFormula}
\arguments{
 \item{\dots}{objects, or lists of objects, from which a formula can be
   extracted.}
 \item{omit}{an optional character vector with the names of variables to
   be omitted from the returned formula. Defaults to c(".", "pi").}
}
\description{
  The names of all variables used in the formulas extracted from the
  objects defined in \code{...} are converted into a single linear
  formula, with the variables names separated by \code{+}. 
}
\value{
  a one-sided linear formula with all variables named in the formulas
  extracted from the objects in \code{...}, except the ones listed in
  \code{omit}. 
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{formula}}, \code{\link{all.vars}}}

\examples{
asOneFormula(y ~ x + z | g, list(~ w, ~ t * sin(2 * pi)))
}
\keyword{models}
