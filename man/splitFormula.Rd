% $Id: splitFormula.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
\name{splitFormula}
\title{Split a Formula}
\usage{
splitFormula(frm, nm)
}
\alias{splitFormula}
\arguments{
 \item{form}{a \code{formula} object.}
 \item{sep}{an optional character string specifying the separator to be
   used for splitting the formula. Defaults to \code{"/"}. }
}
\description{
  Splits the right hand side of \code{form} into a list of subformulas
  according to the presence of \code{sep}. The left hand side of
  \code{form}, if present, will be ignored. The length of the returned
  list will be equal to the number of occurrences of \code{sep} in
  \code{form} plus one. 
}
\value{
  a list of formulas, corresponding to the split of \code{form}
  according to \code{sep}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{formula}}}

\examples{
splitFormula(~ g1/g2/g3)
}
\keyword{models}
