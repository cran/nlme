% $Id: getCovariate.data.frame.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{getCovariate.data.frame}
\title{Extract Data Frame Covariate}
\usage{
getCovariate(object, form)
}
\alias{getCovariate.data.frame}
\arguments{
 \item{object}{an object inheriting from class \code{data.frame}.}
 \item{form}{an optional formula specifying the covariate to be
   evaluated in \code{object}. Defaults to \code{formula(object)}.}
}
\description{
  The right hand side of \code{form}, stripped of any conditioning
  expression (i.e. an expression following a \code{|} operator), is
  evaluated in \code{object}.  
}
\value{
  the value of the right hand side of \code{form}, stripped of
  any conditional expression, evaluated in \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getCovariateFormula}}}

\examples{
data(Orthodont)
getCovariate(Orthodont)
}
\keyword{models}
