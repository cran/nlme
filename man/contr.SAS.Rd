\name{contr.SAS}
\alias{contr.SAS}
\title{Contrast Matrices}
\description{
  Return a matrix of contrasts corresponding to the parameterization for
  a factor used in SAS.
}
\usage{
contr.SAS(n, contrasts=TRUE)
}
\arguments{
  \item{n}{a vector of levels for a factor, or the number of levels.}
  \item{contrasts}{a logical indicating whether contrasts should be
    computed.}
}
\details{
  This function is used for creating contrast matrices for use in
  fitting analysis of variance and regression models.  The columns of
  the resulting matrix contain contrasts which can be used for coding
  a factor with \code{n} levels.  The returned value contains the
  computed contrasts.  If the argument \code{contrasts} is \code{FALSE}
  then a square indicator matrix is returned.

  These contrasts are similar to those from \code{contr.treatment}
  except that they represent shifts relative to the last level of the
  factor, not the first level.  This generates the same parameterization for a
  linear model as is used in SAS (see \url{http://www.sas.com}).
}
\value{
  A matrix with \code{n} rows and \code{k} columns, with \code{k=n-1} if
  \code{contrasts} is \code{TRUE} and \code{k=n} if \code{contrasts} is
  \code{FALSE}.
}
\author{Jose Pinheiro \email{jcp@research.bell-labs.com}, Douglas Bates
  \email{bates@stat.wisc.edu}} 
\seealso{\code{\link[base:contrast]{contr.treatment}}}
\examples{
contr.SAS(4)
}
\keyword{design}
\keyword{regression}
\keyword{array}
