% $Id: Names.formula.Rd,v 1.1 1999/10/13 00:50:10 saikat Exp $
\name{Names.formula}
\title{Extract Names from a formula}
\usage{
Names(object, data, exclude)
}
\alias{Names.formula}
\arguments{
 \item{object}{an object inheriting from class \code{formula}.}
 \item{data}{an optional data frame containing the variables specified
   in \code{object}. By default the variables are taken from the
   environment from which \code{Names.formula} is called.}
 \item{exclude}{an optional character vector with names to be excluded
   from the returned value. Default is \code{c("pi",".")}.}
}
\description{
  This method function returns the names of the terms corresponding to
  the right hand side of \code{object} (treated as a linear formula),
  obtained as the column names of the corresponding
  \code{model.matrix}. 
}
\value{
  a character vector with the column names of the \code{model.matrix}
  corresponding to the right hand side of \code{object} which are not
  listed in \code{excluded}.  
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{model.matrix}}, \code{\link{terms}},
  \code{\link{Names}}
}

\examples{
data(Orthodont)
Names(distance ~ Sex * age, data = Orthodont)
}
\keyword{models}
