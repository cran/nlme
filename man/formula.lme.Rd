% $Id: formula.lme.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{formula.lme}
\title{Extract lme Object Formula}
\usage{
formula(object)
}
\alias{formula.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
}
\description{
  This method function extracts the fixed effects model formula
  associated with \code{object}. 
}
\value{
  a two-sided linear formula specifying the fixed effects model used to
  obtain \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}}

\examples{
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
formula(fm1)
}
\keyword{models}
