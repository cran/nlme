% $Id: formula.lme.Rd,v 1.4 2000/07/03 18:22:46 bates Exp $
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
