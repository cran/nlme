% $Id: formula.lme.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{formula.lme}
\title{Extract lme Object Formula}
\usage{
\method{formula}{lme}(x, \dots)
}
\alias{formula.lme}
\arguments{
 \item{x}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the fixed effects model formula
  associated with \code{x}. 
}
\value{
  a two-sided linear formula specifying the fixed effects model used to
  obtain \code{x}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{lme}}}

\examples{
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
formula(fm1)
}
\keyword{models}
