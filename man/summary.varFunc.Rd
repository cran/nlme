% $Id: summary.varFunc.Rd,v 1.7 2002/03/05 14:59:39 bates Exp $
\name{summary.varFunc}
\title{Summarize varFunc Object}
\usage{
\method{summary}{varFunc}(object, structName, \dots)
}
\alias{summary.varFunc}
\alias{summary.varComb}
\alias{summary.varConstPower}
\alias{summary.varExp}
\alias{summary.varFixed}
\alias{summary.varIdent}
\alias{summary.varPower}

\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
 \item{structName}{an optional character string with a description of
   the \code{varFunc} class. Default depends on the method function:
 \code{"Combination of variance functions"} for \code{varComb},
 \code{"Constant plus power of covariate"} for \code{varConstPower},
 \code{"Exponential of variance covariate"} for \code{varExp},
 \code{"Different standard deviations per stratum"} for \code{varIdent},
 \code{"Power of variance covariate"} for \code{varPower}, and
 \code{data.class(object)} for \code{varFunc}.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  A \code{structName} attribute, with the value of corresponding
  argument, is appended to \code{object} and its class is changed to
  \code{summary.varFunc}.  
}
\value{
  an object similar to \code{object}, with an additional attribute
  \code{structName}, inheriting from class \code{summary.varFunc}. 
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}


\seealso{\code{\link{print.summary.varFunc}}}

\examples{
data(Orthodont)
vf1 <- varPower(0.3, form = ~age)
vf1 <- initialize(vf1, Orthodont)
summary(vf1)
}
\keyword{models}
