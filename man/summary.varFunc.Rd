% $Id: summary.varFunc.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{summary.varFunc}
\title{Summarize varFunc Object}
\usage{
\method{summary}{varFunc}(object, structName)
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

\author{Jose Pinheiro and Douglas Bates}


\seealso{\code{\link{print.summary.varFunc}}}

\examples{
data(Orthodont)
vf1 <- varPower(0.3, form = ~age)
vf1 <- initialize(vf1, Orthodont)
summary(vf1)
}
\keyword{models}
