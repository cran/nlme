% File nlme/man/summary.varFunc.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{summary.varFunc}
\title{Summarize "varFunc" Object}
\usage{
\method{summary}{varFunc}(object, structName, \dots)
}
\alias{summary.varFunc}
\alias{summary.varComb}
\alias{summary.varConstPower}
\alias{summary.varConstProp}
\alias{summary.varExp}
\alias{summary.varFixed}
\alias{summary.varIdent}
\alias{summary.varPower}

\arguments{
 \item{object}{an object inheriting from class \code{"\link{varFunc}"},
   representing a variance function structure.}

 \item{structName}{an optional character string with a description of
   the \code{varFunc} class.  Default depends on the method function:
   \describe{
     \item{for \code{varComb}:}{\code{"Combination of variance functions"},}
     \item{for \code{varConstPower}:}{\code{"Constant plus power of variance covariate"},}
     \item{for \code{varConstProp}:}{\code{"Constant plus proportion of variance covariate"},}
     \item{for \code{varExp}:}{\code{"Exponential of variance covariate"},}
     \item{for \code{varIdent}:}{\code{"Different standard deviations per stratum"},}
     \item{for \code{varPower}:}{\code{"Power of variance covariate"},}
     \item{for \code{varFunc}:}{\code{data.class(object)}.}
   }
 }
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

\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{
  \code{\link{varClasses}},
  \code{\link{varFunc}}
}
\examples{
vf1 <- varPower(0.3, form = ~age)
vf1 <- Initialize(vf1, Orthodont)
summary(vf1)
}
\keyword{models}
