% File nlme/man/gls-internal.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{gls-internal}
\title{Auxiliary functions used by gls}
\alias{glsApVar}
\alias{glsEstimate}
\description{
  These are functions used by \code{gls} to call its compiled C
  code. They are exported to allow experimentation with modified versions.
}
\usage{
glsApVar(glsSt, sigma, conLin = attr(glsSt, "conLin"),
         .relStep = .Machine$double.eps^(1/3), minAbsPar = 0, natural = TRUE)
glsEstimate(object, conLin = attr(object, "conLin"),
            control = list(singular.ok = FALSE))
}
\arguments{
  \item{glsSt, object}{An object inheriting from  class
    \code{"\link{glsStruct}"}.}
  \item{sigma}{the estimated residual standard error: see
    \code{\link{glsObject}}.}
  \item{conLin}{A \sQuote{condensed linear model}: see
    \code{\link{logLik.glsStruct}}.}
  \item{.relStep, minAbsPar, natural}{Control values: see 
    \code{\link{glsControl}}.}
  \item{control}{The relevant part of a \code{\link{glsControl}} return value.}
}
\keyword{internal}
