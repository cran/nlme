% $Id: corClasses.Rd,v 1.3 2000/07/03 18:22:46 bates Exp $
\name{corClasses}
\title{Correlation Structure Classes}
\alias{corClasses}
\alias{corIdent}
\description{
  Standard classes of correlation structures (\code{corStruct})
  available in the \code{nlme} library. 
}
\value{
  Available standard classes:
  \item{corAR1}{autoregressive process of order 1.}
  \item{corARMA}{autoregressive moving average process, with arbitrary
    orders for the autoregressive and moving average components.}
  \item{corCAR1}{continuous autoregressive process (AR(1) process for a
    continuous time covariate).}
  \item{corCompSymm}{compound symmetry structure corresponding to a
    constant correlation.}
  \item{corExp}{exponential spatial correlation.}
  \item{corGaus}{Gaussian spatial correlation.}
  \item{corLin}{linear spatial correlation.}
  \item{corRatio}{Rational quadratics spatial correlation.}
  \item{corSpher}{spherical spatial correlation.}
  \item{corSymm}{general correlation matrix, with no additional
    structure.} 
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\note{
  Users may define their own \code{corStruct} classes by specifying a
  \code{constructor} function and, at a minimum, methods for the
  functions \code{corMatrix} and \code{coef}. For
  examples of these functions, see the methods for classes \code{corSymm}
  and \code{corAR1}. 
}
\seealso{\code{\link{corAR1}}, \code{\link{corARMA}},
  \code{\link{corCAR1}}, \code{\link{corCompSymm}},
  \code{\link{corExp}}, \code{\link{corGaus}}, \code{\link{corLin}},
  \code{\link{corRatio}}, \code{\link{corSpher}}, \code{\link{corSymm}}}
\keyword{models}
