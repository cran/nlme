% $Id: summary.corStruct.Rd,v 1.7 2001/10/30 20:51:15 bates Exp $
\name{summary.corStruct}
%\title{Summarize \code{corStruct} Object}
\title{Summarize a corStruct Object}
\usage{
\method{summary}{corStruct}(object, structName, \dots)
}
\alias{summary.corStruct}
\alias{summary.corAR1}
\alias{summary.corARMA}
\alias{summary.corCAR1}
\alias{summary.corCompSymm}
\alias{summary.corExp}
\alias{summary.corGaus}
\alias{summary.corHF}
\alias{summary.corIdent}
\alias{summary.corLin}
\alias{summary.corNatural}
\alias{summary.corRatio}
\alias{summary.corSpher}
\alias{summary.corSymm}

\arguments{
 \item{object}{an object inheriting from class \code{corStruct},
   representing a correlation structure.}
 \item{structName}{an optional character string defining the type of
   correlation structure associated with \code{object}, to be used in
   the \code{print.summary} method. Defaults to
   \code{class(object)[1]}.} 
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function prepares \code{object} to be printed using the
  \code{print.summary} method, by changing its class and adding a
  \code{structName} attribute to it.
}
\value{
  an object identical to \code{object}, but with its class changed to
  \code{summary.corStruct} and an additional attribute
  \code{structName}. The returned value inherits from the same classes
  as \code{object}.
}
\author{Jose Pinheiro and Douglas Bates }

\seealso{\code{\link{print.summary.corStruct}}}

\examples{
cs1 <- corAR1(0.2)
summary(cs1)
}
\keyword{models}
