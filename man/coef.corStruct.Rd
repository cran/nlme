% $Id: coef.corStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{coef.corStruct}
\title{Coefficients of a corStruct Object}
\usage{
coef(object, unconstrained)
coef(object) <- value
}
\alias{coef.corStruct}
\alias{coef.corAR1}
\alias{coef.corARMAd}   
\alias{coef.corCAR1}   
\alias{coef.corCompSymm}
\alias{coef.corHF}
\alias{coef.corIdent}  
\alias{coef.corLin}
\alias{coef.corSpatial}
\alias{coef.corSpher}
\alias{coef.corSymm}
\alias{coef<-.corStruct}
\alias{coef<-.corAR1}
\alias{coef<-.corARMAd}   
\alias{coef<-.corCAR1}   
\alias{coef<-.corCompSymm}
\alias{coef<-.corHF}
\alias{coef<-.corIdent}  
\alias{coef<-.corLin}
\alias{coef<-.corSpatial}
\alias{coef<-.corSpher}
\alias{coef<-.corSymm}

\arguments{
 \item{object}{an object inheriting from class \code{corStruct},
   representing a correlation structure.}
 \item{unconstrained}{a logical value. If \code{TRUE} the coefficients
   are returned in unconstrained form (the same used in the optimization
   algorithm). If \code{FALSE} the coefficients are returned in
   "natural", possibly constrained, form. Defaults to \code{TRUE}.}
 \item{value}{a vector with the replacement values for the coefficients
   associated with \code{object}. It must be a vector with the same length
   of \code{coef{object}} and must be given in unconstrained form.}
}
\description{
  This method function extracts the coefficients associated with the
  correlation structure represented by \code{object}. 
}
\value{
  a vector with the coefficients corresponding to \code{object}.
}
\author{Jose Pinheiro and Douglas Bates }
\section{SIDE EFFECTS}{
  On the left side of an assignment, sets the values of the coefficients
  of \code{object} to \code{value}. \code{Object} must be initialized (using
  \code{initialize}) before new values can be assigned to its
  coefficients. 
}

\seealso{\code{\link{corAR1}}, \code{\link{corARMA}},
    \code{\link{corCAR1}}, \code{\link{corCompSymm}},
    \code{\link{corExp}}, \code{\link{corGaus}}, \code{\link{corLin}},
    \code{\link{corRatio}}, \code{\link{corSpatial}}, \code{\link{corSpher}},
    \code{\link{corSymm}},\code{\link{initialize}}}
\examples{
library(nlme)
cst1 <- corARMA(p = 1, q = 1)
coef(cst1)
}
\keyword{models}
