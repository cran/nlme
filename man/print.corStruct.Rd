% $Id: print.corStruct.Rd,v 1.7 2001/10/30 20:51:15 bates Exp $
\name{print.corStruct}
\title{Print a corStruct Object}
\usage{
\method{print}{corStruct}(x, \dots)
\method{print}{corSymm}(x, \dots)
%\method{print}{corHF}(x, \dots)
}
\alias{print.corStruct}
\alias{print.corSymm}
%\alias{print.corHF}
\arguments{
 \item{x}{an object inheriting from class \code{corStruct}, representing
   a correlation structure.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  If \code{x} has been initialized, its coefficients are printed in
  constrained form. 
}
\value{
  the printed coefficients of \code{x} in constrained form.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{print.default}}, \code{\link{coef.corStruct}}}
\examples{
cs1 <- corAR1(0.3)
print(cs1)
}
\keyword{models}
