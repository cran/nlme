% $Id: fitted.gnls.Rd,v 1.7 2002/03/05 14:59:39 bates Exp $
\name{fitted.gnls}
\title{Extract gnls Fitted Values}
\usage{
\method{fitted}{gnls}(object, \dots)
}
\alias{fitted.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlienar least squares fitted model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The fitted values for the nonlinear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the fitted values for the nonlinear model represented by
  \code{object}. 
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gnls}}, \code{\link{residuals.gnls}} }

\examples{
data(Soybean)
%\dontrun{
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
fitted(fm1)
%}
}
\keyword{models}
