% $Id: VarCorr.Rd,v 1.8 2002/03/05 14:59:39 bates Exp $
\name{VarCorr}
\title{Extract variance and correlation components}
\usage{
VarCorr(x, sigma, rdig)
}
\alias{VarCorr}
\alias{VarCorr.lme}
\alias{VarCorr.pdMat}
\alias{VarCorr.pdBlocked}
\alias{summary.VarCorr.lme}
\alias{print.VarCorr.lme}
\arguments{
 \item{x}{a fitted model object, usually an object inheriting from
     class \code{lme}. 
 }
 \item{sigma}{an optional numeric value used as a multiplier for the
     standard deviations. Default is \code{1}.
 }
 \item{rdig}{an optional integer value specifying the number of digits
     used to represent correlation estimates. Default is \code{3}.
 }
}
\description{
  This function calculates the estimated variances, standard
  deviations, and correlations between the random-effects terms in a
  linear mixed-effects model, of class \code{lme}, or a nonlinear
  mixed-effects model, of class \code{nlme}. The within-group error
  variance and standard deviation are also calculated. 
}
\value{
    a matrix with the estimated variances, standard deviations, and
    correlations for the random effects. The first two columns, named
    \code{Variance} and \code{StdDev}, give, respectively, the variance
    and the standard deviations. If there are correlation components in
    the random effects model, the third column, named \code{Corr},
    and the remaining unnamed columns give the estimated correlations
    among random effects within the same level of grouping. The
    within-group error variance and standard deviation are included as
    the last row in the matrix. 
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{lme}}, \code{\link{nlme}}}
\examples{
data(Orthodont)
fm1 <- lme(distance ~ age, data = Orthodont, random = ~age)
VarCorr(fm1)
}
\keyword{models}
