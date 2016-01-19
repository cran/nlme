% File nlme/man/VarCorr.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{VarCorr}
\title{Extract variance and correlation components}
\alias{VarCorr}
\alias{VarCorr.lme}
\alias{VarCorr.pdMat}
\alias{VarCorr.pdBlocked}
\alias{print.VarCorr.lme}
\usage{
VarCorr(x, sigma = 1, \dots)
\method{VarCorr}{lme}(x, sigma = 1, rdig = 3, \dots)
## and identical signature for classes 'pdMat' and 'pdBlocked'
}
\arguments{
  \item{x}{a fitted model object, usually an object inheriting from
    class \code{"\link{lme}"}.}
  \item{sigma}{an optional numeric value used as a multiplier for the
    standard deviations.  Default is \code{1}.}
  \item{rdig}{an optional integer value specifying the number of digits
    used to represent correlation estimates.  Default is \code{3}.}
  \item{\dots}{further optional arguments passed to other methods (none
    for the methods documented here).}
}
\description{
  This function calculates the estimated variances, standard
  deviations, and correlations between the random-effects terms in a
  linear mixed-effects model, of class \code{"\link{lme}"}, or a nonlinear
  mixed-effects model, of class \code{"\link{nlme}"}.  The within-group error
  variance and standard deviation are also calculated.
}
\value{
    a matrix with the estimated variances, standard deviations, and
    correlations for the random effects.  The first two columns, named
    \code{Variance} and \code{StdDev}, give, respectively, the variance
    and the standard deviations.  If there are correlation components in
    the random effects model, the third column, named \code{Corr},
    and the remaining unnamed columns give the estimated correlations
    among random effects within the same level of grouping.  The
    within-group error variance and standard deviation are included as
    the last row in the matrix.
}
\references{
  Pinheiro, J.C., and Bates, D.M. (2000)
  \emph{Mixed-Effects Models in S and S-PLUS}, Springer, esp. pp. 100, 461.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{lme}}, \code{\link{nlme}}}
\examples{
fm1 <- lme(distance ~ age, data = Orthodont, random = ~age)
VarCorr(fm1)
}
\keyword{models}
