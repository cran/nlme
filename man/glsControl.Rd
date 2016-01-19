% File nlme/man/glsControl.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{glsControl}
\title{Control Values for gls Fit}
\usage{
glsControl(maxIter, msMaxIter, tolerance, msTol, msVerbose,
           singular.ok, returnObject, apVar, .relStep,
           opt=c("nlminb", "optim"), optimMethod,
           minAbsParApVar, natural, sigma = NULL)
}
\alias{glsControl}
\arguments{
 \item{maxIter}{maximum number of iterations for the \code{gls}
   optimization algorithm. Default is 50.}
 \item{msMaxIter}{maximum number of iterations
   for the optimization step inside the \code{gls}
   optimization. Default is 50.}
 \item{tolerance}{tolerance for the convergence criterion in the
   \code{gls} algorithm. Default is 1e-6.}
 \item{msTol}{tolerance for the convergence criterion of the first outer
   iteration when \code{optim} is used.  Default is 1e-7.}
 \item{msVerbose}{a logical value passed as the \code{trace} argument to
   \code{ms} (see documentation on that function). Default is
   \code{FALSE}.}
 \item{singular.ok}{a logical value indicating whether non-estimable
   coefficients (resulting from linear dependencies among the columns of
   the regression matrix) should be allowed. Default is \code{FALSE}.}
 \item{returnObject}{a logical value indicating whether the fitted
   object should be returned when the maximum number of iterations is
   reached without convergence of the algorithm. Default is
   \code{FALSE}.}
 \item{apVar}{a logical value indicating whether the approximate
   covariance matrix of the variance-covariance parameters should be
   calculated. Default is \code{TRUE} unless \code{sigma} is fixed (to
   a positive number) in which case the covariance computation is not
   available.}
 \item{.relStep}{relative step for numerical derivatives
   calculations. Default is \code{.Machine$double.eps^(1/3)}.}
 \item{opt}{the optimizer to be used, either \code{"\link{nlminb}"} (the
   current default) or \code{"\link{optim}"} (the previous default).}
 \item{optimMethod}{character - the optimization method to be used with
   the \code{\link{optim}} optimizer. The default is
   \code{"BFGS"}.  An alternative is \code{"L-BFGS-B"}.}
 \item{minAbsParApVar}{numeric value - minimum absolute parameter value
   in the approximate variance calculation.  The default is \code{0.05}.}
 \item{natural}{logical.  Should the natural parameterization be used
   for the approximate variance calculations?  Default is \code{TRUE}.}
 \item{sigma}{optionally a positive number to fix the residual error at.
   If \code{NULL}, as by default, sigma is estimated.}
}
\description{
  The values supplied in the function call replace the defaults and a
  list with all possible arguments is returned. The returned list is
  used as the \code{control} argument to the \code{gls} function.
}
\value{
  a list with components for each of the possible arguments.
}

\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}; the
  \code{sigma} option: Siem Heisterkamp and Bert van Willigen.}

\seealso{\code{\link{gls}}}

\examples{
# decrease the maximum number iterations in the optimization call and
# request that information on the evolution of the ms iterations be printed
glsControl(msMaxIter = 20, msVerbose = TRUE)
}
\keyword{models}
