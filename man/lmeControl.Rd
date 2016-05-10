% File nlme/man/lmeControl.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see ../LICENCE.note

\name{lmeControl}
\title{Specifying Control Values for lme Fit}
\alias{lmeControl}
\description{
  The values supplied in the \code{lmeControl()} call replace the
  defaults, and a \code{\link{list}} with all settings (i.e., values for
  all possible arguments) is returned.  The returned list is
  used as the \code{control} argument to the \code{lme} function.
}
\usage{
lmeControl(maxIter = 50, msMaxIter = 50, tolerance = 1e-6, niterEM = 25,
           msMaxEval = 200,
	   msTol = 1e-7, msVerbose = FALSE,
           returnObject = FALSE, gradHess = TRUE, apVar = TRUE,
	   .relStep = .Machine$double.eps^(1/3), minAbsParApVar = 0.05,
           opt = c("nlminb", "optim"),
	   optimMethod = "BFGS", natural = TRUE,
	   sigma = NULL, \dots)
}
\arguments{
 \item{maxIter}{maximum number of iterations for the \code{lme}
   optimization algorithm.  Default is \code{50}.}
 \item{msMaxIter}{maximum number of iterations
   for the optimization step inside the \code{lme} optimization.
   Default is \code{50}.}
 \item{tolerance}{tolerance for the convergence criterion in the
   \code{lme} algorithm. Default is \code{1e-6}.}
 \item{niterEM}{number of iterations for the EM algorithm used to refine
   the initial estimates of the random effects variance-covariance
   coefficients.  Default is \code{25}.}
 \item{msMaxEval}{maximum number of evaluations of the objective
   function permitted for nlminb.  Default is \code{200}.}
 \item{msTol}{tolerance for the convergence criterion on the first
   iteration when \code{optim} is used.  Default is \code{1e-7}.}
 \item{msVerbose}{a logical value passed as the \code{trace} argument to
   \code{\link{nlminb}} or \code{\link{optim}}.  Default is \code{FALSE}.}
 \item{returnObject}{a logical value indicating whether the fitted
   object should be returned when the maximum number of iterations is
   reached without convergence of the algorithm.  Default is
   \code{FALSE}.}
 \item{gradHess}{a logical value indicating whether numerical gradient
   vectors and Hessian matrices of the log-likelihood function should
   be used in the internal optimization.  This option is only available
   when the correlation structure (\code{corStruct}) and the variance
   function structure (\code{varFunc}) have no "varying" parameters and
   the \code{pdMat} classes used in the random effects structure are
   \code{pdSymm} (general positive-definite), \code{pdDiag} (diagonal),
   \code{pdIdent} (multiple of the identity),  or
   \code{pdCompSymm} (compound symmetry).  Default is \code{TRUE}.}
 \item{apVar}{a logical value indicating whether the approximate
   covariance matrix of the variance-covariance parameters should be
   calculated.  Default is \code{TRUE}.}
 \item{.relStep}{relative step for numerical derivatives
   calculations.  Default is \code{.Machine$double.eps^(1/3)}.}
 \item{opt}{the optimizer to be used, either \code{"\link{nlminb}"} (the
   default) or \code{"\link{optim}"}.}
 \item{optimMethod}{character - the optimization method to be used with
   the \code{\link{optim}} optimizer.  The default is
   \code{"BFGS"}.  An alternative is \code{"L-BFGS-B"}.}
 \item{minAbsParApVar}{numeric value - minimum absolute parameter value
   in the approximate variance calculation.  The default is \code{0.05}.}
 \item{natural}{a logical value indicating whether the \code{pdNatural}
   parametrization should be used for general positive-definite matrices
   (\code{pdSymm}) in \code{reStruct}, when the approximate covariance
   matrix of the estimators is calculated.  Default is \code{TRUE}.}
 \item{sigma}{optionally a positive number to fix the residual error at.
   If \code{NULL}, as by default, or \code{0}, sigma is estimated.}
 \item{\dots}{further named control arguments to be passed, depending on
   \code{opt}, to \code{\link{nlminb}} (those from \code{abs.tol} down)
   or \code{\link{optim}} (those except \code{trace} and \code{maxit};
   \code{reltol} is used only from the second iteration).}
}
\value{
  a list with components for each of the possible arguments.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}; the
  \code{sigma} option: Siem Heisterkamp and Bert van Willigen.}
\seealso{\code{\link{lme}}, \code{\link{nlminb}}, \code{\link{optim}}}
\examples{
# decrease the maximum number iterations in the ms call and
# request that information on the evolution of the ms iterations be printed
str(lCtr <- lmeControl(msMaxIter = 20, msVerbose = TRUE))
## This should always work:
do.call(lmeControl, lCtr)
}
\keyword{models}
