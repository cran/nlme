% File nlme/man/nlmeControl.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{nlmeControl}
\title{Control Values for nlme Fit}
\usage{
nlmeControl(maxIter, pnlsMaxIter, msMaxIter, minScale,
            tolerance, niterEM, pnlsTol, msTol,
            returnObject, msVerbose, msWarnNoConv,
            gradHess, apVar, .relStep, minAbsParApVar = 0.05,
            opt = c("nlminb", "nlm"), natural = TRUE, sigma = NULL, \dots)
}
\alias{nlmeControl}
\arguments{
 \item{maxIter}{maximum number of iterations for the \code{nlme}
   optimization algorithm.  Default is 50.}
 \item{pnlsMaxIter}{maximum number of iterations
   for the \code{PNLS} optimization step inside the \code{nlme}
   optimization.  Default is 7.}
 \item{msMaxIter}{maximum number of iterations for \code{\link{nlminb}}
   (\code{iter.max}) or the \code{\link{nlm}} (\code{iterlim}, from the
   10-th step) optimization step inside the \code{nlme}
   optimization.  Default is 50 (which may be too small for e.g. for
   overparametrized cases).}
 \item{minScale}{minimum factor by which to shrink the default step size
   in an attempt to decrease the sum of squares in the \code{PNLS} step.
   Default \code{0.001}.}
 \item{tolerance}{tolerance for the convergence criterion in the
   \code{nlme} algorithm.  Default is \code{1e-6}.}
 \item{niterEM}{number of iterations for the EM algorithm used to refine
   the initial estimates of the random effects variance-covariance
   coefficients.  Default is 25.}
 \item{pnlsTol}{tolerance for the convergence criterion in \code{PNLS}
   step.  Default is \code{1e-3}.}
 \item{msTol}{tolerance for the convergence criterion in \code{nlm},
   passed as the \code{gradtol} argument to the function (see
   documentation on \code{nlm}).  Default is \code{1e-7}. }
 \item{returnObject}{a logical value indicating whether the fitted
   object should be returned when the maximum number of iterations is
   reached without convergence of the algorithm.  Default is
   \code{FALSE}.}
 \item{msVerbose}{a logical value passed as the \code{trace} to
   \code{\link{nlminb}(.., control= list(trace = *, ..))} or
   as argument \code{print.level} to \code{\link{nlm}()}.  Default is
   \code{FALSE}.}
 \item{msWarnNoConv}{logical indicating if a \code{\link{warning}}
   should be signalled whenever the minimization (by \code{opt}) in the
   LME step does not converge; defaults to \code{TRUE}.}
 \item{gradHess}{a logical value indicating whether numerical gradient
   vectors and Hessian matrices of the log-likelihood function should
   be used in the \code{nlm} optimization. This option is only available
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
 \item{minAbsParApVar}{numeric value - minimum absolute parameter value
   in the approximate variance calculation.  The default is \code{0.05}.}
 \item{opt}{the optimizer to be used, either \code{"\link{nlminb}"} (the
   default) or \code{"\link{nlm}"}.}
 \item{natural}{a logical value indicating whether the \code{pdNatural}
   parametrization should be used for general positive-definite matrices
   (\code{pdSymm}) in \code{reStruct}, when the approximate covariance
   matrix of the estimators is calculated.  Default is \code{TRUE}.}
 \item{sigma}{optionally a positive number to fix the residual error at.
   If \code{NULL}, as by default, or \code{0}, sigma is estimated.}
 \item{\dots}{Further, named control arguments to be passed to
   \code{\link{nlminb}} (apart from \code{trace} and \code{iter.max}
   mentioned above), where used (\code{eval.max} and those from
   \code{abs.tol} down).}
}
\description{
  The values supplied in the function call replace the defaults and a
  list with all possible arguments is returned.  The returned list is
  used as the \code{control} argument to the \code{nlme} function.
}
\value{
  a list with components for each of the possible arguments.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}; the
  \code{sigma} option: Siem Heisterkamp and Bert van Willigen.}
\seealso{\code{\link{nlme}}, \code{\link{nlm}}, \code{\link{optim}},
  \code{\link{nlmeStruct}}}
\examples{
# decrease the maximum number of iterations and request tracing
nlmeControl(msMaxIter = 20, msVerbose = TRUE)
}
\keyword{models}
