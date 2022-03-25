% File nlme/man/gnlsControl.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{gnlsControl}
\alias{gnlsControl}
\title{Control Values for gnls Fit}
\usage{
gnlsControl(maxIter = 50, nlsMaxIter = 7, msMaxIter = 50, minScale = 0.001,
            tolerance = 1e-6, nlsTol = 0.001, msTol = 1e-7,
            returnObject = FALSE, msVerbose = FALSE,
            apVar = TRUE, .relStep =,
            opt = c("nlminb", "optim"), optimMethod = "BFGS",
            minAbsParApVar = 0.05, sigma = NULL)
}
\arguments{
 \item{maxIter}{maximum number of iterations for the \code{gnls}
   optimization algorithm.  Default is 50.}
 \item{nlsMaxIter}{maximum number of iterations
   for the \code{nls} optimization step \emph{inside} the \code{gnls}
   optimization.  Default is 7.}
 \item{msMaxIter}{maximum number of iterations
   for the \code{opt}imization step inside the \code{gnls}
   optimization.  Default is 50.}
 \item{minScale}{minimum factor by which to shrink the default step size
   in an attempt to decrease the sum of squares in the \code{nls} step.
   Default 0.001.}
 \item{tolerance}{tolerance for the convergence criterion in the
   \code{gnls} algorithm. Default is 1e-6.}
 \item{nlsTol}{tolerance for the convergence criterion in \code{nls}
   step. Default is 1e-3.}
 \item{msTol}{tolerance for the convergence criterion of the first outer
   iteration when \code{optim} is used.  Default is 1e-7.}
 \item{returnObject}{a logical value indicating whether the fitted
   object should be returned with a \code{\link{warning}} (instead of an
   error via \code{\link{stop}()}) when the maximum number of
   iterations is reached without convergence of the algorithm.}
 \item{msVerbose}{a logical value passed as the \code{trace} argument to
   the optimizer chosen by \code{opt}; see documentation on that.
   Default is \code{FALSE}.}
 \item{apVar}{a logical value indicating whether the approximate
   covariance matrix of the variance-covariance parameters should be
   calculated.  Default is \code{TRUE}.}
 \item{.relStep}{relative step for numerical derivatives
   calculations.  Default is \code{.Machine$double.eps^(1/3)} (about 6e-6).}
 \item{opt}{the optimizer to be used, either \code{"\link{nlminb}"} (the
   current default) or \code{"\link{optim}"} (the previous default).}
 \item{optimMethod}{character - the optimization method to be used with
   the \code{\link{optim}} optimizer. The default is
   \code{"BFGS"}.  An alternative is \code{"L-BFGS-B"}.}
 \item{minAbsParApVar}{numeric value - minimum absolute parameter value
   in the approximate variance calculation.  The default is \code{0.05}.}
 \item{sigma}{optionally a positive number to fix the residual error at.
   If \code{NULL}, as by default, or \code{0}, sigma is estimated.}
}
\description{
  The values supplied in the function call replace the defaults and a
  list with all possible arguments is returned. The returned list is
  used as the \code{control} argument to the \code{gnls} function.
}
\value{
  a list with components for each of the possible arguments.
}

\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}; the
  \code{sigma} option: Siem Heisterkamp and Bert van Willigen.}

\seealso{\code{\link{gnls}}}

\examples{
# decrease the maximum number of iterations and request tracing
gnlsControl(msMaxIter = 20, msVerbose = TRUE)
}
\keyword{models}
