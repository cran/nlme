% File nlme/man/gnlsControl.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{gnlsControl}
\alias{gnlsControl}
\title{Control Values for gnls Fit}
\usage{
gnlsControl(maxIter, nlsMaxIter, msMaxIter, minScale, tolerance,
            nlsTol, msTol, returnObject, msVerbose,
            apVar, .relStep,
            opt = c("nlminb", "optim"), optimMethod,
            minAbsParApVar)
}
\arguments{
 \item{maxIter}{maximum number of iterations for the \code{gnls}
   optimization algorithm. Default is 50.} 
 \item{nlsMaxIter}{maximum number of iterations
   for the \code{nls} optimization step inside the \code{gnls}
   optimization. Default is 7.}
 \item{msMaxIter}{maximum number of iterations
   for the \code{ms} optimization step inside the \code{gnls}
   optimization. Default is 50.}
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
   object should be returned when the maximum number of iterations is
   reached without convergence of the algorithm.  Default is
   \code{FALSE}.} 
 \item{msVerbose}{a logical value passed as the \code{trace} argument to
   \code{ms} (see documentation on that function). Default is
   \code{FALSE}.}
 \item{apVar}{a logical value indicating whether the approximate
   covariance matrix of the variance-covariance parameters should be
   calculated. Default is \code{TRUE}.}
 \item{.relStep}{relative step for numerical derivatives
   calculations. Default is \code{.Machine$double.eps^(1/3)}.}
 \item{opt}{the optimizer to be used, either \code{"\link{nlminb}"} (the
   current default) or \code{"\link{optim}"} (the previous default).}
 \item{optimMethod}{character - the optimization method to be used with
   the \code{\link{optim}} optimizer. The default is
   \code{"BFGS"}.  An alternative is \code{"L-BFGS-B"}.}
 \item{minAbsParApVar}{numeric value - minimum absolute parameter value
   in the approximate variance calculation.  The default is \code{0.05}.}
}
\description{
  The values supplied in the function call replace the defaults and a
  list with all possible arguments is returned. The returned list is
  used as the \code{control} argument to the \code{gnls} function.
}
\value{
  a list with components for each of the possible arguments.
}

\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gnls}}}

\examples{
# decrease the maximum number iterations in the ms call and
# request that information on the evolution of the ms iterations be printed
gnlsControl(msMaxIter = 20, msVerbose = TRUE)
}
\keyword{models}
