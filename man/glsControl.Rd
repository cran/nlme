% $Id: glsControl.Rd,v 1.5 2000/07/03 18:22:47 bates Exp $
\name{glsControl}
\title{Control Values for gls Fit}
\usage{
glsControl(maxIter, msMaxIter, tolerance, msTol, msScale, msVerbose,
           singular.ok, qrTol, returnObject, apVar, nlmStepMax, .relStep) 
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
 \item{msTol}{tolerance for the convergence criterion in \code{ms},
   passed as the \code{rel.tolerance} argument to the function (see
   documentation on \code{ms}). Default is 1e-7.}
 \item{msScale}{scale function passed as the \code{scale} argument to
   the \code{ms} function (see documentation on that function). Default
   is \code{lmeScale}.}
 \item{msVerbose}{a logical value passed as the \code{trace} argument to
   \code{ms} (see documentation on that function). Default is
   \code{FALSE}.}
 \item{singular.ok}{a logical value indicating whether non-estimable
   coefficients (resulting from linear dependencies among the columns of
   the regression matrix) should be allowed. Default is \code{FALSE}.}
 \item{qrTol}{a tolerance for detecting linear dependencies among the
   columns of the regression matrix in its QR decomposition. Default is
   \code{.Machine$single.eps}.}
 \item{returnObject}{a logical value indicating whether the fitted
   object should be returned when the maximum number of iterations is
   reached without convergence of the algorithm. Default is
   \code{FALSE}.} 
 \item{apVar}{a logical value indicating whether the approximate
   covariance matrix of the variance-covariance parameters should be
   calculated. Default is \code{TRUE}.}
 \item{nlmStepMax}{stepmax value to be passed to nlm. See
     \code{\link{nlm}} for details. Default is 100.0}
 \item{.relStep}{relative step for numerical derivatives
   calculations. Default is \code{.Machine$double.eps^(1/3)}.}
}
\description{
  The values supplied in the function call replace the defaults and a
  list with all possible arguments is returned. The returned list is
  used as the \code{control} argument to the \code{gls} function.
}
\value{
  a list with components for each of the possible arguments.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}, \code{\link{lmeScale}}}

\examples{
# decrease the maximum number iterations in the optimization call and
# request that information on the evolution of the ms iterations be printed
glsControl(msMaxIter = 20, msVerbose = TRUE)

}
\keyword{models}
