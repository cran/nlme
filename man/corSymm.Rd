% $Id: corSymm.Rd,v 1.4 2000/07/03 18:22:46 bates Exp $
\name{corSymm}
\title{General Correlation Structure}
\usage{
corSymm(value, form, fixed)
}
\alias{corSymm}
\arguments{
 \item{value}{an optional vector with the parameter values. Default is
   \code{numeric(0)}, which results in a vector of zeros of appropriate
   dimension being assigned to the parameters when \code{object} is
   initialized (corresponding to an identity correlation structure).} 
 \item{form}{a one sided formula of the form \code{~ t}, or \code{~ t |
     g}, specifying a time covariate \code{t} and,  optionally, a
   grouping factor \code{g}. A covariate for this correlation structure
   must be integer valued. When a grouping factor is present in
   \code{form}, the correlation structure is assumed to apply only
   to observations within the same grouping level; observations with
   different grouping levels are assumed to be uncorrelated. Defaults to
   \code{~ 1}, which corresponds to using the order of the observations
   in the data as a covariate, and no groups.}  
 \item{fixed}{an optional logical value indicating whether the
   coefficients should be allowed to vary in the optimization, or kept
   fixed at their initial value. Defaults to \code{FALSE}, in which case
   the coefficients are allowed to vary.}
}
\description{
  This function is a constructor for the \code{corSymm} class,
  representing a general correlation structure. The internal
  representation of this structure, in terms of unconstrained
  parameters, uses the spherical parametrization defined in Pinheiro and
  Bates (1996).  Objects created using this constructor must later be
  initialized using the  appropriate \code{initialize} method. 
}
\value{
  an object of class \code{corSymm} representing a general correlation
  structure.
}
\references{ 
Pinheiro, J.C. and Bates., D.M.  (1996) "Unconstrained
Parametrizations for Variance-Covariance Matrices", Statistics and
Computing, 6, 289-296.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{initialize.corSymm}}}

\examples{
## covariate is observation order and grouping factor is Subject
cs1 <- corSymm(form = ~ 1 | Subject)
}
\keyword{models}
