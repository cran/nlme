% $Id: corCompSymm.Rd,v 1.5 2002/03/05 14:59:39 bates Exp $
\name{corCompSymm}
\title{Compound Symmetry Correlation Structure}
\usage{
corCompSymm(value, form, fixed)
}
\alias{corCompSymm}
\arguments{
 \item{value}{the correlation between any two correlated
   observations. Defaults to 0.} 
 \item{form}{a one sided formula of the form \code{~ t}, or \code{~ t |
     g}, specifying a time covariate \code{t} and,  optionally, a
   grouping factor \code{g}. When a grouping factor is present in
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
  This function is a constructor for the \code{corCompSymm} class,
  representing a compound symmetry structure corresponding to uniform
  correlation. Objects created using this constructor must later be
  initialized using the appropriate \code{Initialize} method. 
}
\value{
  an object of class \code{corCompSymm}, representing a compound
  symmetry correlation structure. 
}
\references{
  Milliken, G. A. and Johnson, D. E. (1992) "Analysis of Messy Data,
  Volume I: Designed Experiments", Van Nostrand Reinhold.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{Initialize.corStruct}}}
\examples{
## covariate is observation order and grouping factor is Subject
cs1 <- corCompSymm(0.5, form = ~ 1 | Subject)
}
\keyword{models}

