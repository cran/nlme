% $Id: corARMA.Rd,v 1.4 2000/07/03 18:22:46 bates Exp $
\name{corARMA}
\title{ARMA(p,q) Correlation Structure}
\usage{
corARMA(value, form, p, q, fixed)
}
\alias{corARMA}
\alias{coef.corARMA}
\arguments{
  \item{value}{a vector with the values of the autoregressive and moving
   average parameters, which must have length \code{p + q} and all
   elements between -1 and 1. Defaults to a vector of zeros,
   corresponding to uncorrelated observations.}
 \item{form}{a one sided formula of the form \code{~ t}, or \code{~ t |
     g}, specifying a time covariate \code{t} and,  optionally, a
   grouping factor \code{g}. A covariate for this correlation structure
   must be integer valued. When a grouping factor is present in
   \code{form}, the correlation structure is assumed to apply only
   to observations within the same grouping level; observations with
   different grouping levels are assumed to be uncorrelated. Defaults to
   \code{~ 1}, which corresponds to using the order of the observations
   in the data as a covariate, and no groups.}  
 \item{p, q}{non-negative integers specifying respectively the
   autoregressive order and the moving average order of the \code{ARMA}
   structure. Both default to 0.}
 \item{fixed}{an optional logical value indicating whether the
   coefficients should be allowed to vary in the optimization, or kept
   fixed at their initial value. Defaults to \code{FALSE}, in which case
   the coefficients are allowed to vary.}
} 
\description{
  This function is a constructor for the \code{corARMA} class,
  representing an autocorrelation-moving average correlation structure
  of order (p, q). Objects created using this constructor must later
  be initialized using the appropriate \code{initialize} method. 
}
\value{
  an object of class \code{corARMA}, representing an
  autocorrelation-moving average correlation structure. 
}

\references{
  Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
  Analysis: Forecasting and Control", 3rd Edition, Holden-Day. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{initialize.corStruct}}}
\examples{
## ARMA(1,2) structure, with observation order as a covariate and
## Mare as grouping factor
cs1 <- corARMA(c(0.2, 0.3, -0.1), form = ~ 1 | Mare, p = 1, q = 2)
}
\keyword{models}
