% $Id: corAR1.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{corAR1}
\title{AR(1) Correlation Structure}
\usage{
corAR1(value, form, fixed)
}
\alias{corAR1}
\arguments{
 \item{value}{the value of the lag 1 autocorrelation, which must be
   between -1 and 1. Defaults to 0 (no autocorrelation).}
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
  This function is a constructor for the \code{corAR1} class,
  representing an autocorrelation structure of order 1. Objects
  created using this constructor must later be initialized using the
  appropriate \code{initialize} method. 
}
\value{
  an object of class \code{corAR1}, representing an autocorrelation
   structure of order 1. 
}
\references{
  Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
  Analysis: Forecasting and Control", 3rd Edition, Holden-Day. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{initialize.corStruct}}}
\examples{
## covariate is observation order and grouping factor is Mare
cs1 <- corAR1(0.2, form = ~ 1 | Mare)
}
\keyword{models}
