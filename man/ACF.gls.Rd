% $Id: ACF.gls.Rd,v 1.4 2000/07/03 18:22:45 bates Exp $
\name{ACF.gls}
\title{Autocorrelation Function for gls Residuals}
\usage{
ACF(object, maxLag, resType, form, na.action)
}
\alias{ACF.gls}
\arguments{
  \item{object}{an object inheriting from class \code{gls}, representing
    a generalized least squares fitted  model.}
  \item{maxLag}{an optional integer giving the maximum lag for which the
    autocorrelation should be calculated. Defaults to maximum lag in the
    residuals.}  
  \item{resType}{an optional character string specifying the type of
    residuals to be used. If \code{"response"}, the "raw" residuals
    (observed - fitted) are used; else, if \code{"pearson"}, the
    standardized residuals (raw residuals divided by the corresponding
    standard errors) are used; else, if \code{"normalized"}, the
    normalized residuals (standardized residuals pre-multiplied by the
    inverse square-root factor of the estimated error correlation
    matrix) are used. Partial matching of arguments is used, so only the
    first character needs to be provided. Defaults to \code{"pearson"}.} 
  \item{form}{an optional one sided formula of the form \code{~ t}, or
    \code{~ t | g}, specifying a time covariate \code{t} and,  optionally, a
    grouping factor \code{g}. The time covariate must be integer
    valued. When a grouping factor is present in 
    \code{form}, the autocorrelations are calculated using residual pairs
    within the same group. Defaults to \code{~ 1}, which corresponds to
    using the order of the observations in the data as a covariate, and
    no groups.}   
  \item{na.action}{a function that indicates what should happen when the
    data contain \code{NA}s.  The default action (\code{na.fail}) causes
    \code{ACF.gls} to print an error message and terminate if there are any
    incomplete observations.}
}
\description{
  This method function calculates the empirical autocorrelation function
  for the residuals from an \code{gls} fit. If a grouping variable is
  specified in \code{form}, the autocorrelation values
  are calculated using pairs of residuals within the same group;
  otherwise all possible residual pairs are used. The autocorrelation
  function is useful for investigating serial correlation models for
  equally spaced data.  
}
\value{
  a data frame with columns \code{lag} and \code{ACF} representing,
  respectively, the lag between residuals within a pair and the corresponding
  empirical autocorrelation. The returned value inherits from class
  \code{ACF}.  
}
\references{
  Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
  Analysis: Forecasting and Control", 3rd Edition, Holden-Day. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{ACF.gls}}, \code{\link{plot.ACF}}}
\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary)
ACF(fm1, form = ~ 1 | Mare)
}
\keyword{models}
