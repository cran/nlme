% $Id: lme.lmList.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
\name{lme.lmList}
\title{LME fit from lmList Object}
\usage{
lme(fixed, data, random, correlation, weights, subset, method,
    na.action, control)
}
\alias{lme.lmList}
\arguments{
  \item{fixed}{an object inheriting from class \code{lmList},
    representing a list of \code{lm} fits with a common model.}
  \item{data}{this argument is included for consistency with the generic
    function. It is ignored in this method function.}
 \item{random}{an optional one-sided linear formula with no conditioning
   expression, or a \code{pdMat} object with a \code{formula}
   attribute. Multiple levels of grouping are not allowed with this
   method function.  Defaults to a formula consisting of the right hand
   side of \code{formula(fixed)}.} 
 \item{other arguments}{identical to the arguments in the generic
   function call. See the documentation on \code{\link{lme}}.}
}
\description{
  If the random effects names defined in \code{random} are a subset of
  the \code{lmList} object coefficient names, initial estimates for the
  covariance matrix of the random effects are obtained (overwriting any
  values given in \code{random}). \code{formula(fixed)} and the
  \code{data} argument in the calling sequence used to obtain
  \code{fixed} are passed as the \code{fixed} and \code{data} arguments
  to \code{lme.formula}, together with any other additional arguments in
  the function call. See the documentation on \code{lme.formula} for a
  description of that function. 
}
\value{
  an object of class \code{lme} representing the linear mixed-effects
  model fit. Generic functions such as \code{print}, \code{plot} and
  \code{summary} have methods to show the results of the fit. See
  \code{lmeObject} for the components of the fit. The functions
  \code{resid}, \code{coef}, \code{fitted}, \code{fixed.effects}, and
  \code{random.effects}  can be used to extract some of its components.
}
\references{
 The computational methods are described in Bates, D.M. and Pinheiro
 (1998) and follow on the general framework of Lindstrom, M.J. and Bates,
 D.M. (1988). The model formulation is described in Laird, N.M. and Ware,
 J.H. (1982).  The variance-covariance parametrizations are described in
 <Pinheiro, J.C. and Bates., D.M.  (1996).   The different correlation
 structures available for the \code{correlation} argument are described
 in Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994), Littel, R.C.,
 Milliken, G.A., Stroup, W.W., and Wolfinger, R.D. (1996), and Venables,
 W.N. and Ripley, B.D. (1997). The use of variance functions for linear
 and nonlinear mixed effects models is presented in detail in Davidian,
 M. and Giltinan, D.M. (1995). 

 Bates, D.M. and Pinheiro, J.C. (1998) "Computational methods for
 multilevel models" available in PostScript or PDF formats at
 http://franz.stat.wisc.edu/pub/NLME/
 
 Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
 Analysis: Forecasting and Control", 3rd Edition, Holden-Day. 

 Davidian, M. and Giltinan, D.M. (1995) "Nonlinear Mixed Effects Models
 for Repeated Measurement Data", Chapman and Hall.

 Laird, N.M. and Ware, J.H. (1982) "Random-Effects Models for
 Longitudinal Data", Biometrics, 38, 963-974.  

 Lindstrom, M.J. and Bates, D.M. (1988) "Newton-Raphson and EM
 Algorithms for Linear Mixed-Effects Models for Repeated-Measures
 Data", Journal of the American Statistical Association, 83,
 1014-1022. 

 Littel, R.C., Milliken, G.A., Stroup, W.W., and Wolfinger, R.D. (1996)
 "SAS Systems for Mixed Models", SAS Institute.

 Pinheiro, J.C. and Bates., D.M.  (1996) "Unconstrained
 Parametrizations for Variance-Covariance Matrices", Statistics and
 Computing, 6, 289-296.

 Venables, W.N. and Ripley, B.D. (1997) "Modern Applied Statistics with
 S-plus", 2nd Edition, Springer-Verlag.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}, \code{\link{lmList}},
  \code{\link{lmeObject}} } 

\examples{
data(Orthodont)
fm1 <- lmList(Orthodont)
fm2 <- lme(fm1)
}
\keyword{models}
