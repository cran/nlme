% $Id: lme.groupedData.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{lme.groupedData}
\title{LME fit from groupedData Object}
\usage{
lme(fixed, data, random, correlation, weights, subset, method,
    na.action, control)
}
\alias{lme.groupedData}
\arguments{
  \item{fixed}{a data frame inheriting from class \code{groupedData}.}
  \item{data}{this argument is included for consistency with the generic
    function. It is ignored in this method function.}
 \item{other arguments}{identical to the arguments in the generic
   function call. See the documentation on \code{\link{lme}}.}
}
\description{
  The response variable and primary covariate in \code{formula(fixed)}
  are used to construct the fixed effects model formula. This formula
  and the   \code{groupedData} object are passed as the \code{fixed} and
  \code{data} arguments to \code{lme.formula}, together with any other
  additional arguments in the function call. See the documentation on
  \code{lme.formula} for a description of that function.
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
 Pinheiro, J.C. and Bates., D.M.  (1996).   The different correlation
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

\seealso{\code{\link{lme}}, \code{\link{groupedData}},
  \code{\link{lmeObject}} } 

\examples{
data(Orthodont)
fm1 <- lme(Orthodont)
}
\keyword{models}
