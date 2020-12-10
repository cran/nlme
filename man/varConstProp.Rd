\name{varConstProp}
\title{Constant Plus Proportion Variance Function}
\usage{
varConstProp(const, prop, form, fixed)
}
\alias{varConstProp}
\arguments{
  \item{const, prop}{optional numeric vectors, or lists of numeric
    values, with, respectively, the coefficients for the constant
    and the proportional error terms. Both arguments must have length one,
    unless a grouping factor is specified in \code{form}. If either argument
    has length greater than one, it must have names which identify its
    elements to the levels of the grouping factor defined in
    \code{form}. If a grouping factor is present in
    \code{form} and the argument has length one, its value will be
    assigned to all grouping levels. Only positive values are allowed
    for \code{const}. Default is 0.1 for both \code{const} and
    \code{prop}. }
  \item{form}{an optional one-sided formula of the form \code{~ v}, or
    \code{~ v | g}, specifying a variance covariate \code{v} and,
    optionally, a grouping factor \code{g} for the coefficients. The
    variance covariate must evaluate to a numeric vector and may involve
    expressions using \code{"."}, representing  a fitted model object
    from which fitted values (\code{fitted(.)}) and residuals
    (\code{resid(.)}) can be extracted (this allows the variance
    covariate to be updated during the optimization of an object
    function). When a grouping factor is present in \code{form},
    a different coefficient value is used for each of its
    levels. Several grouping variables may be
    simultaneously specified, separated by the \code{*} operator, as
    in \code{~ v | g1 * g2 * g3}. In this case, the levels of each
    grouping variable are pasted together and the resulting factor is
    used to group the observations. Defaults to \code{~ fitted(.)}
    representing a variance covariate given by the fitted values of a
    fitted model object and no grouping factor. }
  \item{fixed}{an optional list with components \code{const} and/or
    \code{power}, consisting of numeric vectors, or lists of numeric
    values, specifying the values at which some or all of the
    coefficients in the variance function should be fixed. If a grouping
    factor is specified in \code{form}, the components of \code{fixed}
    must have names identifying which coefficients are to be
    fixed. Coefficients included in \code{fixed} are not allowed to vary
    during the optimization of an objective function. Defaults to
    \code{NULL}, corresponding to no fixed coefficients.}
}
\description{
  This function is a constructor for the \code{varConstProp} class,
  representing a variance function structure corresponding to
  a two-component error model (additive and proportional error). Letting
  \eqn{v} denote the variance covariate and \eqn{\sigma^2(v)}{s2(v)}
  denote the variance function evaluated at \eqn{v}, the two-component variance
  function is defined as
  \eqn{\sigma^2(v) = a^2 + b^2 * v^{2}}{s2(v) = a^2 + b^2 * v^2)}, where a is
  the additive component and b is the relative error component.  In order
  to avoid overparameterisation of the model, it is recommended to use
  the possibility to fix sigma, preferably to a value of 1 (see examples).
}
\note{
  The error model underlying this variance function structure can be understood
  to result from two uncorrelated sequences of standardized random variables
  (Lavielle(2015), p. 55) and has been proposed for use in analytical chemistry
  (Werner et al. (1978), Wilson et al. (2004)) and chemical degradation
  kinetics (Ranke and Meinecke (2019)).  Note that the two-component error
  model proposed by Rocke and Lorenzato (1995) assumed a log-normal
  distribution of residuals at high absolute values, which is not
  compatible with the \code{\link{varFunc}} structures in package \pkg{nlme}.
}
\value{
  a \code{varConstProp} object representing a constant plus proportion variance
  function structure, also inheriting from class \code{varFunc}.
}

\references{
  Lavielle, M. (2015)
  \emph{Mixed Effects Models for the Population Approach: Models, Tasks,
    Methods and Tools}, Chapman and Hall/CRC.
  \doi{10.1201/b17203}

  Pinheiro, J.C., and Bates, D.M. (2000)
  \emph{Mixed-Effects Models in S and S-PLUS}, Springer.
  \doi{10.1007/b98882}

  Ranke, J., and Meinecke, S. (2019)
  Error Models for the Kinetic Evaluation of Chemical Degradation Data.
  \emph{Environments} \bold{6}(12), 124
  \doi{10.3390/environments6120124}

  Rocke, David M. and Lorenzato, Stefan (1995)
  A Two-Component Model for Measurement Error in Analytical Chemistry.
  \emph{Technometrics} \bold{37}(2), 176--184.
  \doi{10.1080/00401706.1995.10484302}

  Werner, Mario, Brooks, Samuel H., and Knott, Lancaster B. (1978)
  Additive, Multiplicative, and Mixed Analytical Errors.
  \emph{Clinical Chemistry} \bold{24}(11), 1895--1898.
  \doi{10.1093/clinchem/24.11.1895}

  Wilson, M.D., Rocke, D.M., Durbin, B. and Kahn, H.D (2004)
  Detection Limits and Goodness-of-Fit Measures for the Two-Component Model
  of Chemical Analytical Error.
  \emph{Analytica Chimica Acta} 2004, 509, 197--208
  \doi{10.1016/j.aca.2003.12.047}
}
\author{José Pinheiro and Douglas Bates (for \code{\link{varConstPower}}) and
  Johannes Ranke (adaptation to \code{varConstProp()}).
}
\seealso{
  \code{\link{varClasses}},
  \code{\link{varWeights.varFunc}},
  \code{\link{coef.varFunc}}}

\examples{
# Generate some synthetic data using the two-component error model and use
# different variance functions, also with fixed sigma in order to avoid
# overparameterisation in the case of a constant term in the variance function
times <- c(0, 1, 3, 7, 14, 28, 56, 120)
pred <- 100 * exp(- 0.03 * times)
sd_pred <- sqrt(3^2 + 0.07^2 * pred^2)
n_replicates <- 2

set.seed(123456)
syn_data <- data.frame(
  time = rep(times, each = n_replicates),
  value = rnorm(length(times) * n_replicates,
    rep(pred, each = n_replicates),
    rep(sd_pred, each = n_replicates)))
syn_data$value <- ifelse(syn_data$value < 0, NA, syn_data$value)

f_const <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
  data = syn_data, na.action = na.omit,
  start = list(parent_0 = 100, lrc = -3))
f_varPower <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
  data = syn_data, na.action = na.omit,
  start = list(parent_0 = 100, lrc = -3),
  weights = varPower())
f_varConstPower <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
  data = syn_data, na.action = na.omit,
  start = list(parent_0 = 100, lrc = -3),
  weights = varConstPower())
f_varConstPower_sf <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
  data = syn_data, na.action = na.omit,
  control = list(sigma = 1),
  start = list(parent_0 = 100, lrc = -3),
  weights = varConstPower())
f_varConstProp <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
  data = syn_data, na.action = na.omit,
  start = list(parent_0 = 100, lrc = -3),
  weights = varConstProp())
f_varConstProp_sf <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
  data = syn_data, na.action = na.omit,
  start = list(parent_0 = 100, lrc = -3),
  control = list(sigma = 1),
  weights = varConstProp())

AIC(f_const, f_varPower, f_varConstPower, f_varConstPower_sf,
  f_varConstProp, f_varConstProp_sf)

# The error model parameters 3 and 0.07 are approximately recovered
intervals(f_varConstProp_sf)
}
\keyword{models}
