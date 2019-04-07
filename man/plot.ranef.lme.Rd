% File nlme/man/plot.ranef.lme.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{plot.ranef.lme}
\title{Plot a ranef.lme Object}
\description{
  Plots (class \code{"Trellis"} from package \CRANpkg{lattice}) of the
  random effects from linear mixed effects model, i.e., the result of
  \code{\link{ranef}(\link{lme}(*))} (of class \code{"\link{ranef.lme}"}).
}
\usage{
\method{plot}{ranef.lme}(x, form = NULL, omitFixed = TRUE, level = Q,
     grid = TRUE, control, xlab, ylab, strip,
     \dots)
}
\alias{plot.ranef.lme}
\arguments{
 \item{x}{an object inheriting from class
   \code{"\link{ranef.lme}"}, representing the estimated coefficients or
   estimated random effects for the \code{lme} object from which it was
   produced.}
 \item{form}{an optional formula specifying the desired type of plot.
   \itemize{
     \item If given as a one-sided formula, a \code{\link[lattice]{dotplot}()}
     of the estimated random effects (coefficients) grouped according to
     all combinations of the levels of the factors named in \code{form}
     is returned.

     \item If given as a two-sided formula (or by default, \code{NULL}),
     an \code{\link[lattice]{xyplot}()} Trellis display of the random
     effect (coefficient) versus the named covariates is returned.  In
     \code{NULL} case the row names of the random effects (coefficients)
     are used (as covariates).
   }
   See also \sQuote{Details:}.
 }
 \item{omitFixed}{an optional logical value indicating whether
   columns with values that are constant across groups should be
   omitted. Default is \code{TRUE}.}
 \item{level}{an optional integer value giving the level of grouping
   to be used for \code{x}. Only used when \code{x} is a list
   with different components for each grouping level. Defaults to the
   highest or innermost level of grouping.}
 \item{grid}{an optional logical value indicating whether a grid should
   be added to plot.  Only applies to plots associated with two-sided
   formulas in \code{form}.  Default is \code{TRUE}.}
 \item{control}{an optional list with control values for the
     plot, when \code{form} is given as a two-sided formula. The control
     values are referenced by name in the \code{control} list and  only
     the ones to be modified from the default need to be
     specified. Available values include: \code{drawLine}, a logical
     value indicating whether a \code{loess} smoother should be added to
     the scatter plots and a line connecting the medians should be added
     to the boxplots (default is \code{TRUE}); \code{span.loess}, used
     as the \code{span} argument in the call to \code{panel.loess}
     (default is \code{2/3}); \code{degree.loess}, used as the
     \code{degree} argument in the call to \code{panel.loess} (default
     is \code{1}); \code{cex.axis}, the character expansion factor for
     the x-axis (default is \code{0.8}); \code{srt.axis}, the rotation
     factor for the x-axis (default is \code{0}); and \code{mgp.axis}, the
     margin parameters for the x-axis (default is \code{c(2, 0.5, 0)}).}
 \item{xlab, ylab}{axis labels, each with a sensible default.}
 \item{strip}{a \code{\link{function}} or \code{FALSE}, see
   \code{\link[lattice]{dotplot}()} from package \CRANpkg{lattice}.}
 \item{\dots}{optional arguments passed to the Trellis \code{dotplot}
   function.}
}
\details{
  If \code{form} is missing, or is given as a one-sided formula, a
  Trellis dot-plot (via \code{\link[lattice]{dotplot}()} from pkg
  \CRANpkg{lattice}) of the random effects is generated, with a different
  panel for each random effect (coefficient).  Rows in the dot-plot are
  determined by the \code{form} argument (if not missing) or by the row
  names of the random effects (coefficients).  Single factors (\code{~g})
  or crossed factors (\code{~g1*g2}) are allowed.  For a single factor,
  its levels determine the dot-plot rows (with possibly multiple dots
  per row); otherwise, if \code{form}
  specifies a crossing of factors, the dot-plot rows are determined by
  all combinations of the levels of the individual factors in the
  formula.

  If \code{form} is a two-sided formula, the left hand side must
  be a single random effect (coefficient) and the right hand side is
  formed by covariates in \code{x} separated by \code{+}.
  An \code{\link[lattice]{xyplot}()} Trellis display is generated,
  with a different panel for each variable listed in the right hand side
  of \code{form}.  Scatter plots are generated for numeric variables and
  boxplots are generated for categorical (\code{factor} or
  \code{ordered}) variables.
}
\value{
  a Trellis plot of the estimated random-effects (coefficients) versus
  covariates, or groups.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{ranef.lme}}, \code{\link{lme}},
  \code{\link{dotplot}}.
}
\examples{
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
plot(ranef(fm1))
fm1RE <- ranef(fm1, aug = TRUE)
plot(fm1RE, form = ~ Sex)
plot(fm1RE, form = age ~ Sex) # "connected" boxplots
}
\keyword{models}
