% $Id: print.anova.lme.Rd,v 1.8 2002/03/05 14:59:39 bates Exp $
\name{print.anova.lme}
\title{Print an anova.lme Object}
\usage{
\method{print}{anova.lme}(x, verbose, \dots)
}
\alias{print.anova.lme}
\arguments{
 \item{x}{an object inheriting from class \code{anova.lme},
   generally obtained by applying the \code{anova.lme} method to an
   \code{lme} object.}
 \item{verbose}{an optional logical value. If \code{TRUE}, the calling
   sequences for each fitted model object are printed with the rest of
   the output, being omitted if \code{verbose = FALSE}. Defaults to
   \code{attr(x, "verbose")}.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  When only one fitted model object is used in the call to
  \code{anova.lme}, a data frame with the estimated values, the
  approximate standard errors, the z-ratios, and the  p-values for the
  fixed effects is printed. Otherwise, when multiple fitted objects are
  being compaerd, a data frame with the degrees of freedom, the
  (restricted) log-likelihood, the Akaike Information Criterion (AIC),
  and the Bayesian Information Criterion (BIC) of each fitted model
  object is printed. If included in \code{x}, likelihood ratio
  statistics, with associated p-values, are included in the output.   
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{anova.lme}}, \code{\link{lme}}, \code{\link[base]{print}}}
\examples{
data(Orthodont)
fm1 <- lme(distance ~ age * Sex, Orthodont, random = ~ age | Subject)
fm2 <- update(fm1, random = ~ 1 | Subject)
print(anova(fm1))                # single argument form
print(anova(fm1, fm2))           # multiple arguments
}
\keyword{models}
