\name{nlme-internal}
\alias{allVarsRec}
\alias{createConLin}
\alias{getFixDF}
\alias{getParsGnls}
\alias{getParsNlme}
\alias{glsApVar}
\alias{glsEstimate}
\alias{gnlsApVar}
\alias{lmeApVar}
\alias{nonlinModel}
\alias{print.correlation}
\title{Internal nlme functions}
\description{
  Internal nlme functions.
}
\usage{
allVarsRec(object)
createConLin(fixed,data,random,\dots)
getFixDF(X, grps, ngrps, assign, terms)
getParsGnls(plist, pmap, beta, N)
getParsNlme(plist, fmap, rmapRel, bmap, groups, beta, bvec, b, level, N)
glsApVar(glsSt, sigma, conLin, .relStep, minAbsPar, natural)
glsEstimate(object, conLin, control)
gnlsApVar(gnlsSt, lsigma, conLin, .relStep, minAbsPar, natural)
lmeApVar(lmeSt, sigma, conLin, .relStep, minAbsPar, natural)
nonlinModel(modelExpression, env, paramNames)
\method{print}{correlation}(x, title, rdig, \dots)
}
\details{
  These are not to be called by the user.
}
\keyword{internal}
