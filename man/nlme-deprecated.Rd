% File nlme/man/nlme-deprecated.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{nlme-deprecated}
\alias{nlme-deprecated}
\title{Deprecated Functions in Package \pkg{nlme}}

%% deprecation announced in 2018-08, with warning since 2023-01
\alias{nfGroupedData}
\alias{nmGroupedData}
%% deprecated 2023-01:
\alias{corIdent}

\usage{
## internal functions wrapped by groupedData()
nfGroupedData(formula, data, order.groups, FUN, outer, inner,
              labels, units)
nmGroupedData(formula, data, order.groups, FUN, outer, inner,
              labels, units)

## trivial "identity" correlation, never usefully implemented:
corIdent(form = NULL)
}

\description{
  These functions are provided for compatibility with older versions of
  \pkg{nlme} only, and may be defunct as soon as the next release.
}

\examples{
assertDeprecation <- function(expr)
  tools::assertCondition(expr, "deprecatedWarning", verbose = TRUE)

assertDeprecation(
  nlme::nfGroupedData(height ~ age | Subject, as.data.frame(Oxboys))
)
assertDeprecation( csId <-  corIdent(~ 1 | Subject) )
assertDeprecation( csI. <- Initialize(csId, data = Orthodont) )
assertDeprecation( corMatrix(csI.) )  # actually errors
}

\keyword{internal}
