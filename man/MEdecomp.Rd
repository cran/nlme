\name{MEdecomp}
\alias{MEdecomp}
\title{Decompose a mixed-effects model}
\description{
  Decompose the model matrices and response from a mixed-effects model
  to provide a more compact representation.
}
\usage{
MEdecomp(conLin)
}
\arguments{
  \item{conLin}{a condensed linear model object representing the model
    matrices and response from a mixed-effects model.}
}
\value{
  another condensed linear model object that provides an equivalent
  model fit with fewer rows.
}
\author{Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{MEestimate}}}
\keyword{models}
