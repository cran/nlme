% File nlme/man/plot.intervals.lmList.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{plot.intervals.lmList}
\title{Plot lmList Confidence Intervals}
\usage{
\method{plot}{intervals.lmList}(x, xlab = "", ylab = attr(x, "groupsName"),
     strip = function(...) strip.default(..., style = 1),
     \dots)
}
\alias{plot.intervals.lmList}
\arguments{
 \item{x}{an object inheriting from class \code{"\link{intervals.lmList}"},
   representing confidence intervals and estimates for the coefficients
   in the \code{lm} components of the \code{lmList} object used to
   produce \code{x}.
 }
 \item{xlab, ylab}{axis labels, each with a sensible default.}
 \item{strip}{a \code{\link{function}} or \code{FALSE}, see
   \code{\link[lattice]{dotplot}()} from package \CRANpkg{lattice}.}
 \item{\dots}{optional arguments passed to the \code{dotplot} function
   (see above).}
}
\description{
  A Trellis dot-plot of the confidence intervals on the linear model
  coefficients is generated, with a different panel for each
  coefficient. Rows in the dot-plot correspond to the names of the
  \code{lm} components of the \code{lmList} object used to produce
  \code{x}. The lower and upper confidence limits are connected by
  a line segment and the estimated coefficients are marked with a
  \code{"+"}. 

  This is based on function \code{\link[lattice]{dotplot}()} from
  package \CRANpkg{lattice}.
}

\value{
  a Trellis plot with the confidence intervals on the coefficients of
  the individual \code{lm} components of the \code{lmList} that
  generated \code{x}.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{intervals.lmList}},
  \code{\link{lmList}},
  \code{\link{dotplot}}}
\examples{
fm1 <- lmList(distance ~ age | Subject, Orthodont)
plot(intervals(fm1))
}
\keyword{models}
