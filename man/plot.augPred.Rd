% $Id: plot.augPred.Rd,v 1.4.2.1 2000/12/02 21:06:08 bates Exp $
\name{plot.augPred}
\title{Plot an augPred Object}
\usage{
plot(x, key, grid, ...)
}
\alias{plot.augPred}
\arguments{
 \item{x}{an object of class \code{augPred}.}
 \item{key}{an optional logical value, or list. If \code{TRUE}, a legend
 is included at the top of the plot indicating which symbols (colors)
 correspond to which prediction levels. If \code{FALSE}, no legend
 is included. If given as a list, \code{key} is passed down as an
 argument to the \code{trellis} function generating the plots
 (\code{xyplot}). Defaults to \code{TRUE}.}
 \item{grid}{an optional logical value indicating whether a grid should
   be added to plot. Default is \code{FALSE}.}
 \item{\dots}{optional arguments passed down to the \code{trellis}
   function generating the plots.}
}
\description{
  A Trellis \code{xyplot} of predictions versus the primary covariate is
  generated, with a different panel for each value of the grouping
  factor. Predicted values are joined by lines, with  
  different line types (colors) being used for each level of
  grouping. Original observations are represented by circles.
}
\value{
  A Trellis plot of predictions versus the primary covariate, with panels
  determined by the grouping factor.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{augPred}}, \code{\link{xyplot}}}

\examples{
\dontrun{
data(Orthodont)
fm1 <- lme(Orthodont)
plot(augPred(fm1, level = 0:1, length.out = 2))
}
}
\keyword{models}
