% $Id: plot.compareFits.Rd,v 1.8 2002/03/05 14:59:39 bates Exp $
\name{plot.compareFits}
\title{Plot a compareFits Object}
\usage{
\method{plot}{compareFits}(x, subset, key, mark, ...)
}
\alias{plot.compareFits}
\arguments{
 \item{x}{an object of class \code{compareFits}.}
 \item{subset}{an optional logical or integer vector specifying which
   rows of \code{x} should be used in the plots. If missing, all
   rows are used.}
 \item{key}{an optional logical value, or list. If \code{TRUE}, a legend
 is included at the top of the plot indicating which symbols (colors)
 correspond to which objects being compared. If \code{FALSE}, no legend
 is included. If given as a list, \code{key} is passed down as an
 argument to the \code{trellis} function generating the plots
 (\code{dotplot}). Defaults to \code{TRUE}.}
 \item{mark}{an optional numeric vector, of length equal to the number of
   coefficients being compared, indicating where vertical lines should
   be drawn in the plots. If missing, no lines are drawn.} 
 \item{\dots}{optional arguments passed down to the \code{trellis}
   function generating the plots.}
}
\description{
  A Trellis \code{dotplot} of the values being compared, with different
  rows per group, is generated, with a different panel for each
  coefficient. Different symbols (colors) are used for each object being
  compared.
}
\value{
  A Trellis \code{dotplot} of the values being compared, with rows
  determined by the groups and panels by the coefficients.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{compareFits}}, \code{\link{pairs.compareFits}},
  \code{\link{dotplot}} }
\examples{
\dontrun{
data(Orthodont)
fm1 <- lmList(Orthodont)
fm2 <- lme(Orthodont)
plot(compareFits(coef(fm1), coef(fm2)))
}
}
\keyword{models}
