% $Id: plot.intervals.lmList.Rd,v 1.4.2.1 2000/12/02 21:06:08 bates Exp $
\name{plot.intervals.lmList}
\title{Plot lmList Confidence Intervals}
\usage{
plot(object, ...)
}
\alias{plot.intervals.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{intervals.lmList},
   representing confidence intervals and estimates for the coefficients
   in the \code{lm} components of the \code{lmList} object used to
   produce \code{object}.
 }
 \item{\dots}{optional arguments passed to the Trellis \code{dotplot}
   function.
 } 
}
\description{
  A Trellis dot-plot of the confidence intervals on the linear model
  coefficients is generated, with a different panel for each
  coefficient. Rows in the dot-plot correspond to the names of the
  \code{lm} components of the \code{lmList} object used to produce
  \code{object}. The lower and upper confidence limits are connected by
  a line segment and the estimated coefficients are marked with a
  \code{"+"}. The Trellis function \code{dotplot} is used in this method
  function.  
}

\value{
  a Trellis plot with the confidence intervals on the coefficients of
  the individual \code{lm} components of the \code{lmList} that
  generated \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{intervals.lmList}}, \code{\link{lmList}},
  \code{\link{dotplot}}} 
\examples{
\dontrun{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
plot(intervals(fm1))
}
}
\keyword{models}
