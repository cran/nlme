% $Id: reStruct.Rd,v 1.6 2002/03/05 14:59:39 bates Exp $
\name{reStruct}
\title{Random Effects Structure}
\usage{
reStruct(object, pdClass, REML, data)
}
\alias{reStruct}
\alias{[.reStruct}
\arguments{
 \item{object}{any of the following: (i) a one-sided formula of the form
   \code{~x1+...+xn | g1/.../gm}, with \code{x1+...+xn} specifying the
   model for the random effects and \code{g1/.../gm} the grouping
   structure (\code{m} may be equal to 1, in which case no \code{/} is
   required). The random effects formula will be repeated for all levels
   of grouping, in the case of multiple levels of grouping; (ii) a list of
   one-sided formulas of the form \code{~x1+...+xn | g}, with possibly
   different random effects models for each grouping level. The order of
   nesting will be assumed the same as the order of the elements in the
   list; (iii) a one-sided formula of the form \code{~x1+...+xn}, or a
   \code{pdMat} object with a formula (i.e. a non-\code{NULL} value for
   \code{formula(object)}), or a list of such formulas or \code{pdMat}
   objects. In this case, the grouping structure formula will be derived
   from the data used to to fit the mixed-effects model, which should
   inherit from class \code{groupedData}; (iv) a named list of formulas or
   \code{pdMat} objects as in (iii), with the grouping factors as
   names. The order of nesting will be assumed the same as the order of
   the order of the elements in the list; (v) an \code{reStruct} object.}
 \item{pdClass}{an optional character string with the name of the
   \code{pdMat} class to be used for the formulas in
   \code{object}. Defaults to \code{"pdSymm"} which corresponds to a
   general positive-definite matrix.}
 \item{REML}{an optional logical value. If \code{TRUE}, the associated
   mixed-effects model will be fitted using restricted maximum
   likelihood; else, if \code{FALSE}, maximum likelihood will be
   used. Defaults to \code{FALSE}.}
 \item{data}{an optional data frame in which to evaluate the variables
   used in the random effects formulas in \code{object}. It is used to
   obtain the levels for \code{factors}, which affect the dimensions and
   the row/column names of the underlying \code{pdMat} objects. If
   \code{NULL}, no attempt is made to obtain information on
   \code{factors} appearing in the formulas. Defaults to the
   parent frame from which the function was called.}
}
\description{
  This function is a constructor for the \code{reStruct} class,
  representing a random effects structure and consisting of a list of
  \code{pdMat} objects, plus a \code{settings} attribute containing
  information for the optimization algorithm used to fit the associated
  mixed-effects model.
}
\value{
  an object inheriting from class \code{reStruct}, representing a random
  effects structure.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{pdMat}}, \code{\link{lme}},
  \code{\link{groupedData}}} 

\examples{
data(Pixel)
rs1 <- reStruct(list(Dog = ~day, Side = ~1), data = Pixel)
rs1
}
\keyword{models}
