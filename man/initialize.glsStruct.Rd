% $Id: initialize.glsStruct.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{initialize.glsStruct}
\title{Initialize a glsStruct Object}
\usage{
initialize(object, data, control)
}
\alias{initialize.glsStruct}
\arguments{
 \item{object}{an object inheriting from class \code{glsStruct},
   representing a list of linear model components, such as
   \code{corStruct} and \code{varFunc} objects.}
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{formula(object)}.}
 \item{control}{an optional list with control parameters for the
   initialization and optimization algorithms used in
   \code{gls}. Defaults to \code{list(singular.ok = FALSE, qrTol =
     .Machine$single.eps)}, implying that linear dependencies are not
   allowed in the model and that the tolerance for detecting linear
   dependencies among the columns of the regression matrix is
   \code{.Machine$single.eps}.}
}

\description{
  The individual linear model components of the \code{glsStruct} list
  are initialized.  
}
\value{
  a \code{glsStruct} object similar to \code{object}, but with
  initialized model components. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}, \code{\link{initialize.corStruct}} ,
  \code{\link{initialize.varFunc}}} 

\keyword{models}
