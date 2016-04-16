% File nlme/man/pdLogChol.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{pdLogChol}
\title{General Positive-Definite Matrix}
\usage{
pdLogChol(value, form, nam, data)
}
\alias{pdLogChol}
\arguments{
 \item{value}{an optional initialization value, which can be any of the
   following: a \code{pdMat} object, a positive-definite
   matrix, a one-sided linear formula (with variables separated by
   \code{+}), a vector of character strings, or a numeric
   vector.  Defaults to \code{numeric(0)}, corresponding to an
   uninitialized object.}
 \item{form}{an optional one-sided linear formula specifying the
   row/column names for the matrix represented by \code{object}.  Because
   factors may be present in \code{form}, the formula needs to be
   evaluated on a data frame to resolve the names it defines.  This
   argument is ignored when \code{value} is a one-sided
   formula.  Defaults to \code{NULL}.}
 \item{nam}{an optional character vector specifying the row/column names
   for the matrix represented by object.  It must have length equal to
   the dimension of the underlying positive-definite matrix and
   unreplicated elements.  This argument is ignored when
   \code{value} is a character vector.  Defaults to \code{NULL}.}
 \item{data}{an optional data frame in which to evaluate the variables
   named in \code{value} and \code{form}.  It is used to obtain the
   levels for \code{factors}, which affect the dimensions and the
   row/column names of the underlying matrix.  If \code{NULL}, no
   attempt is made to obtain information on \code{factor}s appearing in
   the formulas.  Defaults to the parent frame from which the function
   was called.}
}

\description{
  This function is a constructor for the \code{pdLogChol} class,
  representing a general positive-definite matrix.  If the matrix
  associated with \code{object} is of dimension \eqn{n}, it is
  represented by \eqn{n(n+1)/2}{n*(n+1)/2} unrestricted parameters,
  using the log-Cholesky parametrization described in Pinheiro and
  Bates (1996).
  \itemize{
    \item When \code{value} is \code{numeric(0)}, an uninitialized \code{pdMat}
    object, a one-sided formula, or a character vector, \code{object} is
    returned as an \emph{uninitialized} \code{pdLogChol} object (with just
    some of its attributes and its class defined) and needs to have its
    coefficients assigned later, generally using the \code{coef} or
    \code{matrix} replacement functions.
    \item If \code{value} is an \emph{initialized} \code{pdMat} object,
    \code{object} will be constructed from \code{as.matrix(value)}.
    \item Finally, if \code{value} is a numeric vector, it is assumed to
    represent the unrestricted coefficients of the matrix-logarithm
    parametrization of the underlying positive-definite matrix.
  }
}
\details{
  Internally, the \code{pdLogChol} representation of a symmetric
  positive definite matrix is a vector starting with the logarithms of
  the diagonal of the Choleski factorization of that matrix followed by
  its upper triangular portion.
}
\value{
  a \code{pdLogChol} object representing a general positive-definite
  matrix, also inheriting from class \code{pdMat}.
}
\references{
  Pinheiro, J.C. and Bates., D.M. (1996)
  Unconstrained Parametrizations for Variance-Covariance Matrices,
  \emph{Statistics and Computing} \bold{6}, 289--296.

  Pinheiro, J.C., and Bates, D.M. (2000)
  \emph{Mixed-Effects Models in S and S-PLUS}, Springer.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{as.matrix.pdMat}},
  \code{\link{coef.pdMat}},
  \code{\link{pdClasses}},
  \code{\link{matrix<-.pdMat}}}
\examples{
(pd1 <- pdLogChol(diag(1:3), nam = c("A","B","C")))

(pd4 <- pdLogChol(1:6))
(pd4c <- chol(pd4)) # -> upper-tri matrix with off-diagonals  4 5 6
pd4c[upper.tri(pd4c)]
log(diag(pd4c)) # 1 2 3
}
\keyword{models}
