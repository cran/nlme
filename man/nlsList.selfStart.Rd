% $Id: nlsList.selfStart.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{nlsList.selfStart}
\title{nlsList Fit from a selfStart Function}
\usage{
nlsList(model, data, start, control, level, na.action, pool)
}
\alias{nlsList.selfStart}
\arguments{
  \item{model}{a \code{selfStart} model function, which calculates
    initial estimates for the model parameters from \code{data}.}
  \item{data}{a data frame in which to interpret the variables in
    \code{model}. Because no grouping factor can be specified in
    \code{model}, \code{data} must inherit from class
    \code{groupedData}. 
  }
 \item{other arguments}{identical to the arguments in the generic
   function call. See the documentation on \code{\link{nlsList}}.
 }
}
\description{
  The response variable and primary covariate in \code{formula(data)}
  are used together with \code{model} to construct the nonlinear model
  formula. This is used in the \code{nls} calls and, because a
  selfStarting model function can calculate initial estimates for its
  parameters from the data, no starting estimates need to be provided.
}
\value{
  a list of \code{nls} objects with as many components as the number of
  groups defined by the grouping factor. A \code{NULL} value is assigned
  to the components corresponding to clusters for which the \code{nls}
  algorithm failed to converge. Generic functions such as \code{coef},
  \code{fixed.effects}, \code{lme}, \code{pairs}, \code{plot},
  \code{predict}, \code{random.effects}, \code{summary}, and
  \code{update} have methods that can be applied to an \code{nlsList}
  object.  
}

\seealso{
  \code{\link{selfStart}}, \code{\link{groupedData}}, \code{\link{nls}},
  \code{\link{nlme.nlsList}}, \code{\link{nlsList.formula}}
}
\examples{
library(nlme)
data(CO2)
fm1 <- nlsList(SSasympOff, CO2)
fm1
}
\keyword{models}
