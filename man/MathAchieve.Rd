\name{MathAchieve}
\alias{MathAchieve}
\non_function{}
\title{Mathematics achievement scores}
\description{
    The \code{MathAchieve} data frame has 7185 rows and 6 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{School}{
	    an ordered factor identifying the school that the student attends
	}
	\item{Minority}{
	    a factor with levels
	    \code{No} 
	    \code{Yes}
	    indicating if the student is a member of a minority racial group.
	}
	\item{Sex}{
	    a factor with levels
	    \code{Male} 
	    \code{Female} 
	}
	\item{SES}{
	    a numeric vector of socio-economic status.
	}
	\item{MathAch}{
	    a numeric vector of mathematics achievement scores.
	}
	\item{MEANSES}{
	    a numeric vector of the mean SES for the school.
	}
    }
}
\details{
  Each row in this data frame contains the data for one student.
}
%\source{}
\examples{
data(MathAchieve)
summary(MathAchieve)
}
\keyword{datasets}
