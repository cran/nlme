\name{Muscle}
\alias{Muscle}
\non_function{}
\title{Contraction of heart muscle sections}
\description{
    The \code{Muscle} data frame has 60 rows and 3 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Strip}{
	    an ordered factor indicating the strip of muscle being measured.
	}
	\item{conc}{
	    a numeric vector giving the concentration of CaCl2
	}
	\item{length}{
	    a numeric vector giving the shortening ofthe heart muscle strip.
	}
    }
}
\details{
    Baumann and Waldvogel (1963) describe data on the shortening of heart
    muscle strips dipped in a CaCl\eqn{_2}{2} solution.  The muscle strips are
    taken from the left auricle of a rat's heart.}
\source{
    Baumann, F. and Waldvogel, F. (1963), La restitution pastsystolique
    de la contraction de l'oreillette gauche du rat. Effets de divers
    ions et de l'acetylcholine, \emph{Helvetica Physiologica Acta},
    \bold{21}.
}
\examples{
data(Muscle)
}
\keyword{datasets}
