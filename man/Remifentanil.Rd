% File nlme/man/Remifentanil.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{Remifentanil}
\title{Pharmacokinetics of Remifentanil}
\alias{Remifentanil}
\docType{data}
\description{
  Intravenous infusion of remifentanil (a strong analgesic) in different
  rates over varying time periods was applied to a total of 65 patients.
  Concentration measurements of remifentanil were taken along with several
  covariates resulting in the \code{Remifentanil} data frame with 2107 rows
  and 12 columns.
}
\usage{
data("Remifentanil", package = "nlme")
}
\format{
  This data frame (of class \code{"\link{groupedData}"}, specifically
  \code{"nfnGroupedData"}) contains the following columns:
  \describe{
    \item{\code{ID}:}{numerical (patient) IDs.}
    \item{\code{Subject}:}{an \code{\link{ordered}} factor with 65 levels (of the
      \code{ID}s): \code{30} < \code{21} < \code{25} < \code{23} <
      \code{29} < \dots \dots < \code{36} < \code{6} < \code{5} < \code{10} < \code{9}.}
    \item{\code{Time}:}{time from beginning of infusion in minutes (\code{\link{numeric}}).}
    \item{\code{conc}:}{remifentanil concentration in [ng / ml] (numeric).}
    \item{\code{Rate}:}{infusion rate in [µg / min].}
    \item{\code{Amt}:}{amount of remifentanil given in the current time interval in [µg].}
    \item{\code{Age}:}{age of the patient in years.}
    \item{\code{Sex}:}{gender of the patient, a \code{\link{factor}} with levels
      \code{Female} and \code{Male}.}
    \item{\code{Ht}:}{height of the patient in cm.}
    \item{\code{Wt}:}{weight of the patient in kg.}
    \item{\code{BSA}:}{body surface area (DuBois and DuBois 1916): \eqn{%
	BSA := Wt^{0.425} \cdot Ht^{0.725} \cdot 0.007184}{%
	BSA := Wt^0.425      *  Ht^0.725      *  0.007184}.}
    \item{\code{LBM}:}{lean body mass (James 1976), with slightly different
      formula for
        men \eqn{LBM_m := 1.1  Wt - 128 (Wt/Ht)^2}, and
      women \eqn{LBM_f := 1.07 Wt - 148 (Wt/Ht)^2}.}
  }
}
\author{of this help page: Niels Hagenbuch and Martin Maechler, SfS ETH Zurich.}
\source{
  Pinheiro, J. C. and Bates, D. M. (2000).
  \emph{Mixed-Effects Models in S and S-PLUS}, Springer, New York.
}
\references{
  Minto CF, Schnider TW, Egan TD, Youngs E, Lemmens HJM, Gambus PL, Billard
  V, Hoke JF, Moore KHP, Hermann DJ, Muir KT, Mandema JW, Shafer SL (1997).
  Influence of age and gender on the pharmacokinetics and pharmacodynamics
  of remifentanil: I. Model development.
  \emph{Anesthesiology} \bold{86} 1, 10--23.
  \doi{10.1097/00000542-199701000-00004}

  Charles F. Minto, Thomas W. Schnider and Steven L. Shafer (1997).
  Pharmacokinetics and Pharmacodynamics of Remifentanil: II. Model Application.
  \emph{Anesthesiology} \bold{86} 1, 24--33.
  \doi{10.1097/00000542-199701000-00005}
}
\examples{
plot(Remifentanil, type = "l", lwd = 2) # shows the 65 patients' remi profiles

## The same on  log-log  scale  (*more* sensible for modeling ?):
plot(Remifentanil, type = "l", lwd = 2, scales = list(log=TRUE))

str(Remifentanil)
summary(Remifentanil)

plot(xtabs(~Subject, Remifentanil))
summary(unclass(table(Remifentanil$Subject)))
## between 20 and 54 measurements per patient (median: 24; mean: 32.42)

## Only first measurement of each patient :
dim(Remi.1 <- Remifentanil[!duplicated(Remifentanil[,"ID"]),]) #  65 x 12

LBMfn <- function(Wt, Ht, Sex) ifelse(Sex == "Female",
                                        1.07 * Wt - 148*(Wt/Ht)^2,
                                        1.1  * Wt - 128*(Wt/Ht)^2)
with(Remi.1,
     stopifnot(all.equal(BSA, Wt^{0.425} * Ht^{0.725} * 0.007184,
                         tolerance = 1.5e-5),
               all.equal(LBM, LBMfn(Wt, Ht, Sex),
                         tolerance = 7e-7)
))

## Rate: typically  3 µg / kg body weight, but :
sunflowerplot(Rate ~ Wt, Remifentanil)
abline(0,3, lty=2, col=adjustcolor("black", 0.5))
}
\keyword{datasets}
