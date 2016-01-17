library(nlme)
data(petrol, package = 'MASS')
Petrol <- petrol
Petrol[, 2:5] <- scale(Petrol[, 2:5], scale = FALSE)
pet3.lme <- lme(Y ~ SG + VP + V10 + EP,
                random = ~ 1 | No, data = Petrol, method="ML")
upet3 <- update(pet3.lme, Y ~ SG + VP + V10)
upet3
vc3 <- VarCorr(upet3)
upet2 <- lme(Y ~ SG + VP + V10, random = ~ 1 | No, data = Petrol, method = "ML")
stopifnot(
    all.equal(upet3, upet2, tol = 1e-15)
    ,
    all.equal(fixef(upet3),
	      c("(Intercept)" = 19.659375, SG = 0.125045632,
		VP = 2.27818601, V10 = 0.0672413592), tol = 1e-8)# 1e-9
    ,
    all.equal(as.numeric(vc3[,"StdDev"]),
	      c(0.00029397, 9.69657845), tol=1e-6)
)


