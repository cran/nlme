library(nlme)
data(Assay)
as1 <- lme(logDens~sample*dilut, data=Assay,
           random=pdBlocked(list(
                     pdIdent(~1),
                     pdIdent(~sample-1),
                     pdIdent(~dilut-1))))

as1s <- update(as1, random=pdCompSymm(~sample-1))
(an.1s <- anova(as1, as1s)) # non significant
stopifnot(
    all.equal(drop(data.matrix(an.1s[2,-1])),
              c(Model = 2, df = 33, AIC = -10.958851, BIC = 35.280663,
                logLik = 38.479425, Test = 2,
                L.Ratio = 0.11370211, `p-value` = 0.73596807), tol=8e-8))

as1S <- update(as1, . ~ sample+dilut) # dropping FE interaction
tools::assertWarning(anova(as1, as1S))# REML not ok for different FE.
as1M  <- update(as1,  method = "ML")
as1SM <- update(as1S, method = "ML")
(anM <- anova(as1M, as1SM)) # anova() OK: comparing MLE fits
## ==> significant: P ~= 0.0054
stopifnot(
    all.equal(drop(data.matrix(anM[2,])[,-(1:2)]),
	      c(df = 14, AIC = -169.588248, BIC = -140.267424, logLik = 98.7941241,
		Test = 2, L.Ratio = 39.7345188, `p-value` = 0.0053958561),
	      tol = 8e-8)
)

