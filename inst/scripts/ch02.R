#-*- R -*-

library( nlme )
library( SASmixed )   # needed for the PBIB data
options( width = 65, digits = 5 )
options( contrasts = c(unordered = "contr.helmert", ordered = "contr.poly") )

# Chapter 2    Theory and Computational Methods for Linear Mixed-Effects Models

# 2.2   Likelihood Estimation for LME Models

Xmat <- matrix( c(1, 1, 1, 1, 8, 10, 12, 14), ncol = 2 )
Xmat
Xqr <- qr( Xmat )               # creates a QR structure
qr.R( Xqr )                     # returns R
qr.Q( Xqr )                     # returns Q-truncated
qr.Q( Xqr, complete = TRUE )    # returns the full Q

data( Rail )
fm1Rail.lme <- lme( travel ~ 1, data = Rail, random = ~ 1 | Rail,
       control = list( msVerbose = TRUE ) )
fm1Rail.lme <- lme( travel ~ 1, data = Rail, random = ~ 1 | Rail,
   control = list( msVerbose = TRUE, niterEM = 0 ))

data( Machines )
fm1Machine <-
  lme( score ~ Machine, data = Machines, random = ~ 1 | Worker )
fm2Machine <- update( fm1Machine, random = ~ 1 | Worker/Machine )
anova( fm1Machine, fm2Machine )

data( Orthodont )
OrthoFem <- Orthodont[ Orthodont$Sex == "Female", ]
fm1OrthF <- lme( distance ~ age, data = OrthoFem,
    random = ~ 1 | Subject )
fm2OrthF <- update( fm1OrthF, random = ~ age | Subject )
orthLRTsim <- simulate.lme( fm1OrthF, fm2OrthF, nsim = 1000 )
#plot( orthLRTsim, df = c(1, 2) )    # produces Figure 2.3

machineLRTsim <- simulate.lme(fm1Machine, fm2Machine, nsim= 1000)
#plot( machineLRTsim, df = c(0, 1),      # produces Figure 2.4
# layout = c(4,1), between = list(x = c(0, 0.5)) )

data( ergoStool )
stoolLRTsim <-
  simulate.lme( m1 = list(fixed = effort ~ 1, data = ergoStool,
                          random = ~ 1 | Subject),
                m2 = list(fixed = effort ~ Type),
                method = "ML", nsim = 1000 )
#plot( stoolLRTsim, df = c(3, 4) )    # Figure 2.5
data( PBIB )
pbibLRTsim <-
    simulate.lme( m1 = list( fixed = response ~ 1, data = PBIB,
                    random = ~ 1 | Block ),
                 m2 = list( fixed = response ~ Treatment ),
                 method = "ML", nsim = 1000 )
#plot( pbibLRTsim, df = c(14,16,18), weights = FALSE )    # Figure 2.6

summary( fm2Machine )

data( PBIB )
fm1PBIB <- lme( response ~ Treatment, data = PBIB, random = ~ 1 )
anova( fm1PBIB )
fm2PBIB <- update( fm1PBIB, method = "ML" )
fm3PBIB <- update( fm2PBIB, response ~ 1 )
anova( fm2PBIB, fm3PBIB )
anova( fm2Machine )

# cleanup

proc.time()
q()
