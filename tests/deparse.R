# Tests to make sure that the deparse(as.vector(x)) construction is not
# tripping us up again.

library(nlme)
data(Loblolly)
fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
            data = Loblolly,
            fixed = Asym + R0 + lrc ~ 1,
            random = Asym ~ 1,
            start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
fm1

model <- height ~ SSasymp(age, Asym, R0, lrc)
fixed <- Asym + R0 + lrc ~ 1
random <- Asym ~ 1
start <- c(Asym = 103, R0 = -8.5, lrc = -3.3)
fm2 <- nlme(model,
            data = Loblolly,
            fixed = fixed,
            random = random,
            start = start)
fm2
