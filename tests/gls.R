## reported by simon bond <shug0131@yahoo.co.uk> to R-help 2007-03-16

library(nlme)
x <- rnorm(10, 0.1, 1)
try(gls(x ~ 0))  # segfaulted in 3.1-79
