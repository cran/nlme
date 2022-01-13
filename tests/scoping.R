library("nlme")

## PR#18157
mygnls <- function (mydata)
    gnls(weight ~ SSlogis(Time, Asym, xmid, scal), data = mydata)
fm1 <- mygnls(Soybean) # failed in 3.1-153 with
## Error in stats::nls(formula = weight ~ SSlogis(Time, Asym, xmid, scal),  : 
##   object 'mydata' not found

## similarly, each of the following calls of
## nlme.formula(), nlsList.selfStart(), nlme.nlsList()
## using a self-starting model with local data would fail in 3.1-153 with
## Error in is.data.frame(data) : object 'mydata' not found
local({
    mydata <- subset(Loblolly, Seed < "307")
    fm2 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
                data = mydata, random = Asym ~ 1)
    fml <- nlsList(SSasymp, data = mydata)
    fm3 <- nlme(fml, random = Asym ~ 1)
})
