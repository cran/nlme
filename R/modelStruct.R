###         modelStruct - a virtual class of model structures
###
### Copyright 1997-2003  Jose C. Pinheiro <Jose.Pinheiro@pharma.novartis.com>,
###                      Douglas M. Bates <bates@stat.wisc.edu>

### Constructor
### There is no constructor function for this class (i.e. no function
### called modelStruct) because the class is virtual.
### Objects inheriting from this class are required to have a "conLin"
### (condensed linear model) attribute and a "pmap" (parameter map)
### attribute

###*# Methods for standard generics

coef.modelStruct <-
  function(object, unconstrained = TRUE, ...)
{
  unlist(lapply(object, coef, unconstrained))
}

"coef<-.modelStruct" <-
  function(object, ..., value)
{
  value <- as.numeric(value)
  parMap <- attr(object, "pmap")
  for(i in names(object)) {
    if (any(parMap[,i])) {
      coef(object[[i]]) <- value[parMap[,i]]
    }
  }
  object
}

formula.modelStruct <-
  function(x, ...)
{
  lapply(x, formula)
}

needUpdate.modelStruct <-
  function(object) any(unlist(lapply(object, needUpdate)))

print.modelStruct <-
  function(x, ...)
{
  for(i in names(x)) {
    if ((length(aux <- coef(x[[i]]))) > 0) {
      cat(paste(i, " parameters:\n"))
      print(aux)
    }
  }
}

print.summary.modelStruct <-
  function(x, ...)
{
  lapply(x, print, ...)
}

recalc.modelStruct <-
  function(object, conLin = attr(object, "conLin"), ...)
{
  for(i in rev(seq_along(object))) {
    conLin <- recalc(object[[i]], conLin)
    NULL
  }
  conLin
}

summary.modelStruct <-
  function(object, ...)
{
  value <- lapply(object, summary)
  class(value) <- "summary.modelStruct"
  value
}

## will not work as it is. fitted needs more than one argument!
update.modelStruct <-
  function(object, data, ...)
{
  if (needUpdate(object)) {
    object[] <- lapply(object, update, c(list("." = object), as.list(data)))
  }
  object
}

### Local Variables:
### mode:S
### End:

