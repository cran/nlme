require(nls)

nonlinModel <- function( modelExpression, env,
                        paramNames = get(".parameters", envir = env)) {
  modelExpression <- modelExpression[[2]]
  thisEnv <- environment()
  offset <- 0
  ind <- vector("list", length(paramNames))
  names(ind) <- paramNames
  for( i in paramNames ) {
    ind[[ i ]] <- offset + seq( along = get(i, envir = env))
    offset <- offset + length( get(i, envir = env) )
  }
  modelValue <- eval(modelExpression, env)
  on.exit(remove(i, offset, paramNames))
  function( newPars) {
    if(!missing(newPars)) {
      for( i in names(ind) ) {
        assign( i, clearNames(newPars[ ind[[i]] ]), envir = env)
      }
      assign("modelValue", eval(modelExpression, env),
             envir = thisEnv)
    }
    modelValue
  }
}

