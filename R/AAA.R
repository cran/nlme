### $Id: AAA.R,v 1.5 2001/01/10 19:53:04 bates Exp $
###
### nlme for R
###
### Copyright 1999-2001 Douglas M. Bates <bates@stat.wisc.edu>,
###                     Saikat DebRoy <saikat@stat.wisc.edu>
###
### This file is part of the nlme library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

require( "nls" )

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
