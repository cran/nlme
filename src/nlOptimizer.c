/* $Id: nlOptimizer.c,v 1.1 1999/11/04 16:40:30 saikat Exp $

   Implementation of eval_model() and spread() for R.

   Copyright 1999 Saikat DebRoy <saikat@stat.wisc.edu>

   This file is part of the nlme library for S and related languages
   and is made available under the terms of the GNU General Public
   License, version 2, or at your option, any later version,
   incorporated herein by reference.

   This program is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public
   License along with this program; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
   MA 02111-1307, USA
 
*/

#include "nlOptimizer.h"

#ifndef R_S_H
#include "nonlin.h"
#endif /* R_S_H */

int
evaluate(double *param, longint nParam aMOD, double **value aSEV)
{
#ifdef R_S_H
  SEXP newPars, result;
  int i, nResult;

  PROTECT(newPars = allocVector(REALSXP, nParam));
  for(i = 0; i < nParam; i++)
    REAL(newPars)[i] = param[i];
  PROTECT(result = eval(lang2(model, newPars), R_GlobalEnv));
  nResult = LENGTH(result);
  if(value[0] == (double *) 0) {
    UNPROTECT(2);
    return(nResult);
  }
  for(i = 0; i < nResult; i++)
    value[0][i] = REAL(result)[i];
  UNPROTECT(2);
#else
  spread(param, nParam SEV);
  eval_model(TRUE SEV);
  value[0] = nl_results[0];
#endif
  return(-1);
}
