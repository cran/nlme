/* $Id: nlmefit.h,v 1.1 1999/11/04 16:40:30 saikat Exp $

   header file for the nlme package

   Copyright 1999 Saikat DebRoy <saikat@stat.wisc.edu>

   This file is part of the nlme library for R and related languages
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

#ifndef NLME_NLMEFIT_H
#define NLME_NLMEFIT_H
#include "base.h"

typedef struct dim_struct {
  longint
    N,				/* number of observations in original data */
    ZXrows,			/* number of rows in ZXy  */
    ZXcols,			/* number of columns in ZXy */
    Q,				/* number of levels of random effects */
    Srows,			/* number of rows in decomposition */
    *q,				/* dimensions of the random effects */
    *ngrp,			/* numbers of groups at each level */
    *DmOff,			/* offsets into the DmHalf array */
    *ncol,			/* no. of columns decomposed at each level */
    *nrot,			/* no. of columns rotated at each level */
    **ZXoff,			/* offsets into ZXy */
    **ZXlen,			/* groups lengths */
    **SToff,			/* offsets into storage */
    **DecOff,			/* offsets into decomposition */
    **DecLen;			/* decomposition group lengths */
} *dimPTR;

extern dimPTR dims(longint *);
extern void dimFree(dimPTR);
extern double internal_loglik(dimPTR, double *, double *, longint *,
			      double *, double *);
#ifndef R_S_H
extern void mixed_combined(double *, longint *, double *, longint *,
			   longint *, longint *, double *, double *,
			   double *, longint *);
#endif
extern void mixed_decomp(double *, longint *);
extern void mixed_loglik(double *, longint *, double *, longint *,
			 double *, double *);
extern void internal_estimate(dimPTR, double *);
extern void mixed_estimate(double *, longint *, double *, longint *,
			   double *, double *, longint *);
extern void mixed_EM(double *, longint *, double *, longint *,
		     longint *, longint *, double *, double *, double *);
extern void mixed_calcf(longint *, double *, longint *, double *,
			longint *, double *, void (*)(void));
extern void mixed_calcgh(longint *, double *, longint *, double *,
			 double *, longint *, double *, void (*)(void));
extern void gls_loglik(double *, longint *, double *, double *);
extern void gls_estimate(double *, longint *, double *, double *,
			 double *, double *, longint *, longint *);

#endif /* NLME_NLMEFIT_H */
