/* $Id: nlOptimizer.h,v 1.3 2000/07/03 18:22:49 bates Exp $

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

#ifndef NLME_NLOPTIMIZER_H
#define NLME_NL_OPTIMIZER_H

#include "base.h"

#ifdef R_S_H
#include "Rinternals.h"
#endif

extern int evaluate(double *, longint aMOD, double ** aSEV);
extern void fit_gnls(double *, longint *, double *, double *, longint
		     *, double *, double *, longint *, longint * aMOD);
#endif /* NLME_NLOPTIMIZER_H */
