/*
   header file for the nlme package

   Copyright 1999-2001  Saikat DebRoy,
			Douglas Bates <bates@stat.wisc.edu>
   Copyright 2007-2013  The R Core Team

   This file is part of the nlme package for S and related languages
   and is made available under the terms of the GNU General Public
   License, version 2, or at your option, any later version,
   incorporated herein by reference.

   This program is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, a copy is available at
   http://www.r-project.org/Licenses/

*/

#ifndef NLME_BASE_H
#define NLME_BASE_H

#include <stdlib.h>
#include <string.h> // for memcpy
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#ifdef ENABLE_NLS
# include <libintl.h>
# define _(String) dgettext ("nlme", String)
#else
# define _(String) (String)
#endif

#define DNULLP (double *) 0

extern double sqrt_eps;
extern double xlower;

#endif /* NLME_BASE_H */
