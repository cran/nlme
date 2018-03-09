/*
   Copyright 2018  The R Core Team

   This file is part of the nlme package for R and related languages
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

/*
  Replaces FORTRAN version in rs.f which did not handle special values
  such as NaN.
  
  return value in an arg for maximal portability -- see WRE section 6.6
*/

#include <R.h>
#include <math.h>

void F77_SUB(hypot)(double *a, double *b, double *p)
{
    *p = hypot(*a, *b);
}
