/*
   header file for the nlme package

   Copyright 1999-2001  Saikat DebRoy
   Copyright 2007-2016  The R Core Team

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

#ifndef NLME_MATRIX_H
#define NLME_MATRIX_H
#include "base.h"
void F77_NAME(chol)(double *a, int *lda, int *n, double *v, int *info);
int F77_NAME(rs)(int *nm, int *n, double *a, double *w,
		 int *matz, double *z, double *fv1, double *fv2, int *ierr);
#include <R_ext/Linpack.h>

typedef struct QR_struct {
  double *mat, *qraux;
  int *pivot, rank, ldmat, nrow, ncol;
} *QRptr;

extern void d_axpy(double *, double, double *, int);
extern double d_dot_prod(double *, int, double *, int, int);
extern double d_sum_sqr( double *, int);
extern double *copy_mat(double *, int, double *, int, int, int);
extern double *copy_trans(double *, int, double *, int, int, int);
extern double *mult_mat(double *, int, double *, int, int, int, double *,
			int, int);
extern QRptr QR(double *, int, int, int);
extern void QRfree(QRptr);
extern int QRqty(QRptr, double *, int, int);
extern int QRsolve(QRptr, double *, int, int, double *, int);
extern double QRlogAbsDet(QRptr);
extern void QRstoreR(QRptr, double *, int);
extern int QR_and_rotate(double *, int, int, int,
			 double *, int, int, double *,
			 double *, int);
#endif /* NLME_MATRIX_H */
