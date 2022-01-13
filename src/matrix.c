/*
   Basic matrix manipulations and QR decomposition

   Copyright 1997-2005  Douglas M. Bates <bates@stat.wisc.edu>,
			Jose C. Pinheiro,
			Saikat DebRoy
   Copyright 2007-2022  The R Core Team

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

#include <float.h> // for DBL_EPSILON

#include "matrix.h"
/* find qr decomposition, dqrdc2() is basis of R's qr(), also used by nlme */
void F77_NAME(dqrdc2)(double *x, int *ldx, int *n, int *p,
		      double *tol, int *rank,
		      double *qraux, int *pivot, double *work);
void F77_NAME(dqrls)(double *x, int *n, int *p, double *y, int *ny,
		     double *tol, double *b, double *rsd,
		     double *qty, int *k,
		     int *jpvt, double *qraux, double *work);

void
d_axpy(double *y, double a, double *x, int n)
{				/* y <- a * x + y  */
  while (n-- > 0) { *y++ += a * *x++; }
}

double
d_sum_sqr( double *x, int n )
{				/* sum(x * x) */
  double accum = 0.0;
  while (n-- > 0) { accum += *x * *x; x++; }
  return accum;
}

double
d_dot_prod( double *x, int incx, double *y, int incy, int n )
{				/* sum(x * y) */
  double accum = 0.0;
  while (n-- > 0) { accum += *x * *y; x +=incx; y += incy; }
  return accum;
}

double *
copy_mat(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{				/* y <- x */
  double * ret = y;
  while (ncol-- > 0) { Memcpy(y, x, nrow); y += ldy; x += ldx; }
  return ret;
}

double *
copy_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
 /* y <- t(x) */
{
  double * ret = y;
  int i, j;
  for (i = 0L; i < nrow; i++) {
    for (j = 0L; j < ncol; j++) { y[j] = x[i + j * ldx]; }
    y += ldy;
  }
  return ret;
}

double *
mult_mat(double *z, int ldz,
	 double *x, int ldx, int xrows, int xcols,
	 double *y, int ldy, int ycols)
{				/* z <- x %*% y */
  double *t, *tmp = R_Calloc((size_t)(xrows * ycols), double);
  int i, j;			/* use tmp so z can be either x or y */

  t = tmp;
  for (i = 0; i < ycols; i++) {
    for (j = 0; j < xcols; j++) {
      d_axpy(t, y[j], x + j * ldx, xrows);
    }
    t += xrows;
    y += ldy;
  }
  copy_mat(z, ldz, tmp, xrows, xrows, ycols);
  R_Free(tmp);
  return z;
}

static void
zero_mat(double *y, int ldy, int nrow, int ncol)
{				/* y[,] <- 0 */
  while (ncol-- > 0) {
    int i;
    for (i = 0; i < nrow; i++) { y[i] = 0.0; }
    y += ldy;
  }
}

QRptr
QR(double *mat, int ldmat, int nrow, int ncol)
{				/* Constructor for a QR object */
  QRptr value = R_Calloc((size_t) 1, struct QR_struct);
  int j;  double *work;

  if (sqrt_eps == 0.) { sqrt_eps = sqrt(DBL_EPSILON); }
  value->mat = mat;
  value->ldmat = ldmat;
  value->nrow = nrow;
  value->ncol = ncol;
  value->qraux = R_Calloc((size_t) ncol, double);
  value->pivot = R_Calloc((size_t) ncol, int);
  for (j = 0; j < ncol; j++) { (value->pivot)[j] = j; }
  work = R_Calloc( 2 * ncol, double );
  F77_CALL(dqrdc2) (mat, &ldmat, &nrow, &ncol, &sqrt_eps, &(value->rank),
		    value->qraux, value->pivot, work);
  R_Free(work);
  return value;
}

void
QRfree(QRptr this)
{				/* destructor for a QR object*/
  R_Free(this->pivot);
  R_Free(this->qraux);
  R_Free(this);
}

int
QRqty(QRptr this, double *ymat, int ldy, int ycol)
{				/* ymat <- qr.qty(this, ymat) */
  int j, info, task = 1000L;
  for (j = 0; j < ycol; j++) {
    double *col = ymat + j * ldy;
    F77_CALL(dqrsl) (this->mat, &(this->ldmat), &(this->nrow), &(this->ncol),
		     this->qraux, col, DNULLP, col, DNULLP, DNULLP, DNULLP,
		     &task, &info);
  }
  return info;
}

int
QRsolve( QRptr this, double *ymat, int ldy, int ycol, double *beta, int ldbeta )
{				/* beta <- qr.beta(this, ymat) */
  int j, info, task = 1100L;
  double *qty = R_Calloc( this->nrow, double ),
    *bb = R_Calloc( this->ncol, double );

  for (j = 0; j < ycol; j++) {
    Memcpy( qty, ymat, this->nrow );
    F77_CALL(dqrsl) (this->mat, &(this->ldmat), &(this->nrow), &(this->ncol),
		     this->qraux, qty, DNULLP, qty, bb, DNULLP,
		     DNULLP, &task, &info);
    Memcpy( beta, bb, this->ncol );
    ymat += ldy;
    beta += ldbeta;
  }
  R_Free( qty ); R_Free( bb );
  return info;
}

double
QRlogAbsDet(QRptr this)
{				/* log(abs(det(upper triangle))) */
  int j;
  double accum = 0.0;
  for (j = 0; j < this->rank; j++)
    accum += log(fabs(this->mat[j * (this->ldmat + 1L)]));
  return accum;
}

void
QRstoreR(QRptr this, double *dest, int ldDest)
{				/* store the R part into dest */
  int i;
  for (i = 0; i < this->ncol; i++) {
    Memcpy(dest + this->pivot[i] * ldDest, this->mat + i * this->ldmat,
	   ((i + 1) > this->rank) ? this->rank : i + 1);
  }
}

int
QR_and_rotate(double *mat, int ldmat, int nrow, int ncol,
	      double *DmHalf, int qi, int ndecomp,
	      double *logdet, double *store, int ldstr)
     /* Append DmHalf to the bottom of mat and take a QR decomposition
	of the first ndecomp columns.  Apply the rotations to the other
	columns.  Return the rank and increment log(abs(det(R11))). */
{
  int rank, arow = nrow + qi,  /* number of rows in augmented matrix */
    ndrow = ((arow > ndecomp) ? ndecomp : arow);
  double *aug = R_Calloc((size_t) arow * ncol, double);
  QRptr aQR;

  copy_mat(aug, arow, mat, ldmat, nrow, ncol);
  copy_mat(aug + nrow, arow, DmHalf, qi, qi, qi);
  aQR = QR(aug, arow, arow, ndecomp);
  if (logdet != DNULLP) { *logdet += QRlogAbsDet(aQR); }
  QRqty(aQR, aug + ndecomp * arow, arow, ncol - ndecomp);
  if (ldstr > 0) {
    QRstoreR(aQR, store, ldstr);
    copy_mat(store + ndecomp * ldstr, ldstr, aug + ndecomp * arow,
	     arow, ndrow, ncol - ndecomp);
  }
  if (qi < ndecomp) { zero_mat(mat, ldmat, nrow, ncol); }
  copy_mat(mat + ndecomp * ldmat, ldmat, aug + ndecomp * (arow + 1L),
	   arow, arow - ndrow, ncol - ndecomp);
  rank = aQR->rank;
  QRfree(aQR); R_Free(aug);
  return rank;
}
