/*
   Routines for fitting gnls models

   Copyright 1997-2005 Douglas M. Bates <bates@stat.wisc.edu>,
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

#include "nlOptimizer.h"
#include "matrix.h"
#include "nlmefit.h"

extern void corStruct_recalc(double *, int *, int *, double *);

/* gnls functions and variables     */

typedef struct gnls_struct {	/* Generalized nonlinear least squares structure */
  double *residuals, *gradient, *corFactor, *varWeights, minFactor,
    tolerance, *newtheta, *theta, *incr, *add_ons,
    new_objective, objective;
  double *result[1];
  int corOpt, varOpt, npar, ncol, N, nrdof, maxIter, *corDims;
  SEXP model;
  int conv_failure;
} *gnlsPtr;

static gnlsPtr
gnls_init(double *ptheta, int *dims, double *corFactor, double *varWeights,
	  int *corDims, double *settings, double *additional,
	  int corOpt, int varOpt, SEXP model)
{
  int nResult,
      npar = dims[0], // == p = pLen
      N = dims[1]; // == NReal == length(additional)
  gnlsPtr gnls = R_Calloc(1, struct gnls_struct);
  gnls->theta = ptheta;
  gnls->corFactor = corFactor;
  gnls->varWeights = varWeights;
  gnls->corDims = corDims;
  gnls->npar = npar;
  gnls->N = N;
  gnls->nrdof = N - npar;
  gnls->ncol = npar + 1;
  gnls->maxIter = (int) settings[0];
  gnls->minFactor = settings[1];
  gnls->tolerance = settings[2];
  gnls->newtheta = R_Calloc(npar, double);
  gnls->incr = R_Calloc(npar, double);
  gnls->varOpt = varOpt;
  gnls->corOpt = corOpt;
  gnls->add_ons = additional;
  gnls->model = model;
  gnls->result[0] = DNULLP;
  nResult = evaluate(ptheta, npar, model, gnls->result);
  gnls->result[0] = R_Calloc(nResult, double);
  return gnls;
}

static void
gnlsFree( gnlsPtr gnls )
{
  R_Free(gnls->newtheta);
  R_Free(gnls->incr);
  R_Free(gnls->result[0]);
  R_Free(gnls);
}

static double
gnls_objective(gnlsPtr gnls)
{
  int i, j;
  if(gnls->varOpt) {			/* variance function correction */
    for(i = 0; i < gnls->N; i++) {
      for(j = 0; j < gnls->ncol; j++) {
	*(gnls->result[0] + i + j * gnls->N) *= gnls->varWeights[i];
      }
    }
  }
  if(gnls->corOpt) {			/* correlation structure correction */
    corStruct_recalc(gnls->result[0], gnls->corDims, &gnls->ncol, gnls->corFactor);
  }
  gnls->residuals = gnls->result[0] + gnls->npar * gnls->N;
  gnls->gradient = gnls->result[0];
  return(d_sum_sqr(gnls->residuals, gnls->N));
}

static double
gnls_increment(gnlsPtr gnls)
{
  if (sqrt_eps == 0.0) sqrt_eps = sqrt(DBL_EPSILON);
  double* auxRes = R_Calloc(gnls->N, double);
  Memcpy(auxRes, gnls->residuals, gnls->N);
  QRptr aQR = QR(gnls->gradient, gnls->N, gnls->N, gnls->npar);
  QRsolve(aQR, gnls->residuals, gnls->N, 1L, gnls->incr, gnls->npar);
  QRqty(aQR, auxRes, gnls->N, 1L);
  double regSS = 0.;
  for(int i=0; i < gnls->npar; i++) {
    regSS += auxRes[i] * auxRes[i];
  }
  QRfree(aQR);
  R_Free(auxRes);
  return(sqrt(((double) gnls->nrdof) * regSS /
	      ((double) gnls->npar) * (gnls->new_objective - regSS)));
}

static int
gnls_iterate(gnlsPtr gnls)
{
  double factor, criterion;
  int iteration;
  SEXP model = gnls->model;

  Memcpy(gnls->newtheta, gnls->theta, gnls->npar);
  evaluate(gnls->theta, gnls->npar , model, gnls->result);
  gnls->new_objective = gnls->objective = gnls_objective(gnls);
  gnls->conv_failure = 0;
  for (factor = 1.0, iteration = 1; iteration <= gnls->maxIter;
       iteration++) {		/* outer iteration loop */
				/* increment and convergence criterion */
    criterion = gnls_increment(gnls);
    if (gnls->conv_failure) return(iteration); /* Unable to make increment */
    if (criterion < gnls->tolerance) return(iteration); /* successful completion */
    do {			/* inner loop for acceptable step size */
      if (factor < gnls->minFactor) {
	gnls->conv_failure = 1;
	return(iteration);
      }
      Memcpy(gnls->newtheta, gnls->theta, gnls->npar);
      d_axpy(gnls->newtheta, factor, gnls->incr, gnls->npar);
      evaluate(gnls->newtheta, gnls->npar , model, gnls->result);
      gnls->new_objective = gnls_objective(gnls);
      if (gnls->conv_failure) return(iteration); /* unable to evaluate objective */
      factor /= 2.0;
    } while (gnls->new_objective >= gnls->objective);
    factor *= 4.0;
    if (factor > 1.0)
      factor = 1.0;
    gnls->objective = gnls->new_objective;
    Memcpy(gnls->theta, gnls->newtheta, gnls->npar);
  }
  gnls->conv_failure = 2;	/* Maximum number of iterations exceeded */
  return(iteration - 1);
}

static void
gnls_wrapup(gnlsPtr gnls)
{
  SEXP model = gnls->model;
  evaluate(gnls->theta, gnls->npar , model, gnls->result);
  Memcpy(gnls->add_ons, gnls->result[0] + gnls->npar * gnls->N, gnls->N);
  gnls->objective = gnls_objective(gnls);
}

void
fit_gnls(double *ptheta, int *pdims, // = Dims
	 double *pcorFactor, double *pvarWeights, int *pcorDims,
	 double *settings, double *additional,
	 int *pcorOpt, int *pvarOpt, SEXP model)
{
  gnlsPtr gnls;

  PROTECT(model);
  if(sqrt_eps == 0.0) sqrt_eps = sqrt(DBL_EPSILON);
  gnls = gnls_init(ptheta, pdims, pcorFactor, pvarWeights, pcorDims,
		   settings, additional, *pcorOpt, *pvarOpt, model);
  settings[4] = (double) gnls_iterate(gnls);
  gnls_wrapup(gnls);
  settings[3] = gnls->conv_failure;
  settings[5] = gnls->objective;
  gnlsFree(gnls);
  UNPROTECT(1);
}
