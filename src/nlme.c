/* $Id: nlme.c,v 1.3 2000/07/03 18:22:49 bates Exp $

   Routines for fitting nlme models

   Copyright 1997, 1999 Douglas M. Bates <bates@stat.wisc.edu>,
                        Jose C. Pinheiro <jcp@research.bell-labs.com>,
			Saikat DebRoy <saikat@stat.wisc.edu>

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

#include "nlOptimizer.h"
#include "matrix.h"
#include "nlmefit.h"

extern void corStruct_recalc(double *, longint *, longint *, double *);

#ifdef R_S_H
#define is_na_DOUBLE(x) ISNA(*(x))
#endif /* R_S_H */

/* nlme functions and variables     */

typedef struct nlme_struct {	/* Nonlinear mixed-effects structure */
  double *residuals, *gradient, *DmHalf, *corFactor, *varWeights,
    *newtheta, *theta, *incr, *add_ons, new_objective, objective;
  longint corOpt, varOpt, nparTot, nrdof, *sgroups, *corDims, *npar;
  double *result[1];
  dimPTR dd;
#ifdef R_S_H
  SEXP model;
#endif /* R_S_H */  
  int conv_failure;
} *nlmePtr;

double sqrt_eps = 0.0;

static longint *
make_sequential(longint *dest, longint *src, longint n)
{                    
  /*  copy the pattern from src to dest */
  /*  but in sequential values starting */
  /*  from 0 */
  longint val = 0, *ret = dest, sval;
  if (n <= 0) return dest;
  sval = *src++; *dest++ = val;
  while (--n) {
    if (*src != sval) {sval = *src; val++;}
    src++;
    *dest++ = val;
  }
  return ret;
}

static nlmePtr
nlme_init(double *ptheta, double *pDmHalf, longint *pgroups, 
	  longint *pdims, double *pcorFactor, double *pvarWeights, 
	  longint *pcorDims, double *additional, longint *pcorOpt,
	  longint *pvarOpt aMOD)
{
  longint i, *src, nResult;
  nlmePtr nlme = Calloc(1, struct nlme_struct);
  nlme->DmHalf = pDmHalf;
  nlme->corFactor = pcorFactor;
  nlme->varWeights = pvarWeights;
  nlme->corDims = pcorDims;
  nlme->dd = dims(pdims);
  nlme->npar = Calloc(nlme->dd->Q + 1, longint);
  for(i = 0, nlme->nparTot = 0; i <= nlme->dd->Q; i++) {
    nlme->npar[i] = (nlme->dd->ncol)[i] * (nlme->dd->ngrp)[i];
    nlme->nparTot += nlme->npar[i];
  }
  nlme->nrdof = nlme->dd->N - nlme->nparTot;
  nlme->sgroups = pgroups;
  for(i = 0, src = nlme->sgroups; i < nlme->dd->Q; i++) {
    make_sequential(src, src, nlme->dd->N);
    src += nlme->dd->N;
  }
  nlme->corOpt = *pcorOpt;
  nlme->varOpt = *pvarOpt;
  nlme->theta = ptheta; 
  nlme->add_ons = additional;
  nlme->newtheta = Calloc(nlme->nparTot, double);
  nlme->incr =  Calloc(nlme->nparTot, double);
#ifdef R_S_H
  nlme->model = model;
  nlme->result[0] = DNULLP;
  nResult = evaluate(ptheta, nlme->nparTot MOD, nlme->result SEV);
  nlme->result[0] = Calloc(nResult, double);
#endif /* R_S_H */
  return(nlme);
}

void
nlmeFree(nlmePtr nlme)
{
  Free(nlme->newtheta);
  Free(nlme->incr);
#ifdef R_S_H
  Free(nlme->result[0]);
#endif
  Free(nlme);
}

static double
nlme_objective(nlmePtr nlme)
{
  double RSS, *srcB;
  longint i, j;
  if(nlme->varOpt) {		/* variance function correction */
    for(i = 0; i < nlme->dd->N; i++) {
      for(j = 0; j < nlme->dd->ZXcols; j++) {
	*(nlme->result[0] + i + j * nlme->dd->N) *= nlme->varWeights[i];
      }
    }
  }
  if(nlme->corOpt) {		/* correlation structure correction */
    corStruct_recalc(nlme->result[0], nlme->corDims, &nlme->dd->ZXcols,
		     nlme->corFactor);
  }
  nlme->residuals = nlme->result[0] + (nlme->dd->ZXcols - 1) * nlme->dd->N;
  nlme->gradient = nlme->result[0];
  RSS = d_sum_sqr(nlme->residuals, nlme->dd->N);
  for(i = 0, srcB = nlme->newtheta; i < nlme->dd->Q; i++) { 
    double *work = Calloc(nlme->npar[i], double);
    mult_mat(work, (nlme->dd->ncol)[i], nlme->DmHalf + (nlme->dd->DmOff)[i],
	     (nlme->dd->ncol)[i], (nlme->dd->ncol)[i], (nlme->dd->ncol)[i],
	     srcB, (nlme->dd->ncol)[i], (nlme->dd->ngrp)[i]);
    RSS += d_sum_sqr(work, nlme->npar[i]);
    srcB += nlme->npar[i];
    Free(work);
  }
  return(RSS);
}

static void
nlme_correction(nlmePtr nlme)
{
  longint i, j, nq;
  double *theta = nlme->theta;

  for(i = 0; i < nlme->dd->N; i++) {
    for(j = 0, nq = 0; j < nlme->dd->Q; j++) {
      nlme->residuals[i] +=
	d_dot_prod(nlme->gradient+(nlme->dd->ZXoff)[j][0] + i, 
		   nlme->dd->N, 
		   theta+nq+(nlme->sgroups+j*nlme->dd->N)[i]*nlme->dd->ncol[j],
		   1L, (nlme->dd->ncol)[j]);
      nq += nlme->npar[j];
    }
    nlme->residuals[i] += 
      d_dot_prod(nlme->gradient + (nlme->dd->ZXoff)[nlme->dd->Q][0] + i,
		 nlme->dd->N, theta + nq, 1L, (nlme->dd->ncol)[nlme->dd->Q]);
  }
}

static double
nlme_RegSS(nlmePtr nlme, double *mat)
{
  longint i, j, nq;
  double regSS = 0, aux = 0, *src;
  double *incr = nlme->incr;

  for(i = 0; i < nlme->dd->N; i++) {
    for(j = 0, nq = 0, aux = 0; j < nlme->dd->Q; j++) {
      aux += d_dot_prod(mat + (nlme->dd->ZXoff)[j][0] + i, 
			nlme->dd->N,
			incr+nq+(nlme->sgroups+j*nlme->dd->N)[i]*nlme->dd->ncol[j], 
			1L, nlme->dd->ncol[j]);
      nq += nlme->npar[j];
    }
    aux += d_dot_prod(mat + (nlme->dd->ZXoff)[nlme->dd->Q][0] + i,
		      nlme->dd->N, incr + nq, 1L, (nlme->dd->ncol)[nlme->dd->Q]);
    regSS += aux * aux;
  }
  for(i = 0, src = incr; i < nlme->dd->Q; i++) { 
    double *work = Calloc((nlme->npar)[i], double);
    mult_mat(work, (nlme->dd->ncol)[i], nlme->DmHalf + (nlme->dd->DmOff)[i],
	     (nlme->dd->ncol)[i], (nlme->dd->ncol)[i], (nlme->dd->ncol)[i],
	     src, (nlme->dd->ncol)[i], (nlme->dd->ngrp)[i]);
    regSS += d_sum_sqr(work, (nlme->npar)[i]);
    src += (nlme->npar)[i];
    Free(work);
  }
  return(regSS);
}

static double
nlme_increment(nlmePtr nlme)
{
  double regSS, *dc = Calloc(nlme->dd->Srows * nlme->dd->ZXcols, double),
    ll, *dest, *src, 
    *auxGrad = Calloc(nlme->dd->N * (nlme->dd->ZXcols - 1), double);
  double *incr = nlme->incr;
  double *theta = nlme->theta;

  longint i, j, start, RML = 0;
  if (!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  Memcpy(auxGrad, nlme->gradient, (nlme->dd->ZXcols - 1) * nlme->dd->N);
  nlme_correction(nlme);
  ll = internal_loglik(nlme->dd, nlme->result[0], nlme->DmHalf, &RML, dc, DNULLP);
  internal_estimate(nlme->dd, dc);
  src = dc +  (nlme->dd->ZXcols - 1) * nlme->dd->Srows;
  dest = incr;
  for(i = 0, start = 0; i <= nlme->dd->Q; i++) {
    for(j = 0; j < (nlme->dd->ngrp)[i]; j++) {
      Memcpy(dest, src + ((nlme->dd->SToff)[i][j] - start), (nlme->dd->ncol)[i]);
      dest += (nlme->dd->ncol)[i];
    }
    start += (nlme->dd->ncol)[i] * nlme->dd->Srows;
  }
  for(i = 0; i < nlme->nparTot; i++) {
    incr[i] -= theta[i];
  }
  regSS = nlme_RegSS(nlme, auxGrad);	/* Regression Sum of Squares */
  Free(dc); Free(auxGrad);
  return(sqrt(((double) nlme->nrdof) * regSS / 
	      (((double) nlme->nparTot) * (nlme->new_objective - regSS))));
}

static longint
nlme_iterate(nlmePtr nlme, double *settings aSEV)
{
  double factor, criterion;
#ifdef R_S_H
  SEXP model = nlme->model;
#endif /* R_S_H */
  double *newtheta = nlme->newtheta;
  double *theta = nlme->theta;
  longint iteration;
  long maxIter = (long) settings[0];
  double minFactor = settings[1];
  double tolerance = settings[2];
  int conv_failure;
  S_EVALUATOR

  Memcpy(newtheta, theta, nlme->nparTot);
  evaluate(theta, nlme->nparTot MOD, nlme->result SEV);
  nlme->new_objective = nlme->objective = nlme_objective(nlme);
  settings[3] = conv_failure = 0;
  for (factor = 1.0, iteration = 1; iteration <= maxIter;
       iteration++) {		/* outer iteration loop */
				/* increment and convergence criterion */
    criterion = nlme_increment(nlme); 
    if (conv_failure) return(iteration); /* Unable to make increment  */
    if (criterion < tolerance) return(iteration); /* successful completion */
    do {			/* inner loop for acceptable step size */
      if (factor < minFactor) {
	settings[3] = 1;
	return(iteration);
      }
      Memcpy(newtheta, theta, nlme->nparTot);
      d_axpy(newtheta, factor, nlme->incr, nlme->nparTot);
      evaluate(newtheta, nlme->nparTot MOD, nlme->result SEV);
      nlme->new_objective = nlme_objective(nlme);
      if (conv_failure) return(iteration); /* unable to evaluate objective */
      factor /= 2.0;
    } while (nlme->new_objective >= nlme->objective);
    factor *= 4.0;
    if (factor > 1.0)
      factor = 1.0;
    nlme->objective = nlme->new_objective;
    Memcpy(theta, newtheta, nlme->nparTot);
  }
  settings[3] = 2;  /* Maximum number of iterations exceeded */
  return(iteration - 1);
}

static void
nlme_wrapup(nlmePtr nlme aSEV)
{
#ifdef R_S_H
  SEXP model = nlme->model;
#endif
  S_EVALUATOR
  evaluate(nlme->theta, nlme->nparTot MOD, nlme->result SEV);
  Memcpy(nlme->add_ons, nlme->result[0], nlme->dd->N * nlme->dd->ZXcols);
  nlme->objective = nlme_objective(nlme);
  Free(nlme->npar);
  dimFree(nlme->dd);
}

void
fit_nlme(double *ptheta, double *pDmHalf, longint *pgroups, 
	longint *pdims, double *pcorFactor, double *pvarWeights, 
	longint *pcorDims, double *settings, double *additional, 
	longint *pcorOpt, longint *pvarOpt aMOD)
{
  nlmePtr nlme;
  S_EVALUATOR
#ifdef R_S_H
  PROTECT(model);
#endif /* R_S_H */
  nlme = nlme_init(ptheta, pDmHalf, pgroups, pdims,
    pcorFactor, pvarWeights, pcorDims,
    additional, pcorOpt, pvarOpt MOD);
  if(!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  settings[4] = (double) nlme_iterate(nlme, settings SEV);
  nlme_wrapup(nlme SEV);
  settings[5] = nlme->objective;
  nlmeFree(nlme);
#ifdef R_S_H
  UNPROTECT(1);
#endif /* R_S_H */
}

void 
nlme_one_comp_open (longint *nrow, double *Resp, double *inmat)
{
  longint i, nn = *nrow;
  double ke, ka, tl = 0, delta, C = 0, Ca = 0, interval,
    *Subject, *Time, *Conc, *Dose, *Interval, *V, *Ka, *Ke,
         sl = DOUBLE_EPS;	/* sl is last subject number, usually */
				/* an integer but passed as double. */
				/* It is started at an unlikely value. */
  Subject = inmat;
  Time = inmat + nn;
  Conc = inmat + 2 * nn;
  Dose = inmat + 3 * nn;
  Interval = inmat + 4 * nn;
  V = inmat + 5 * nn;
  Ka = inmat + 6 * nn;
  Ke = inmat + 7 * nn;
  for(i = nn; i--; Resp++, Subject++, Time++, Conc++, Dose++,
                   Interval++, V++, Ka++, Ke++) {
    ke = *Ke; ka = *Ka;
    if (*Subject != sl) {	/* new Subject */
      sl = *Subject;
      tl = *Time;
      *Resp = 0;
      if (!is_na_DOUBLE(Interval)) {		/* steady-state dosing */
	interval = *Interval;
	C = *Dose * ka * (1/(1 - exp(-ke * interval)) -
			  1/(1 - exp(-ka * interval)))/
			    (*V * (ka - ke));
	Ca = *Dose / (*V * (1 - exp(-ka * interval)));
      } else {			/* non-steady-state */
	C = 0;
	Ca = *Dose/ *V;
      }
    } else {			/* same Subject */
      if (!is_na_DOUBLE(Dose)) {
	if (!is_na_DOUBLE(Interval)) {	/* steady-state dosing */
	  interval = *Interval;
	  C = *Dose * ka * (1/(1 - exp(-ke * interval)) -
			    1/(1 - exp(-ka * interval)))/
			      (*V * (ka - ke));
	  Ca = *Dose / (*V * (1 - exp(-ka * interval)));
	} else {		/* non-steady-state */
	  delta = *Time - tl;
	  C = C*exp(-ke * delta) +
	    Ca*ka*(exp(-ke*delta) - exp(-ka*delta))/(ka -ke);
	  Ca = Ca * exp(-ka*delta) + *Dose / *V;
	}
	tl = *Time;
	*Resp = 0;
      } else if (!is_na_DOUBLE(Conc)) {
	delta = *Time - tl;
	*Resp = C * exp(-ke * delta) + Ca * ka *
	  (exp(-ke * delta) - exp(-ka * delta))/(ka - ke);
      } else *Resp = 0;
    }
  }
}

/* Phenobarbital Model */ 

void 
nlme_one_comp_first (longint *nrow, double *Resp, double *inmat)
{
  longint i, j, nn = *nrow, mm = 0;
  double v, cl, *tl = Calloc(nn, double), *ds = Calloc(nn, double),
    *Subject, *Time, *Dose, *V, *Cl,
         sl = DOUBLE_EPS;	/* sl is last subject number, usually */
				/* an integer but passed as double. */
				/* It is started at an unlikely value. */
  Subject = inmat;
  Time = inmat + nn;
  Dose = inmat + 2 * nn;
  V = inmat + 3 * nn;
  Cl = inmat + 4 * nn;
  for(i = nn; i--; Resp++, Subject++, Time++, Dose++, V++, Cl++) {
    v = *V; cl = *Cl;
    *Resp = 0;
    if (*Subject != sl) {	/* new Subject */
      if (is_na_DOUBLE(Dose)) {
	PROBLEM
	  "First observation on an individual must have a dose"
	  RECOVER(NULL_ENTRY);
      }
      sl = *Subject;
      mm = 0;
      tl[mm] = *Time;
      ds[mm] = *Dose;
    } else {			/* same Subject */
      if (!is_na_DOUBLE(Dose)) { /* Dose measurement */
	mm++;
	tl[mm] = *Time;
	ds[mm] = *Dose;
      } else {			/* Concentration measurement */
	for(j = 0; j <= mm; j++) {
	  *Resp += ds[j] * exp(-cl * (*Time - tl[j]) / v) / v;
	}
      }
    }
  }
  Free(ds); Free(tl);
}
