/*
   Routines for fitting nlme models

   Copyright 2007-2022  The R Core Team
   Copyright 1997-2001  Douglas M. Bates <bates@stat.wisc.edu>,
			Jose C. Pinheiro,
			Saikat DebRoy

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

#define is_na_DOUBLE(x) ISNA(*(x))

/* nlme functions and variables     */

typedef struct nlme_struct {	/* Nonlinear mixed-effects structure */
    double *residuals, *gradient, *DmHalf, *corFactor, *varWeights,
	*newtheta, *theta, *incr, *add_ons, new_objective, objective, RSS,
	*sigma; // <- 17-11-2015; Fixed sigma patch; E van Willigen; Quant.Sol.
    int corOpt, varOpt, nparTot, ngrpTot, nrdof, *sgroups, *corDims,
	*npar, *pdClass, *pdims, *ZXoff, *ZXlen;
    double *result[1];
    dimPTR dd, d1;
    SEXP model;
    int conv_failure;
} *nlmePtr;

double sqrt_eps = 0.0;

static int *
make_sequential(int *dest, int *src, int n)
{
    /*  copy the pattern from src to dest */
    /*  but in sequential values starting */
    /*  from 0 */
    int val = 0, *ret = dest, sval;
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
nlme_init(double *ptheta, double *pDmHalf, int *pgroups, int *pdims,
	  int *pdClass, double *pcorFactor, double *pvarWeights,
	  int *pcorDims, double *additional, int *pcorOpt, int *pvarOpt,
	  // 17-11-2015; Fixed sigma patch; E van Willigen; Quantitative Solutions
	  double *sigma, SEXP model)
{
    int i, *src, nResult;
    nlmePtr nlme = R_Calloc(1, struct nlme_struct);
    nlme->pdims = pdims;
    nlme->DmHalf = pDmHalf;
    nlme->pdClass = pdClass;
    nlme->corFactor = pcorFactor;
    nlme->varWeights = pvarWeights;
    nlme->corDims = pcorDims;
    nlme->dd = dims(pdims);
    nlme->npar = R_Calloc(nlme->dd->Q + 1, int);
    nlme->sigma = sigma; // 17-11-2015; Fixed sigma patch; E van Willigen; Quantitative Solutions
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
    nlme->ngrpTot = 0;
    for (i = 0; i < (nlme->dd->Q + 2); i++) { nlme->ngrpTot += nlme->dd->ngrp[i]; }
    nlme->ZXoff = R_Calloc(nlme->ngrpTot, int);
    Memcpy(nlme->ZXoff, nlme->dd->ZXoff[0], nlme->ngrpTot);
    nlme->ZXlen = R_Calloc(nlme->ngrpTot, int);
    Memcpy(nlme->ZXlen, nlme->dd->ZXlen[0], nlme->ngrpTot);
    nlme->newtheta = R_Calloc(nlme->nparTot, double);
    nlme->incr =  R_Calloc(nlme->nparTot, double);
    nlme->model = model;
    nlme->result[0] = DNULLP;
    nResult = evaluate(ptheta, nlme->nparTot, model, nlme->result);
    nlme->result[0] = R_Calloc(nResult, double);
    return(nlme);
}

static void
nlmeFree(nlmePtr nlme)
{
    R_Free(nlme->newtheta);
    R_Free(nlme->incr);
    R_Free(nlme->npar);
    R_Free(nlme->ZXoff);
    R_Free(nlme->ZXlen);
    R_Free(nlme->result[0]);
    R_Free(nlme);
}

static void			/* undo changes in dd from internal_decomp */
restore_dims(nlmePtr nlme)
{
    nlme->dd->ZXrows = nlme->dd->N;
    Memcpy(nlme->dd->ZXoff[0], nlme->ZXoff, nlme->ngrpTot);
    Memcpy(nlme->dd->ZXlen[0], nlme->ZXlen, nlme->ngrpTot);
}

static void
nlme_wtCorrAdj(nlmePtr nlme)
{
    int i, j;
    if(nlme->varOpt) {		/* variance function adjustment */
	for(i = 0; i < nlme->dd->N; i++) {
	    for(j = 0; j < nlme->dd->ZXcols; j++) {
		*(nlme->result[0] + i + j * nlme->dd->N) *= nlme->varWeights[i];
	    }
	}
    }
    if(nlme->corOpt) {		/* correlation structure adjustment */
	corStruct_recalc(nlme->result[0], nlme->corDims, &nlme->dd->ZXcols,
			 nlme->corFactor);
    }
}

static double
nlme_RSS(nlmePtr nlme)
{
    nlme->residuals = nlme->result[0] + (nlme->dd->ZXcols - 1) * nlme->dd->N;
    nlme->gradient = nlme->result[0];
    nlme->RSS = d_sum_sqr(nlme->residuals, nlme->dd->N);
    return(nlme->RSS);
}

static double
nlme_objective(nlmePtr nlme)
{
    int i;
    double RSS, *srcB;

    RSS = nlme->RSS;
    for(i = 0, srcB = nlme->newtheta; i < nlme->dd->Q; i++) {
	double *work = R_Calloc(nlme->npar[i], double);
	mult_mat(work, (nlme->dd->ncol)[i], nlme->DmHalf + (nlme->dd->DmOff)[i],
		 (nlme->dd->ncol)[i], (nlme->dd->ncol)[i], (nlme->dd->ncol)[i],
		 srcB, (nlme->dd->ncol)[i], (nlme->dd->ngrp)[i]);
	RSS += d_sum_sqr(work, nlme->npar[i]);
	srcB += nlme->npar[i];
	R_Free(work);
    }
    return(RSS);
}

static void
nlme_workingRes(nlmePtr nlme)
{
    int i, j, k;
    double *theta = nlme->theta;


    for(j = 0; j < nlme->dd->Q; j++) {
	int nb = nlme->dd->ncol[j];
	double *res =
	    nlme->gradient + nlme->dd->ZXrows * (nlme->dd->ZXcols - 1);
	for(k = 0; k < nlme->dd->ngrp[j]; k++) {
	    double *Zjk = nlme->gradient + nlme->dd->ZXoff[j][k];

	    for(i = 0; i < nlme->dd->ZXlen[j][k]; i++) {
		*res += d_dot_prod(Zjk + i, nlme->dd->ZXrows, theta, 1, nb);
		res++;
	    }
	    theta += nb;
	}
    }
}


static double
nlme_increment(nlmePtr nlme)
{
    double predObj, *dest, *src, logLik, lRSS,
	*Ra = R_Calloc(nlme->dd->DmOff[nlme->dd->Q], double),
	*dc = R_Calloc(nlme->dd->Srows * nlme->dd->ZXcols, double)
/*      , *auxGrad = R_Calloc(nlme->dd->N * (nlme->dd->ZXcols - 1), double) */
	;
    double *incr = nlme->incr;
    double *theta = nlme->theta;
    int i, j, start, RML = 0;

    if (sqrt_eps == 0.0) sqrt_eps = sqrt(DBL_EPSILON);
/*    Memcpy(auxGrad, nlme->gradient, (nlme->dd->ZXcols - 1) * nlme->dd->N); */
    internal_decomp(nlme->dd, nlme->gradient);
    nlme_workingRes(nlme);
    internal_EM(nlme->dd, nlme->gradient, nlme->DmHalf, 20,
		// 17-11-2015; Fixed sigma patch; E van Willigen; Quantitative Solutions
		nlme->pdClass, &RML, &logLik, Ra, &lRSS, nlme->sigma);
    {
	statePTR st = R_Calloc(1, struct state_struct);
	int ntheta = count_DmHalf_pars( nlme->dd, nlme->pdClass ),
	    itrmcd, itncnt, msg, iagflg;
	double epsm,
	    *theta = R_Calloc(ntheta, double),
	    *typsiz = R_Calloc(ntheta, double),
	    *grad = R_Calloc(ntheta, double),
	    *newtheta = R_Calloc(ntheta, double),
	    *a = R_Calloc(ntheta * ntheta, double),
	    *work = R_Calloc(ntheta * 9, double);

	st->dd = nlme->dd;
	st->ZXy = nlme->gradient;
	st->pdClass = nlme->pdClass;
	st->RML = &RML;
	st->sigma = nlme->sigma; // 17-11-2015; Fixed sigma patch; E van Willigen; Quantitative Solutions

	generate_theta(theta, nlme->dd, nlme->pdClass, nlme->DmHalf);

	epsm = DBL_EPSILON;
	msg = 9;		/* don't inhibit checks but suppress output */
	for (i = 0; i < ntheta; i++) { typsiz[i] = 1.0; }
/*	iagflg = 1; */
/*	for (i = 0; i < nlme->dd->Q; i++) { */
/*	    if (nlme->pdClass[i] < 1 || nlme->pdClass[i] == 3 || nlme->pdClass[i] > 4) { */
/*		iagflg = 0; */
/*		break; */
/*	    } */
/*	} */
	iagflg = 0;		/* temporary modification */

	// Call the documented C API interface to R's nlm() :
	optif9(ntheta, ntheta, theta, (fcn_p) mixed_fcn, (fcn_p) mixed_grad,
	       /* no hessian: */ (d2fcn_p) 0,
	       st, typsiz, 1.0 /*fscale*/, 1 /*method*/, 1 /*iexp*/, &msg,
	       -1 /*ndigit*/, 20 /*itnlim*/, iagflg, 0 /*iahflg*/,
	       -1. /*dlt*/, pow(epsm, 1.0/3.0) /*gradtl*/, 0. /*stepmx*/,
	       sqrt(epsm) /*steptl*/, newtheta, &logLik, grad, &itrmcd, a,
	       work, &itncnt);
	if (msg == 0) {
	    generate_DmHalf(nlme->DmHalf, nlme->dd, nlme->pdClass, theta);
	}
	R_Free(work);
	R_Free(a);
	R_Free(newtheta);
	R_Free(grad);
	R_Free(typsiz);
	R_Free(theta);
	R_Free(st);
    }
    nlme->objective = nlme_objective(nlme);
    // 17-11-2015; Fixed sigma patch; E van Willigen; Quantitative Solutions
    internal_loglik(nlme->dd, nlme->result[0], nlme->DmHalf, &RML, dc, DNULLP, nlme->sigma);
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
    for(i = 0; i < (nlme->nparTot - nlme->npar[nlme->dd->Q]); i++) {
	incr[i] -= theta[i];
    }
    predObj = dc[nlme->dd->ZXcols * nlme->dd->ZXrows - 1];
    predObj = predObj * predObj;
    /*    regSS = nlme_RegSS(nlme, auxGrad); */	/* Regression Sum of Squares */
    R_Free(Ra); R_Free(dc);
/*    R_Free(auxGrad); */
    return(sqrt(((double) nlme->nrdof) * (nlme->objective - predObj) /
		(((double) nlme->nparTot) * predObj)));
}

static int
nlme_iterate(nlmePtr nlme, double *settings)
{
    double factor, criterion;
    SEXP model = nlme->model;
    double *newtheta = nlme->newtheta;
    double *theta = nlme->theta;
    int iteration;
    long maxIter = (long) settings[0];
    double minFactor = settings[1];
    double tolerance = settings[2];

    Memcpy(newtheta, theta, nlme->nparTot);
    evaluate(theta, nlme->nparTot , model, nlme->result);
    nlme_wtCorrAdj(nlme);
    nlme_RSS(nlme);
    settings[3] = nlme->conv_failure = 0;
    for (factor = 1.0, iteration = 1; iteration <= maxIter;
	 iteration++) {		/* outer iteration loop */
				/* increment and convergence criterion */
	criterion = nlme_increment(nlme);
	if (nlme->conv_failure) return(iteration); /* Unable to make increment  */
	if (criterion < tolerance) return(iteration); /* successful completion */
	int inner_it = 0;
	do {			/* inner loop for acceptable step size */
	    if (factor < minFactor) {
		settings[3] = 1;
		return(iteration);
	    }
	    Memcpy(newtheta, theta, nlme->nparTot);
	    d_axpy(newtheta, factor, nlme->incr, nlme->nparTot);
	    evaluate(newtheta, nlme->nparTot , model, nlme->result);
	    restore_dims(nlme);
	    nlme_wtCorrAdj(nlme);
	    nlme_RSS(nlme);
	    nlme->new_objective = nlme_objective(nlme);
	    if (nlme->conv_failure)
		return(iteration); /* unable to evaluate objective */
	    if((inner_it++) % 1000 == 0)
		R_CheckUserInterrupt();
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
nlme_wrapup(nlmePtr nlme)
{
    SEXP model = nlme->model;
    evaluate(nlme->theta, nlme->nparTot , model, nlme->result);
    Memcpy(nlme->add_ons, nlme->result[0], nlme->dd->N * nlme->dd->ZXcols);
    nlme->objective = nlme_objective(nlme);
    R_Free(nlme->npar);
    dimFree(nlme->dd);
}

void
fit_nlme(double *ptheta, double *pDmHalf, int *pgroups,
	 int *pdims, int *pdClass, double *pcorFactor,
	 double *pvarWeights, int *pcorDims, double *settings,
	 // 17-11-2015; Fixed sigma patch; E van Willigen; Quantitative Solutions
	 double *additional, int *pcorOpt, int *pvarOpt, double *sigma,
	 SEXP model)
{
    nlmePtr nlme;

    PROTECT(model);
    nlme = nlme_init(ptheta, pDmHalf, pgroups, pdims, pdClass,
		     pcorFactor, pvarWeights, pcorDims,
		     // 17-11-2015; Fixed sigma patch; E van Willigen; Quantitative Solutions
		     additional, pcorOpt, pvarOpt, sigma , model);
    if(sqrt_eps == 0.0) sqrt_eps = sqrt(DBL_EPSILON);
    settings[4] = (double) nlme_iterate(nlme, settings);
    nlme_wrapup(nlme);
    settings[5] = nlme->objective;
    nlmeFree(nlme);
    UNPROTECT(1);
}

void
nlme_one_comp_open (int *nrow, double *Resp, double *inmat)
{
    int i, nn = *nrow;
    double ke, ka, tl = 0, delta, C = 0, Ca = 0, interval,
	*Subject, *Time, *Conc, *Dose, *Interval, *V, *Ka, *Ke,
	sl = DBL_EPSILON;	/* sl is last subject number, usually */
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
nlme_one_comp_first (int *nrow, double *Resp, double *inmat)
{
    int nn = *nrow, mm = 0;
    double v, cl,
	*tl = R_Calloc(nn, double),
	*ds = R_Calloc(nn, double),
	sl = DBL_EPSILON;	/* sl is last subject number, usually */
				/* an integer but passed as double. */
				/* It is started at an unlikely value. */
    double
	*Subject = inmat,
	*Time    = inmat + nn,
	*Dose    = inmat + 2 * nn,
	*V       = inmat + 3 * nn,
	*Cl      = inmat + 4 * nn;
    for(int i = nn; i--; Resp++, Subject++, Time++, Dose++, V++, Cl++) {
	v = *V; cl = *Cl;
	*Resp = 0;
	if (*Subject != sl) {	/* new Subject */
	    if (is_na_DOUBLE(Dose)) {
		error(_("First observation on an individual must have a dose"));
	    }
	    sl = *Subject;
	    mm = 0;
	    tl[mm] = *Time;
	    ds[mm] = *Dose;
	} else {		/* same Subject */
	    if (!is_na_DOUBLE(Dose)) { /* Dose measurement */
		mm++;
		tl[mm] = *Time;
		ds[mm] = *Dose;
	    } else {		/* Concentration measurement */
		for(int j = 0; j <= mm; j++) {
		    *Resp += ds[j] * exp(-cl * (*Time - tl[j]) / v) / v;
		}
	    }
	}
    }
    R_Free(ds); R_Free(tl);
}
