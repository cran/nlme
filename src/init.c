/*
   Copyright 2005-2018  The R Core Team

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

#include <stddef.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "nlmefit.h"
#include "nlOptimizer.h"
#include "pdMat.h"

extern void corStruct_factList(double *, int *, double *, double *);
extern void corStruct_recalc(double *, int *, int *, double *);
extern void symm_fullCorr(double *, int *, double *);
extern void symm_matList(double *, int *, int *, int *, double *);
extern void symm_factList(double *, int *, int *, int *p,
		  double *, double *);
extern void symm_recalc(double *, int *, int *, double *,
	    int *, int *, double *);
extern void nat_fullCorr(double *, int *, double *);
extern void nat_matList(double *, int *, int *, int *, double *);
extern void nat_factList(double *, int *, int *, int *,
	     double *, double *);
extern void nat_recalc(double *, int *, int *l, double *,
	   int *, int *, double *);
extern void AR1_matList(double *, int *, double *);
extern void AR1_factList(double *, int *, double *, double *);
extern void AR1_recalc(double *, int *, int *, double *, double *);
extern void CAR1_matList(double *, double *, int *, double *);
extern void CAR1_factList(double *, double *, int *,  double *, double *);
extern void CAR1_recalc(double *Xy, int *pdims, int *ZXcol,
	    double *par, double *time, double *logdet);
extern void ARMA_constCoef(int *, int *, double *);
extern void ARMA_unconstCoef(int *, int *, double *);
extern void ARMA_matList(double *, int *, int *, int *,
	     int *, int *, double *);
extern void ARMA_factList(double *, int *, int *, int *,
	      int *, int *, double *, double *);
extern void ARMA_recalc(double *, int *, int *, double *,
	    int *, int *, int *, int *, double *);
extern void compSymm_matList(double *, double *, int *, double *);
extern void compSymm_factList(double *, double *, int *,
		  double *, double *);
extern void compSymm_recalc(double *, int *, int *, double *,
		double *, double *);
extern void spatial_matList(double *, int *, double *, int *,
		double *, double *);
extern void spatial_factList(double *, int *, double *, int *,
		 double *, double *, double *);
extern void spatial_recalc(double *, int *, int *, double *,
	       double *, double *, int *, double *);

extern void fit_nlme(double *, double *, int *,
	 int *, int *, double *,
	 double *, int *, double *,
	 double *, int *, int *,
	 double, SEXP model);
extern void nlme_one_comp_first (int *, double *, double *);
extern void nlme_one_comp_open (int *, double *, double *);

extern void inner_perc_table(double *, int *, int *, int *,
		 int *, double *);

extern void natural_pd(double *, int *, double *);

R_CMethodDef CEntries[] = {
    {"corStruct_factList", (DL_FUNC) &corStruct_factList, 4},
    {"corStruct_recalc", (DL_FUNC) &corStruct_recalc, 4},
    {"symm_factList", (DL_FUNC) &symm_factList, 6},
    {"symm_matList", (DL_FUNC) &symm_matList, 5},
    {"symm_fullCorr", (DL_FUNC) &symm_fullCorr, 3},
    {"symm_recalc", (DL_FUNC) &symm_recalc, 7},
    {"nat_factList", (DL_FUNC) &nat_factList, 6},
    {"nat_matList", (DL_FUNC) &nat_matList, 5},
    {"nat_fullCorr", (DL_FUNC) &nat_fullCorr, 3},
    {"nat_recalc", (DL_FUNC) &nat_recalc, 7},
    {"AR1_factList", (DL_FUNC) &AR1_factList, 4},
    {"AR1_matList", (DL_FUNC) &AR1_matList, 3},
    {"AR1_recalc", (DL_FUNC) &AR1_recalc, 5},
    {"CAR1_factList", (DL_FUNC) &CAR1_factList, 5},
    {"CAR1_matList", (DL_FUNC) &CAR1_matList, 4},
    {"CAR1_recalc", (DL_FUNC) &CAR1_recalc, 6},
    {"ARMA_constCoef", (DL_FUNC) &ARMA_constCoef, 3},
    {"ARMA_unconstCoef", (DL_FUNC) &ARMA_unconstCoef, 3},
    {"ARMA_factList", (DL_FUNC) &ARMA_factList, 8},
    {"ARMA_matList", (DL_FUNC) &ARMA_matList, 7},
    {"ARMA_recalc", (DL_FUNC) &ARMA_recalc, 9},
    {"compSymm_factList", (DL_FUNC) &compSymm_factList, 5},
    {"compSymm_matList", (DL_FUNC) &compSymm_matList, 4},
    {"compSymm_recalc", (DL_FUNC) &compSymm_recalc, 6},
    {"spatial_factList", (DL_FUNC) &spatial_factList, 7},
    {"spatial_matList", (DL_FUNC) &spatial_matList, 6},
    {"spatial_recalc", (DL_FUNC) &spatial_recalc, 8},
    {"gls_loglik", (DL_FUNC) &gls_loglik, 5}, // 17-11-2015; Fixed sigma ...

    {"gls_estimate", (DL_FUNC) &gls_estimate, 8},
    {"fit_gnls", (DL_FUNC) &fit_gnls, 10},
    {"inner_perc_table", (DL_FUNC) &inner_perc_table, 6},
    {"mixed_loglik", (DL_FUNC) &mixed_loglik, 7}, // 17-11-2015; Fixed sigma ...
    {"mixed_decomp", (DL_FUNC) &mixed_decomp, 2},
    {"mixed_EM", (DL_FUNC) &mixed_EM, 10},        // 17-11-2015; Fixed sigma ...
    {"nlme_one_comp_first", (DL_FUNC) &nlme_one_comp_first, 3},
    {"nlme_one_comp_open", (DL_FUNC) &nlme_one_comp_open, 3},
    {"fit_nlme", (DL_FUNC) &fit_nlme, 14},        // 17-11-2015; Fixed sigma ...
    {"matrixLog_pd", (DL_FUNC) &matrixLog_pd, 3},
    {"logChol_pd", (DL_FUNC) &logChol_pd, 3},
    {"natural_pd", (DL_FUNC) &natural_pd, 3},
    {"compSymm_pd", (DL_FUNC) &compSymm_pd, 3},
    {"mixed_estimate", (DL_FUNC) &mixed_estimate, 8}, // 17-11-2015; Fixed sigma ...
    {"mixed_combined", (DL_FUNC) &mixed_combined, 11},
    {NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_nlme(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
