#define longint long int
			/* solve a positive definite system */
static longint
posDef_solve(double *mat, longint n, double *y, double *cvrg, longint *info)
{
 longint nn = n, zero = 0L, job;
 double *work;

 if (n == 1) {  
   if (*mat > 0.0) {
     *cvrg = *y * *y / *mat;
     *y = *y / *mat;
     return 0L;
   } else {
     return 1L;
   }
 }
 work = Calloc( (size_t) n, double );
 F77_CALL(chol) (mat, &nn, work, &zero, &zero, info);
 Free( work );
 if (*info < n) { return *info; }
 job = 11;
 F77_CALL(dtrsl) (mat, &nn, &nn, y, &job, info);
 *cvrg = d_dot_prod( y, 1L, y, 1L, nn );
 if (*info != 0) { return *info; }
 job = 1;
 F77_CALL(dtrsl) (mat, &nn, &nn, y, &job, info);
  
 return *info;
}

void
DmHalf( double *Delta, longint *pdims, double *pars, longint *pdClass )
{				/* test the generate_DmHalf function */
  dimPTR dd = dims( pdims );
  generate_DmHalf( Delta, dd, pdClass, pars);
  dimFree( dd );
}

void				/* general pd, logChol parametrization */
logChol_pd(double *L, longint *q, double *l)
{
  longint i, qq = *q;
  double *ll = l + qq;
  L[0] = exp(*l);
  for(i = 1; i < qq; i++) {
    L[i * (qq + 1)] = exp(l[i]);
    Memcpy(L + i * qq, ll, i);
    ll += i;
  }
}

void 
spher_pd(double *L, longint *q, double *l)
{
  longint i, j, qq = *q;
  double *src = l + qq, aux, aux1;
  for(i = 0; i < qq; i++) {
    aux = exp(l[i]);
    for(j = 0; j < i; j++, src++) {
      aux1 = PI * exp(*src)/(1 + exp(*src));
      L[i * qq + j] = aux * cos(aux1);
      aux *= sin(aux1);
    }
    L[i * (qq + 1)] = aux;
  }
}

void 
Givens_pd(double *L, longint *q, double *l)
{
  longint i, j, m, nn, n = *q,  k = (n * (n + 1) / 2) - 1;
  double  *dest, *src, ang, aux0, aux1, aux2;

  if (n > 1) {
    for(nn = n, aux0 = 0, dest = L + n * n - 1, src = l; nn--; 
	dest -= (n + 1), src++)
      {
	aux0 += exp(*src);
	*dest = sqrt(aux0);		/* L <- diag(sqrt(eigenvalues)) */
      }
    /* Transposed Givens rotations */
    for(i = (n - 2); i >= 0; i--) {
      for(j = (n - 1); j > i; j--, k--) {
	ang = PI * exp(l[k]) / (exp(l[k]) + 1); /* converting to */
                                                /* (0, pi)       */
	aux0 = cos(ang);
	aux1 = sin(ang);
	for(m = 0; m < n; m++) {
	  aux2 = L[i * n + m];
	  L[i * n + m] = L[i * n + m] * aux0 - L[j * n + m] *  aux1;  
	  L[j * n + m] = aux2 * aux1 + L[j * n + m] * aux0;
	}
      }
    }
  }
  else *L = exp(*l);
}
