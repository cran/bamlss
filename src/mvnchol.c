#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

/*
Function: log_dmvnchol

Input:
Y observations (matrix)
PAR parameters (matrix)
N sample size (integer)
K dimension of multivariate distribution

Output: Vector of log-densities
*/
SEXP log_dmvncholC(SEXP Y, SEXP PAR, SEXP N, SEXP K)
{
/* Pointers for obs and par */
  double *Yptr = REAL(Y);
  double *PARptr = REAL(PAR);

/* integers from input */
  int n = asInteger(N);
  int k = asInteger(K);

/* loop counters */
  int i, j, l, pos;

/* integer vector for mapping a lower triangular matrix  */
  int lowtri[k*k];

/* residuals for single obs */
  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);

/* inverse of L matrix times residuals for single obs */
  SEXP Linvymu;
  PROTECT(Linvymu = allocVector(REALSXP, k));
  double *Linvymuptr = REAL(Linvymu);

/* return value: vector of log-densities */
  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

/* initialize three terms of log-density */
  double lpi = -0.5 * (float)k * log(2.0 * 3.14159265358979323846);
  double det = 0.0;
  double nrm = 0.0;

/* map lower triangular */
  i = 0;
  for(j = 0; j < k; j++) {
    for(l = 0; l < k; l++) {
      lowtri[l + j * k] = 0;
      if(l > j) {
	i += 1;
	lowtri[l + j * k] = i;
      }
/*      printf(" i = %d, j = %d, l = %d, lowtri = %d \n", i, j, l, lowtri[l + j * k]); */
    }
  }

/* HERE comes the computation of the log-density */
  for(i = 0; i < n; i++) {
/* compute term2 the determinante */
    det = 0.0;
    for(j = 0; j < k; j++) {
      det += log(PARptr[i + n * (j + k)]);
    }

/* compute term3 the norm */
    /* y - mu */
    for(j = 0; j < k; j++) {
      ymuptr[j] = Yptr[i + j * n] - PARptr[i + j * n];
    }

    /* inverse of L times vector of residuals */
    for(j = 0; j < k; j++) {
      Linvymuptr[j] = PARptr[i + n * (j + k)] * ymuptr[j];
      for(l = 0; l < j; l++) {
	pos = k + k + lowtri[j + k * l] - 1;
        Linvymuptr[j] += PARptr[i + n * pos] * ymuptr[l];
      }
    }

    /* L2-norm */
    nrm = 0.0;
    for(j = 0; j < k; j++) {
      nrm += pow(Linvymuptr[j], 2.0);
    }

/* sum up terms of log likelihood */
    dptr[i] = lpi + det - 0.5 * nrm;
  }

/* return */
  UNPROTECT(3);
  return d;
}



/*
Function: mu_score_mvncholC

Input:
Y observations (matrix)
PAR parameters (matrix)
N sample size (integer)
K dimension of multivariate distribution
II index of mu vector for partial derivative

Output: Vector of scores for (II)th element of mu vector
*/
SEXP mu_score_mvncholC(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP II)
{
/* Pointers for obs and par */
  double *Yptr = REAL(Y);
  double *PARptr = REAL(PAR);

/* integers from input */
  int n = asInteger(N);
  int k = asInteger(K);
  int ii = asInteger(II) - 1;

/* loop counters */
  int i, j, l, pos;

/* integer vector for mapping a lower triangular matrix  */
  int lowtri[k*k];
  
  double temp;

/* residuals for single obs */
  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);

/* return value: vector of mu-scores */
  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

/* inverse of L matrix times residuals for single obs */
  SEXP Linvymu;
  PROTECT(Linvymu = allocVector(REALSXP, k));
  double *Linvymuptr = REAL(Linvymu);



/* map lower triangular */
  i = 0;
  for(j = 0; j < k; j++) {
    for(l = 0; l < k; l++) {
      lowtri[l + j * k] = 0;
      if(l > j) {
	i += 1;
	lowtri[l + j * k] = i;
      }
/*      printf(" i = %d, j = %d, l = %d, lowtri = %d \n", i, j, l, lowtri[l + j * k]); */
    }
  }


  for(i = 0; i < n; i++) {
    /* y - mu */
    for(j = 0; j < k; j++) {
      ymuptr[j] = Yptr[i + j * n] - PARptr[i + j * n];
    }

    /* inverse of L times vector of residuals */
    for(j = 0; j < k; j++) {
      Linvymuptr[j] = PARptr[i + n * (j + k)] * ymuptr[j];
      for(l = 0; l < j; l++) {
	pos = k + k + lowtri[j + k * l] - 1;
        Linvymuptr[j] += PARptr[i + n * pos] * ymuptr[l];
      }
    }

    /* now calculate entry ii of L inv transpose times residual vector*/
    temp = Linvymuptr[ii] * PARptr[i + n * (ii + k)];
    for(j = ii + 1; j < k; j++) {
      pos = k + k + lowtri[j + k * ii] - 1;
      temp += PARptr[i + n * pos] * Linvymuptr[j];
    }

    /* assign to vector  */
    dptr[i] = temp;
  }


  /* return */
  UNPROTECT(3);
  return d;
}



/*
Function: lamdiag_score_mvncholC

Input:
Y observations (matrix)
PAR parameters (matrix)
N sample size (integer)
K dimension of multivariate distribution
II index of lamdiag for partial derivative

Output: Vector of scores for (II)th diagonal of L_inv
*/
SEXP lamdiag_score_mvncholC(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP II)
{
/* Pointers for obs and par */
  double *Yptr = REAL(Y);
  double *PARptr = REAL(PAR);

/* integers from input */
  int n = asInteger(N);
  int k = asInteger(K);
  int ii = asInteger(II) - 1;

/* loop counters */
  int i, j, l, pos;

/* integer vector for mapping a lower triangular matrix  */
  int lowtri[k*k];
  
  double temp;

/* residuals for single obs */
  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);

/* return value: vector of lamdiag-scores */
  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);


/* map lower triangular */
  i = 0;
  for(j = 0; j < k; j++) {
    for(l = 0; l < k; l++) {
      lowtri[l + j * k] = 0;
      if(l > j) {
	i += 1;
	lowtri[l + j * k] = i;
      }
/*      printf(" i = %d, j = %d, l = %d, lowtri = %d \n", i, j, l, lowtri[l + j * k]); */
    }
  }

  for(i = 0; i < n; i++) {
    /* y - mu */
    for(j = 0; j < k; j++) {
      ymuptr[j] = Yptr[i + j * n] - PARptr[i + j * n];
    }
   
    temp = PARptr[i + n * (ii + k)] * ymuptr[ii]; 
    for(j = 0; j < ii; j++) {
      pos = k + k + lowtri[ii + k * j] - 1;
      temp += PARptr[i + n * pos] * ymuptr[j];
    }

    /* assign to vector  */
    dptr[i] = 1 - PARptr[i + n * (ii + k)] * ymuptr[ii] * temp; 

  }


  /* return */
  UNPROTECT(2);
  return d;
}


/*
Function: lambda_score_mvncholC

Input:
Y observations (matrix)
PAR parameters (matrix)
N sample size (integer)
K dimension of multivariate distribution
II row index of lambda for partial derivative
JJ column index of lambda for partial derivative

Output: Vector of scores for (II)th diagonal of L_inv
*/
SEXP lambda_score_mvncholC(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP II, SEXP JJ)
{
/* Pointers for obs and par */
  double *Yptr = REAL(Y);
  double *PARptr = REAL(PAR);

/* integers from input */
  int n = asInteger(N);
  int k = asInteger(K);
  int ii = asInteger(II) - 1;
  int jj = asInteger(JJ) - 1;

/* loop counters */
  int i, j, l, pos;

/* integer vector for mapping a lower triangular matrix  */
  int lowtri[k*k];
  
  double temp;

/* residuals for single obs */
  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);

/* return value: vector of lamdiag-scores */
  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);


/* map lower triangular */
  i = 0;
  for(j = 0; j < k; j++) {
    for(l = 0; l < k; l++) {
      lowtri[l + j * k] = 0;
      if(l > j) {
	i += 1;
	lowtri[l + j * k] = i;
      }
/*      printf(" i = %d, j = %d, l = %d, lowtri = %d \n", i, j, l, lowtri[l + j * k]); */
    }
  }

  for(i = 0; i < n; i++) {
    /* y - mu */
    for(j = 0; j < k; j++) {
      ymuptr[j] = Yptr[i + j * n] - PARptr[i + j * n];
    }
   
    temp = PARptr[i + n * (jj + k)] * ymuptr[jj]; 
    for(j = 0; j < jj; j++) {
      pos = k + k + lowtri[jj + k * j] - 1;
      temp += PARptr[i + n * pos] * ymuptr[j];
    }

    /* assign to vector  */
    dptr[i] = 0 - ymuptr[ii] * temp; 

  }


  /* return */
  UNPROTECT(2);
  return d;
}


