#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Complex.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


/* (1) Helper functions. */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt, names;
  PROTECT(elmt = R_NilValue);
  PROTECT(names = getAttrib(list, R_NamesSymbol));

  for(int i = 0; i < length(list); i++) {
	  if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	  }
  }

  UNPROTECT(2);

  return elmt;
}


int getListElement_index(SEXP list, const char *str)
{
  SEXP names;
  PROTECT(names = getAttrib(list, R_NamesSymbol));

  int j = -1;

  for(int i = 0; i < length(list); i++) {
	  if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    j = i;
	    break;
	  }
  }

  UNPROTECT(1);

  return j;
}


void pvec(SEXP vec)
{
  int i, n = length(vec);
  double *vecptr = REAL(vec);
  for(i = 0; i < n; i++, vecptr++)
    Rprintf(" %g", *vecptr);
  Rprintf("\n");
}
    
void pmat(SEXP mat)
{
  int i,j,n = nrows(mat), k = ncols(mat);
  Rprintf("   ");
  for(j = 0; j < k; ++j)
    Rprintf("[%d] ", j);
  Rprintf("\n");
  for(i = 0; i < n; ++i) {
    Rprintf("[%d]", i);
    for(j = 0; j < k; ++j)
      Rprintf(" %g", REAL(mat)[i + j * n]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void pmat_int(SEXP mat)
{
  int i,j,n = nrows(mat), k = ncols(mat);
  Rprintf("   ");
  for(j = 0; j < k; ++j)
    Rprintf("[%d] ", j);
  Rprintf("\n");
  for(i = 0; i < n; ++i) {
    Rprintf("[%d]", i);
    for(j = 0; j < k; ++j)
      Rprintf(" %d", INTEGER(mat)[i + j * n]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void merr()
{
  char *m = "stopped";
  error(m);
}


/* (2) Main IWLS functions. */
SEXP iwls_eval(SEXP fun, SEXP response, SEXP eta, SEXP id, SEXP rho)
{
  SEXP R_fcall, rval;

  PROTECT(R_fcall = lang4(fun, response, eta, id));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


/* Map linear predictor to parameter scale */
SEXP map2par(SEXP fun, SEXP eta, SEXP rho)
{
  SEXP R_fcall, rval;

  PROTECT(R_fcall = lang2(fun, eta));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


SEXP eval_prior(SEXP fun, SEXP theta, SEXP rho)
{
  SEXP R_fcall, rval;

  PROTECT(R_fcall = lang2(fun, theta));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


SEXP get_S_mat(SEXP fun, SEXP theta, SEXP hyper, SEXP rho)
{
  SEXP R_fcall, rval;

  if(isNull(hyper))
    PROTECT(R_fcall = lang2(fun, theta));
  else
    PROTECT(R_fcall = lang3(fun, theta, hyper));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


SEXP eval_dmvnorm_log(SEXP fun, SEXP x, SEXP mu, SEXP sigma, SEXP rho)
{
  SEXP R_fcall, rval;

  PROTECT(R_fcall = lang4(fun, x, mu, sigma));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


/* logPost evaluation */
double lp_eval(SEXP fun, SEXP g, SEXP x,
  SEXP family, SEXP response, SEXP eta,
  SEXP id, SEXP rho)
{
  SEXP R_fcall, t, rval;

  PROTECT(t = R_fcall = allocList(7));
  SET_TYPEOF(R_fcall, LANGSXP);

  SETCAR(R_fcall, fun);
  t = CDR(t); SETCAR(t, g);
  t = CDR(t); SETCAR(t, x);
  t = CDR(t); SETCAR(t, family);
  t = CDR(t); SETCAR(t, response);
  t = CDR(t); SETCAR(t, eta);
  t = CDR(t); SETCAR(t, id);

  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);
  return REAL(rval)[0];
}


/* Univariate slice sampling */
SEXP uni_slice(SEXP g, SEXP x, SEXP family, SEXP response, SEXP eta, SEXP id, SEXP j,
  SEXP W, SEXP M, SEXP LOWER, SEXP UPPER, SEXP logPost, SEXP rho)
{
  int nProtected = 0;
  int jj = INTEGER(j)[0] - 1;

  int m = INTEGER(M)[0] + 1;
  double w = REAL(W)[0];
  double lower = REAL(LOWER)[0];
  double upper = REAL(UPPER)[0];

  SEXP gL, gR;
  PROTECT(gL = duplicate(g));
  ++nProtected;
  PROTECT(gR = duplicate(g));
  ++nProtected;

  double *gLptr = REAL(gL);
  double *gRptr = REAL(gR);
  double *gptr = REAL(g);

  double x0 = gptr[jj];
  double gx0 = lp_eval(logPost, g, x, family, response, eta, id, rho);

  GetRNGstate();
  double logy = gx0 - rexp(1); 
  double u = runif(0.0, w);
  PutRNGstate();

  gLptr[jj] = x0 - u;
  gRptr[jj] = x0 + (w - u);

  if(m > 1) {
    GetRNGstate();
    int J = floor(runif(0.0, m));
    PutRNGstate();
    int K = (m - 1) - J;
    while(J > 0) {
      if(gLptr[jj] <= lower)
        break;
      if(lp_eval(logPost, gL, x, family, response, eta, id, rho) <= logy)
        break;
      gLptr[jj] = gLptr[jj] - w;
      J = J - 1;
    }
    while(K > 0) {
      if(gRptr[jj] >= upper)
        break;
      if(lp_eval(logPost, gR, x, family, response, eta, id, rho) <= logy)
        break;
      gRptr[jj] = gRptr[jj] + w;
      K = K - 1;
    }
  }

  if(gLptr[jj] < lower) {
    gLptr[jj] = lower;
  }
  if(gRptr[jj] > upper) {
    gRptr[jj] = upper;
  }

  int run = 1;
  
  while(run > 0) {
    gptr[jj] = runif(gLptr[jj], gRptr[jj]);

    double gx1 = lp_eval(logPost, g, x, family, response, eta, id, rho);

    if(gx1 >= logy)
      run = 0;

    if(gptr[jj] > x0) {
      gRptr[jj] = gptr[jj];
    } else {
      gLptr[jj] = gptr[jj];
    }
  }

  UNPROTECT(nProtected);
  return g;
}


/* Compute the centroid of a polygon. */
SEXP cpos(SEXP p, SEXP K)
{
  int i, n, k;
  n = INTEGER(K)[0];
  k = n + 1;
    
  double tmp, *pptr, asum, xsum, ysum;
    
  pptr = REAL(p);
    
  asum = 0.0;
  xsum = 0.0;
  ysum = 0.0;
    
  for(i = 0; i < n; i++) {
    tmp = pptr[i] * pptr[i + k + 1] - pptr[i + 1] * pptr[i + k];
    asum += tmp;
    xsum += (pptr[i] + pptr[i + 1]) * tmp;
    ysum += (pptr[i + k] + pptr[i + k + 1]) * tmp;
  }
        
  tmp = 1 / (3 * asum);

  SEXP pos;
  PROTECT(pos = allocVector(REALSXP, 2));

  REAL(pos)[0] = tmp * xsum;
  REAL(pos)[1] = tmp * ysum;

  UNPROTECT(1);
    
  return pos;
}


/* Compute reduced weights and residuals. */
void xbin_fun(SEXP ind, SEXP weights, SEXP e, SEXP xweights, SEXP xrres, SEXP order)
{
  int i;
  int j = 0;
  int n = length(ind);
  int k = 0;

  double *weightsptr = REAL(weights);
  double *eptr = REAL(e);
  double *xweightsptr = REAL(xweights);
  double *xrresptr = REAL(xrres);
  int *indptr = INTEGER(ind);
  int *orderptr = INTEGER(order);

  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }
    k = orderptr[i] - 1;

    if(ISNA(weightsptr[k]))
      weightsptr[k] = 1.490116e-08;
    if(weightsptr[k] < 1e-10)
      weightsptr[k] = 1e-10;
    if(weightsptr[k] < 0.0)
      weightsptr[k] = -1.0 * weightsptr[k];
    if(weightsptr[k] > 1e+10)
      weightsptr[k] = 1e+10;

    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }
}


/* Process derivatives. */
SEXP process_derivs(SEXP x, SEXP w)
{
  int i;
  int n = length(x);

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, n));

  double *xptr = REAL(x);
  double *rvalptr = REAL(rval);

  for(i = 0; i < n; i++) {
    rvalptr[i] = xptr[i];
    if(ISNA(xptr[i]))
      rvalptr[i] = 1.490116e-08;
    if(xptr[i] > 1e+10)
      rvalptr[i] = 1e+10;
    if(LOGICAL(w)[0]) {
      if(xptr[i] < 0.0)
        rvalptr[i] = -1.0 * xptr[i];
    } else {
      if(xptr[i] < -1e+10)
        rvalptr[i] = -1e+10;
    }
  }

  UNPROTECT(1);

  return rval;
}


/* Fast computation of sum of diagonal. */
SEXP sum_diag(SEXP x, SEXP N)
{
  int i;
  int n = INTEGER(N)[0];

  double *xptr = REAL(x);
  double sum = 0.0;

  for(i = 0; i < n; i++) {
    sum += xptr[i + n * i];
  }

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, 1));
  REAL(rval)[0] = sum;
  UNPROTECT(1);

  return rval;
}


SEXP sum_diag2(SEXP x, SEXP y)
{
  int c, d, k;
  int n = ncols(x);

  double *xptr = REAL(x);
  double *yptr = REAL(y);
  double sum1 = 0.0;
  double sum2 = 0.0;

  for(c = 0; c < n; c++) {
    for(d = c; d < n; d++) {
      for(k = 0; k < n; k++) {
        sum1 = sum1 + xptr[c + k * n] * yptr[k + d * n];
      }
      if(c == d)
        sum2 += sum1;
      sum1 = 0.0;
    }
  }

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, 1));
  REAL(rval)[0] = sum2;
  UNPROTECT(1);

  return rval;
}


/* Efficient IWLS sampling. */
SEXP gmcmc_iwls(SEXP family, SEXP theta, SEXP id,
  SEXP eta, SEXP response, SEXP x, SEXP z, SEXP e, SEXP id2, SEXP W, SEXP rho)
{
  int i, j, k, nProtected = 0;
  int n = INTEGER(getListElement(x, "nobs"))[0];
  int fixed = LOGICAL(getListElement(x, "fixed"))[0];
  int fxsp = LOGICAL(getListElement(x, "fxsp"))[0];
  int nW = length(W);
  if(nW > 1) {
    if(nW != n)
      nW = 1;
  }
  double *Wptr = REAL(W);

  SEXP theta2;
  PROTECT(theta2 = duplicate(getListElement(getListElement(theta,
    CHAR(STRING_ELT(id, 0))), CHAR(STRING_ELT(id, 1)))));
  double *thetaptr = REAL(theta2);
  ++nProtected;

  int S_ind = getListElement_index(x, "S");
  int ntau2;
  if(fixed > 0) {
    ntau2 = 0;
  } else {
    ntau2 = length(VECTOR_ELT(x, S_ind));
  }
  int nc = length(theta2);
  if(fixed < 1) {
    nc -= ntau2;
  }

  SEXP gamma0, gamma1, tau2;
  PROTECT(gamma0 = allocVector(REALSXP, nc));
  double *gamma0ptr = REAL(gamma0);
  ++nProtected;

  PROTECT(gamma1 = allocVector(REALSXP, nc));
  double *gamma1ptr = REAL(gamma1);
  ++nProtected;

  PROTECT(tau2 = allocVector(REALSXP, ntau2));
  double *tau2ptr = REAL(tau2);
  ++nProtected;

  /* Evaluate loglik, weights and score vector. */
  SEXP peta;
  PROTECT(peta = map2par(getListElement(family, "map2par"), eta, rho));
  ++nProtected;
  int ll_ind = getListElement_index(family, "loglik");
  double pibeta = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, id2, rho))[0];

  SEXP weights;
  PROTECT(weights = iwls_eval(getListElement(getListElement(family, "hess"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *weightsptr = REAL(weights);
  ++nProtected;

  SEXP score;
  PROTECT(score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *scoreptr = REAL(score);
  ++nProtected;

  SEXP eta2;
  PROTECT(eta2 = duplicate(eta));
  double *etaptr = REAL(getListElement(eta2, CHAR(STRING_ELT(id, 0))));
  ++nProtected;

  /* Create weighted matrix */
  int X_ind = getListElement_index(x, "X");
  int nr = nrows(VECTOR_ELT(x, X_ind));

  /* More pointers needed. */
  double *zptr = REAL(z);
  double *eptr = REAL(e);
  double *xweightsptr = REAL(getListElement(x, "weights"));
  double *xrresptr = REAL(getListElement(x, "rres"));
  double *XWptr = REAL(getListElement(x, "XW"));
  double *XWXptr = REAL(getListElement(x, "XWX"));
  double *Xptr = REAL(VECTOR_ELT(x, X_ind));
  double *Sptr;
  int *idptr = INTEGER(getListElement(getListElement(x, "binning"), "match.index"));
  int *indptr = INTEGER(getListElement(getListElement(x, "binning"), "sorted.index"));
  int *orderptr = INTEGER(getListElement(getListElement(x, "binning"), "order"));

  /* Handling fitted.values. */
  SEXP fitname;
  PROTECT(fitname = allocVector(STRSXP, 1));
  ++nProtected;
  SET_STRING_ELT(fitname, 0, mkChar("fitted.values"));
  double *fitptr = REAL(getAttrib(theta2, fitname));
  double *fitrptr = REAL(getListElement(x, "fit.reduced"));

  /* Start. */
  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  j = 0;
  int jj;
  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      for(jj = 0; jj < nc; jj++) {
        XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
      }
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }
    k = orderptr[i] - 1;

    if(ISNA(weightsptr[k]))
      weightsptr[k] = 1.490116e-08;
    if(weightsptr[k] < 1e-10)
      weightsptr[k] = 1e-10;
    if(weightsptr[k] < 0.0)
      weightsptr[k] = -1.0 * weightsptr[k];
    if(weightsptr[k] > 1e+10)
      weightsptr[k] = 1e+10;
    if(nW > 1)
      weightsptr[k] *= Wptr[k];

    if(ISNA(scoreptr[k]))
      scoreptr[k] = 1.490116e-08;
    if(scoreptr[k] < -1e+10)
      scoreptr[k] = -1e+10;
    if(scoreptr[k] > 1e+10)
      scoreptr[k] = 1e+10;

    zptr[k] = etaptr[k] + scoreptr[k] / weightsptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];

    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
  }

  /* Compute X'WX. */
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, &nc, &nc, &nr, &one,
    XWptr, &nc, Xptr, &nr, &zero, XWXptr, &nc);

  SEXP XWX0;
  PROTECT(XWX0 = duplicate(getListElement(x, "XWX")));
  double *XWX0ptr = REAL(XWX0);
  ++nProtected;

  /* Add penalty matrix and variance parameter. */
  if(fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      tau2ptr[jj] = thetaptr[nc + jj];
      for(i = 0; i < nc; i++) {
        for(j = 0; j < nc; j++) {
          XWXptr[i + nc * j] += 1 / tau2ptr[jj] * Sptr[i + nc * j];
        }
      }
    }
  }

  /* Cholesky decompostion of XWX. */
  SEXP L;
  PROTECT(L = duplicate(getListElement(x, "XWX")));
  double *Lptr = REAL(L);
  ++nProtected;

  for(j = 0; j < nc; j++) { 	/* Zero the lower triangle. */
    for(i = j + 1; i < nc; i++) {
      Lptr[i + nc * j] = 0.0;
    }
  }

  int info;
  F77_CALL(dpotrf)("Upper", &nc, Lptr, &nc, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV;
  PROTECT(PINV = duplicate(L));
  double *PINVptr = REAL(PINV);
  ++nProtected;

  F77_CALL(dpotri)("Upper", &nc, PINVptr, &nc, &info);

  SEXP PINVL;
  PROTECT(PINVL = duplicate(PINV));
  double *PINVLptr = REAL(PINVL);
  ++nProtected;
  F77_CALL(dpotrf)("Upper", &nc, PINVLptr, &nc, &info);

  for(j = 0; j < nc; j++) {
    for(i = j + 1; i < nc; i++) {
      PINVptr[i + j * nc] = PINVptr[j + i * nc];
    }
  }

  /* Compute mu. */
  SEXP mu0;
  PROTECT(mu0 = allocVector(REALSXP, nc));
  double *mu0ptr = REAL(mu0);
  ++nProtected;

  SEXP mu1;
  PROTECT(mu1 = allocVector(REALSXP, nc));
  double *mu1ptr = REAL(mu1);
  ++nProtected;

  int k1 = 1;
  char *transa2 = "T";
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nr, &one,
    Xptr, &nr, xrresptr, &nr, &zero, mu0ptr, &nc);
  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,
    PINVptr, &nc, mu0ptr, &nc, &zero, mu1ptr, &nc);

  /* Sample. */
  double edf1 = 0.0;
  double edf2 = 0.0;
  GetRNGstate();
  for(j = 0; j < nc; j++) {
    gamma0ptr[j] = rnorm(0, 1);
    for(i = j; i < nc; i++) {
      for(k = 0; k < nc; k++) {
        edf1 = edf1 + XWX0ptr[j + k * nc] * PINVptr[k + i * nc];
      }
      if(j == i)
        edf2 += edf1;
      edf1 = 0.0;
    }
  }
  PutRNGstate();
 
  SEXP edf;
  PROTECT(edf = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(edf)[0] = edf2;
  
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nc, &one,
    PINVLptr, &nc, gamma0ptr, &nc, &zero, gamma1ptr, &nc);

  double sdiag0 = 0.0;
  for(j = 0; j < nc; j++) {
    gamma1ptr[j] += mu1ptr[j];
    sdiag0 += log(pow(Lptr[j + nc * j], 2.0));
  }

  /* Log priors. */
  double p1 = 0.0;
  double p2 = 0.0;
  double qbetaprop = 0.0;
  double tsum1 = 0.0;
  double tsum2 = 0.0;
  double gSg = 0.0;
  double a = 0.0;
  double b = 0.0;
  if(ntau2 < 2) {
    if(fixed < 1) {
      a = REAL(getListElement(x, "a"))[0];
      b = REAL(getListElement(x, "b"))[0];
      double *rankptr = REAL(getListElement(x, "rank"));
      for(jj = 0; jj < ntau2; jj++) {
        Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
        gSg = 0.0;
        for(i = 0; i < nc; i++) {
          tsum1 = 0.0;
          tsum2 = 0.0;
          for(j = 0; j < nc; j++) {
            tsum1 += gamma1ptr[j] * Sptr[j + i * nc];
            tsum2 += (gamma1ptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
          }
          gSg += tsum1 * gamma1ptr[i];
          qbetaprop += tsum2 * (gamma1ptr[i] - mu1ptr[i]);
        }
        if(ntau2 < 2) {
          p2 += -log(tau2ptr[jj]) * rankptr[jj] / 2 - 0.5 / tau2ptr[jj] * gSg +
            log(pow(b, a)) - exp(lgamma(a)) + (-a - 1) * log(tau2ptr[jj]) - b / tau2ptr[jj];
        }
      }
    } else {
      for(i = 0; i < nc; i++) {
        tsum1 = 0.0;
        p2 += dnorm(gamma1ptr[i], 0, 1000, 1);
        for(j = 0; j < nc; j++) {
          tsum1 += (gamma1ptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
        }
        qbetaprop += tsum1 * (gamma1ptr[i] - mu1ptr[i]);
      }
    }
    qbetaprop = 0.5 * sdiag0 - 0.5 * qbetaprop;
  } else {
    SEXP P1;
    PROTECT(P1 = eval_prior(getListElement(x, "prior"), theta2, rho));
    ++nProtected;
    p1 = REAL(P1)[0];

    SEXP QBETAPROP;
    PROTECT(QBETAPROP = eval_dmvnorm_log(getListElement(x, "dmvnorm_log"), gamma1, mu1, PINV, rho));
    ++nProtected;
    qbetaprop = REAL(QBETAPROP)[0];
  }

  /* Part 2. */
  /* Obtain new fitted values and update predictor. */
  F77_CALL(dgemm)(transa, transb, &nr, &k1, &nc, &one,
    Xptr, &nr, gamma1ptr, &nr, &zero, fitrptr, &nr);

  for(i = 0; i < n; i++) {
    k = idptr[i] - 1;
    fitptr[i] = fitrptr[k];
    etaptr[i] += fitptr[i];
  }

  /* Evaluate loglik, weights and score vector. */
  peta = map2par(getListElement(family, "map2par"), eta2, rho);
  double pibetaprop = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, id2, rho))[0];

  SEXP weights2;
  PROTECT(weights2 = iwls_eval(getListElement(getListElement(family, "hess"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *weights2ptr = REAL(weights2);
  ++nProtected;

  SEXP score2;
  PROTECT(score2 = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *score2ptr = REAL(score2);
  ++nProtected;

  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  j = 0;
  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      for(jj = 0; jj < nc; jj++) {
        XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
      }
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }
    k = orderptr[i] - 1;

    if(ISNA(weights2ptr[k]))
      weights2ptr[k] = 1.490116e-08;
    if(weights2ptr[k] < 1e-10)
      weights2ptr[k] = 1e-10;
    if(weights2ptr[k] < 0.0)
      weights2ptr[k] = -1.0 * weights2ptr[k];
    if(weights2ptr[k] > 1e+10)
      weights2ptr[k] = 1e+10;
    if(nW > 1)
      weights2ptr[k] *= Wptr[k];

    if(ISNA(score2ptr[k]))
      score2ptr[k] = 1.490116e-08;
    if(score2ptr[k] < -1e+10)
      score2ptr[k] = -1e+10;
    if(score2ptr[k] > 1e+10)
      score2ptr[k] = 1e+10;

    zptr[k] = etaptr[k] + score2ptr[k] / weights2ptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];
    xweightsptr[j] += weights2ptr[k];
    xrresptr[j] += weights2ptr[k] * eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
  }

  /* Compute X'WX. */
  F77_CALL(dgemm)(transa, transb, &nc, &nc, &nr, &one,
    XWptr, &nc, Xptr, &nr, &zero, XWXptr, &nc);

  /* Add penalty matrix and variance parameter. */
  if(fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      for(i = 0; i < nc; i++) {
        for(j = 0; j < nc; j++) {
          XWXptr[i + nc * j] += 1 / tau2ptr[jj] * Sptr[i + nc * j];
        }
      }
    }
  }

  /* Cholesky decompostion of XWX. */
  SEXP L2;
  PROTECT(L2 = duplicate(getListElement(x, "XWX")));
  Lptr = REAL(L2);
  ++nProtected;

  for(j = 0; j < nc; j++) { 	/* Zero the lower triangle. */
    for(i = j + 1; i < nc; i++) {
      Lptr[i + nc * j] = 0.0;
    }
  }

  F77_CALL(dpotrf)("Upper", &nc, Lptr, &nc, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV2;
  PROTECT(PINV2 = duplicate(L2));
  PINVptr = REAL(PINV2);
  ++nProtected;

  F77_CALL(dpotri)("Upper", &nc, PINVptr, &nc, &info);
  F77_CALL(dpotrf)("Upper", &nc, PINVLptr, &nc, &info);

  sdiag0 = 0.0;
  for(j = 0; j < nc; j++) {
    sdiag0 += log(pow(Lptr[j + nc * j], 2.0));
    for(i = j + 1; i < nc; i++) {
      PINVptr[i + j * nc] = PINVptr[j + i * nc];
    }
  }

  /* Compute mu. */
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nr, &one,
    Xptr, &nr, xrresptr, &nr, &zero, mu0ptr, &nc);

  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,
    PINVptr, &nc, mu0ptr, &nc, &zero, mu1ptr, &nc);

  /* Log priors. */
  double qbeta = 0.0;
  if(ntau2 < 2) {
    double a2 = 0.0;
    double b2 = 0.0;
    if(fixed < 1) {
      double *rankptr2 = REAL(getListElement(x, "rank"));
      for(jj = 0; jj < ntau2; jj++) {
        Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
        gSg = 0.0;
        for(i = 0; i < nc; i++) {
          tsum1 = 0.0;
          tsum2 = 0.0;
          for(j = 0; j < nc; j++) {
            tsum1 += thetaptr[j] * Sptr[j + i * nc];
            tsum2 += (thetaptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
          }
          gSg += tsum1 * thetaptr[i];
          qbeta += tsum2 * (thetaptr[i] - mu1ptr[i]);
        }
        if(ntau2 < 2) {
          p1 += -log(tau2ptr[jj]) * rankptr2[jj] / 2 + -0.5 / tau2ptr[jj] * gSg +
            log(pow(b, a)) - exp(lgamma(a)) + (-a - 1) * log(tau2ptr[jj]) - b / tau2ptr[jj];
        }
        if((fxsp < 1) && (ntau2 < 2)) {
          a2 = rankptr2[jj] / 2 + a;
          b2 = 0.5 * gSg + b;
          GetRNGstate();
          tau2ptr[jj] = 1 / rgamma(a2, 1 / b2);
          PutRNGstate();
        }
      }
    } else {
      for(i = 0; i < nc; i++) {
        tsum1 = 0.0;
        p1 += dnorm(thetaptr[i], 0, 1000, 1);
        for(j = 0; j < nc; j++) {
          tsum1 += (thetaptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
        }
        qbeta += tsum1 * (thetaptr[i] - mu1ptr[i]);
      }
    }
    qbeta = 0.5 * sdiag0 - 0.5 * qbeta;
  }

  double cval = 0.0;
  for(j = 0; j < nc; j++) {
    cval = thetaptr[j];
    thetaptr[j] = gamma1ptr[j];
    gamma1ptr[j] = cval;
  }
  if((fxsp < 1) && (fixed < 1) && (ntau2 < 2)) {
    for(jj = 0; jj < ntau2; jj++) {
      thetaptr[nc + jj] = tau2ptr[jj];
    }
  }

  if(ntau2 > 1) {
    SEXP P2;
    PROTECT(P2 = eval_prior(getListElement(x, "prior"), theta2, rho));
    ++nProtected;
    p2 = REAL(P2)[0];

    SEXP QBETA;
    PROTECT(QBETA = eval_dmvnorm_log(getListElement(x, "dmvnorm_log"), gamma1, mu1, PINV, rho));
    ++nProtected;
    qbeta = REAL(QBETA)[0];
  }

  SEXP alpha;
  PROTECT(alpha = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(alpha)[0] = (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1);

  SEXP loglik;
  PROTECT(loglik = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(loglik)[0] = pibetaprop;

/*Rprintf("pibetaprop %g\n", pibetaprop);*/
/*Rprintf("qbeta %g\n", qbeta);*/
/*Rprintf("p2 %g\n", p2);*/
/*Rprintf("pibeta %g\n", pibeta);*/
/*Rprintf("qbetaprop %g\n", qbetaprop);*/
/*Rprintf("p1 %g\n", p1);*/
/*Rprintf("alpha %g\n", exp(REAL(alpha)[0]));*/

  /* Stuff everything together. */
  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 4));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, theta2);
  SET_VECTOR_ELT(rval, 1, alpha);
  SET_VECTOR_ELT(rval, 2, edf);
  SET_VECTOR_ELT(rval, 3, loglik);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 4));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("parameters"));
  SET_STRING_ELT(nrval, 1, mkChar("alpha"));
  SET_STRING_ELT(nrval, 2, mkChar("edf"));
  SET_STRING_ELT(nrval, 3, mkChar("loglik"));
        
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(nProtected);
  return rval;
}


SEXP gmcmc_iwls_gp(SEXP family, SEXP theta, SEXP id,
  SEXP eta, SEXP response, SEXP x, SEXP z, SEXP e, SEXP id2, SEXP W, SEXP rho)
{
  int i, j, k, nProtected = 0;
  int n = INTEGER(getListElement(x, "nobs"))[0];
  int fixed = LOGICAL(getListElement(x, "fixed"))[0];
  int fxsp = LOGICAL(getListElement(x, "fxsp"))[0];
  int nW = length(W);
  if(nW > 1) {
    if(nW != n)
      nW = 1;
  }
  double *Wptr = REAL(W);

  SEXP theta2;
  PROTECT(theta2 = duplicate(getListElement(getListElement(theta,
    CHAR(STRING_ELT(id, 0))), CHAR(STRING_ELT(id, 1)))));
  double *thetaptr = REAL(theta2);
  ++nProtected;

  int S_ind = getListElement_index(x, "S");
  int *penFun = INTEGER(getListElement(x, "penaltyFunction"));
  int ntau2;
  if(fixed > 0) {
    ntau2 = 0;
  } else {
    ntau2 = length(VECTOR_ELT(x, S_ind));
  }
  int nc = length(theta2);
  if(fixed < 1) {
    nc -= ntau2;
  }

  SEXP gamma0, gamma1, tau2;
  PROTECT(gamma0 = allocVector(REALSXP, nc));
  double *gamma0ptr = REAL(gamma0);
  ++nProtected;

  PROTECT(gamma1 = allocVector(REALSXP, nc));
  double *gamma1ptr = REAL(gamma1);
  ++nProtected;

  PROTECT(tau2 = allocVector(REALSXP, ntau2));
  double *tau2ptr = REAL(tau2);
  ++nProtected;

  /* Evaluate loglik, weights and score vector. */
  SEXP peta;
  PROTECT(peta = map2par(getListElement(family, "map2par"), eta, rho));
  ++nProtected;
  int ll_ind = getListElement_index(family, "loglik");
  double pibeta = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, id2, rho))[0];

  SEXP weights;
  PROTECT(weights = iwls_eval(getListElement(getListElement(family, "hess"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *weightsptr = REAL(weights);
  ++nProtected;

  SEXP score;
  PROTECT(score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *scoreptr = REAL(score);
  ++nProtected;

  SEXP eta2;
  PROTECT(eta2 = duplicate(eta));
  double *etaptr = REAL(getListElement(eta2, CHAR(STRING_ELT(id, 0))));
  ++nProtected;

  /* Create weighted matrix */
  int X_ind = getListElement_index(x, "X");
  int nr = nrows(VECTOR_ELT(x, X_ind));

  /* More pointers needed. */
  double *zptr = REAL(z);
  double *eptr = REAL(e);
  double *xweightsptr = REAL(getListElement(x, "weights"));
  double *xrresptr = REAL(getListElement(x, "rres"));
  double *XWptr = REAL(getListElement(x, "XW"));
  double *XWXptr = REAL(getListElement(x, "XWX"));
  double *Xptr = REAL(VECTOR_ELT(x, X_ind));
  double *Sptr;
  int *idptr = INTEGER(getListElement(getListElement(x, "binning"), "match.index"));
  int *indptr = INTEGER(getListElement(getListElement(x, "binning"), "sorted.index"));
  int *orderptr = INTEGER(getListElement(getListElement(x, "binning"), "order"));

  /* Handling fitted.values. */
  SEXP fitname;
  PROTECT(fitname = allocVector(STRSXP, 1));
  ++nProtected;
  SET_STRING_ELT(fitname, 0, mkChar("fitted.values"));
  double *fitptr = REAL(getAttrib(theta2, fitname));
  double *fitrptr = REAL(getListElement(x, "fit.reduced"));

  /* Start. */
  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  j = 0;
  int jj;
  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      for(jj = 0; jj < nc; jj++) {
        XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
      }
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }
    k = orderptr[i] - 1;

    if(ISNA(weightsptr[k]))
      weightsptr[k] = 1.490116e-08;
    if(weightsptr[k] < 0.0)
      weightsptr[k] = -1.0 * weightsptr[k];
    if(weightsptr[k] > 1e+10)
      weightsptr[k] = 1e+10;
    if(nW > 1)
      weightsptr[k] *= Wptr[k];

    if(ISNA(scoreptr[k]))
      scoreptr[k] = 1.490116e-08;
    if(scoreptr[k] < -1e+10)
      scoreptr[k] = -1e+10;
    if(scoreptr[k] > 1e+10)
      scoreptr[k] = 1e+10;

    zptr[k] = etaptr[k] + scoreptr[k] / weightsptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];

    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
  }

  /* Compute X'WX. */
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, &nc, &nc, &nr, &one,
    XWptr, &nc, Xptr, &nr, &zero, XWXptr, &nc);

  SEXP XWX0;
  PROTECT(XWX0 = duplicate(getListElement(x, "XWX")));
  double *XWX0ptr = REAL(XWX0);
  ++nProtected;

  /* Add penalty matrix and variance parameter. */
  if(fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      if(penFun[jj] > 0) {
        Sptr = REAL(get_S_mat(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj), theta2, getListElement(x, "fixed.hyper"), rho));
      } else {
        Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      }
      tau2ptr[jj] = thetaptr[nc + jj];
      for(i = 0; i < nc; i++) {
        for(j = 0; j < nc; j++) {
          XWXptr[i + nc * j] += 1 / tau2ptr[jj] * Sptr[i + nc * j];
        }
      }
    }
  }

  /* Cholesky decompostion of XWX. */
  SEXP L;
  PROTECT(L = duplicate(getListElement(x, "XWX")));
  double *Lptr = REAL(L);
  ++nProtected;

  for(j = 0; j < nc; j++) { 	/* Zero the lower triangle. */
    for(i = j + 1; i < nc; i++) {
      Lptr[i + nc * j] = 0.0;
    }
  }

  int info;
  F77_CALL(dpotrf)("Upper", &nc, Lptr, &nc, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV;
  PROTECT(PINV = duplicate(L));
  double *PINVptr = REAL(PINV);
  ++nProtected;

  F77_CALL(dpotri)("Upper", &nc, PINVptr, &nc, &info);

  SEXP PINVL;
  PROTECT(PINVL = duplicate(PINV));
  double *PINVLptr = REAL(PINVL);
  ++nProtected;
  F77_CALL(dpotrf)("Upper", &nc, PINVLptr, &nc, &info);

  for(j = 0; j < nc; j++) {
    for(i = j + 1; i < nc; i++) {
      PINVptr[i + j * nc] = PINVptr[j + i * nc];
    }
  }

  /* Compute mu. */
  SEXP mu0;
  PROTECT(mu0 = allocVector(REALSXP, nc));
  double *mu0ptr = REAL(mu0);
  ++nProtected;

  SEXP mu1;
  PROTECT(mu1 = allocVector(REALSXP, nc));
  double *mu1ptr = REAL(mu1);
  ++nProtected;

  int k1 = 1;
  char *transa2 = "T";
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nr, &one,
    Xptr, &nr, xrresptr, &nr, &zero, mu0ptr, &nc);
  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,
    PINVptr, &nc, mu0ptr, &nc, &zero, mu1ptr, &nc);

  /* Sample. */
  double edf1 = 0.0;
  double edf2 = 0.0;
  GetRNGstate();
  for(j = 0; j < nc; j++) {
    gamma0ptr[j] = rnorm(0, 1);
    for(i = j; i < nc; i++) {
      for(k = 0; k < nc; k++) {
        edf1 = edf1 + XWX0ptr[j + k * nc] * PINVptr[k + i * nc];
      }
      if(j == i)
        edf2 += edf1;
      edf1 = 0.0;
    }
  }
  PutRNGstate();
 
  SEXP edf;
  PROTECT(edf = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(edf)[0] = edf2;
  
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nc, &one,
    PINVLptr, &nc, gamma0ptr, &nc, &zero, gamma1ptr, &nc);

  for(j = 0; j < nc; j++) {
    gamma1ptr[j] += mu1ptr[j];
  }

  /* Log priors. */
  SEXP P1;
  PROTECT(P1 = eval_prior(getListElement(x, "prior"), theta2, rho));
  ++nProtected;
  double p1 = REAL(P1)[0];

  SEXP QBETAPROP;
  PROTECT(QBETAPROP = eval_dmvnorm_log(getListElement(x, "dmvnorm_log"), gamma1, mu1, PINV, rho));
  ++nProtected;
  double qbetaprop = REAL(QBETAPROP)[0];

  /* Part 2. */
  /* Obtain new fitted values and update predictor. */
  F77_CALL(dgemm)(transa, transb, &nr, &k1, &nc, &one,
    Xptr, &nr, gamma1ptr, &nr, &zero, fitrptr, &nr);

  for(i = 0; i < n; i++) {
    k = idptr[i] - 1;
    fitptr[i] = fitrptr[k];
    etaptr[i] += fitptr[i];
  }

  /* Evaluate loglik, weights and score vector. */
  peta = map2par(getListElement(family, "map2par"), eta2, rho);
  double pibetaprop = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, id2, rho))[0];

  SEXP weights2;
  PROTECT(weights2 = iwls_eval(getListElement(getListElement(family, "hess"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *weights2ptr = REAL(weights2);
  ++nProtected;

  SEXP score2;
  PROTECT(score2 = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *score2ptr = REAL(score2);
  ++nProtected;

  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  j = 0;
  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      for(jj = 0; jj < nc; jj++) {
        XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
      }
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }
    k = orderptr[i] - 1;

    if(ISNA(weights2ptr[k]))
      weights2ptr[k] = 1.490116e-08;
    if(weights2ptr[k] < 0.0)
      weights2ptr[k] = -1.0 * weights2ptr[k];
    if(weights2ptr[k] > 1e+10)
      weights2ptr[k] = 1e+10;
    if(nW > 1)
      weightsptr[k] *= Wptr[k];

    if(ISNA(score2ptr[k]))
      score2ptr[k] = 1.490116e-08;
    if(score2ptr[k] < -1e+10)
      score2ptr[k] = -1e+10;
    if(score2ptr[k] > 1e+10)
      score2ptr[k] = 1e+10;

    zptr[k] = etaptr[k] + score2ptr[k] / weights2ptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];
    xweightsptr[j] += weights2ptr[k];
    xrresptr[j] += weights2ptr[k] * eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
  }

  /* Compute X'WX. */
  F77_CALL(dgemm)(transa, transb, &nc, &nc, &nr, &one,
    XWptr, &nc, Xptr, &nr, &zero, XWXptr, &nc);

  /* Add penalty matrix and variance parameter. */
  if(fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      if(penFun[jj] > 0) {
        Sptr = REAL(get_S_mat(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj), gamma1, getListElement(x, "fixed.hyper"), rho));
      } else {
        Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      }
      for(i = 0; i < nc; i++) {
        for(j = 0; j < nc; j++) {
          XWXptr[i + nc * j] += 1 / tau2ptr[jj] * Sptr[i + nc * j];
        }
      }
    }
  }

  /* Cholesky decompostion of XWX. */
  SEXP L2;
  PROTECT(L2 = duplicate(getListElement(x, "XWX")));
  Lptr = REAL(L2);
  ++nProtected;

  for(j = 0; j < nc; j++) { 	/* Zero the lower triangle. */
    for(i = j + 1; i < nc; i++) {
      Lptr[i + nc * j] = 0.0;
    }
  }

  F77_CALL(dpotrf)("Upper", &nc, Lptr, &nc, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV2;
  PROTECT(PINV2 = duplicate(L2));
  PINVptr = REAL(PINV2);
  ++nProtected;

  F77_CALL(dpotri)("Upper", &nc, PINVptr, &nc, &info);
  F77_CALL(dpotrf)("Upper", &nc, PINVLptr, &nc, &info);

  for(j = 0; j < nc; j++) {
    for(i = j + 1; i < nc; i++) {
      PINVptr[i + j * nc] = PINVptr[j + i * nc];
    }
  }

  /* Compute mu. */
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nr, &one,
    Xptr, &nr, xrresptr, &nr, &zero, mu0ptr, &nc);

  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,
    PINVptr, &nc, mu0ptr, &nc, &zero, mu1ptr, &nc);

  /* Log priors. */
  double cval = 0.0;
  for(j = 0; j < nc; j++) {
    cval = thetaptr[j];
    thetaptr[j] = gamma1ptr[j];
    gamma1ptr[j] = cval;
  }
  if((fxsp < 1) && (fixed < 1) && (ntau2 < 2)) {
    for(jj = 0; jj < ntau2; jj++) {
      thetaptr[nc + jj] = tau2ptr[jj];
    }
  }

  SEXP P2;
  PROTECT(P2 = eval_prior(getListElement(x, "prior"), theta2, rho));
  ++nProtected;
  double p2 = REAL(P2)[0];

  SEXP QBETA;
  PROTECT(QBETA = eval_dmvnorm_log(getListElement(x, "dmvnorm_log"), gamma1, mu1, PINV, rho));
  ++nProtected;
  double qbeta = REAL(QBETA)[0];

  SEXP alpha;
  PROTECT(alpha = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(alpha)[0] = (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1);

  SEXP loglik;
  PROTECT(loglik = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(loglik)[0] = pibetaprop;

  /* Stuff everything together. */
  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 4));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, theta2);
  SET_VECTOR_ELT(rval, 1, alpha);
  SET_VECTOR_ELT(rval, 2, edf);
  SET_VECTOR_ELT(rval, 3, loglik);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 4));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("parameters"));
  SET_STRING_ELT(nrval, 1, mkChar("alpha"));
  SET_STRING_ELT(nrval, 2, mkChar("edf"));
  SET_STRING_ELT(nrval, 3, mkChar("loglik"));
        
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(nProtected);
  return rval;
}


SEXP gmcmc_iwls_gp_diag_lasso(SEXP family, SEXP theta, SEXP id,
  SEXP eta, SEXP response, SEXP x, SEXP z, SEXP e, SEXP id2, SEXP W, SEXP rho)
{
  int i, j, k, nProtected = 0;
  int n = INTEGER(getListElement(x, "nobs"))[0];
  int nW = length(W);
  if(nW > 1) {
    if(nW != n)
      nW = 1;
  }
  double *Wptr = REAL(W);

  SEXP theta2;
  PROTECT(theta2 = duplicate(getListElement(getListElement(theta,
    CHAR(STRING_ELT(id, 0))), CHAR(STRING_ELT(id, 1)))));
  double *thetaptr = REAL(theta2);
  ++nProtected;

  int S_ind = getListElement_index(x, "S");
  int ntau2 = length(VECTOR_ELT(x, S_ind));
  int nc = length(theta2) - ntau2;

  int *iptr = INTEGER(getListElement(getListElement(x, "sparse.setup"), "matrix"));
  int nc_index = ncols(getListElement(getListElement(x, "sparse.setup"), "matrix"));

  double *dfptr = REAL(getListElement(getListElement(x, "lasso"), "df"));

  SEXP gamma0, gamma1;
  PROTECT(gamma0 = allocVector(REALSXP, nc));
  ++nProtected;

  PROTECT(gamma1 = allocVector(REALSXP, nc));
  double *gamma1ptr = REAL(gamma1);
  ++nProtected;

  /* Evaluate loglik, weights and score vector. */
  SEXP peta;
  PROTECT(peta = map2par(getListElement(family, "map2par"), eta, rho));
  ++nProtected;
  int ll_ind = getListElement_index(family, "loglik");
  double pibeta = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, id2, rho))[0];

  SEXP weights;
  PROTECT(weights = iwls_eval(getListElement(getListElement(family, "hess"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *weightsptr = REAL(weights);
  ++nProtected;

  SEXP score;
  PROTECT(score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *scoreptr = REAL(score);
  ++nProtected;

  SEXP eta2;
  PROTECT(eta2 = duplicate(eta));
  double *etaptr = REAL(getListElement(eta2, CHAR(STRING_ELT(id, 0))));
  ++nProtected;

  /* Create weighted matrix */
  int X_ind = getListElement_index(x, "X");
  int nr = nrows(VECTOR_ELT(x, X_ind));

  /* More pointers needed. */
  double *zptr = REAL(z);
  double *eptr = REAL(e);
  double *xweightsptr = REAL(getListElement(x, "weights"));
  double *xrresptr = REAL(getListElement(x, "rres"));
  double *XWXptr = REAL(getListElement(x, "XWX"));
  double *Xptr = REAL(VECTOR_ELT(x, X_ind));
  int *idptr = INTEGER(getListElement(getListElement(x, "binning"), "match.index"));
  int *indptr = INTEGER(getListElement(getListElement(x, "binning"), "sorted.index"));
  int *orderptr = INTEGER(getListElement(getListElement(x, "binning"), "order"));

  /* Handling fitted.values. */
  SEXP fitname;
  PROTECT(fitname = allocVector(STRSXP, 1));
  ++nProtected;
  SET_STRING_ELT(fitname, 0, mkChar("fitted.values"));
  double *fitptr = REAL(getAttrib(theta2, fitname));
  double *fitrptr = REAL(getListElement(x, "fit.reduced"));

  /* Init mu. */
  SEXP mu0;
  PROTECT(mu0 = allocVector(REALSXP, nc));
  double *mu0ptr = REAL(mu0);
  ++nProtected;

  SEXP mu1;
  PROTECT(mu1 = allocVector(REALSXP, nc));
  double *mu1ptr = REAL(mu1);
  ++nProtected;

  /* Start. */
  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  j = 0;
  int jj;

  for(jj = 0; jj < nc; jj++) {
    XWXptr[jj + nc * jj] = 0.0;
    mu0ptr[jj] = 0.0;
    mu1ptr[jj] = 0.0;
  }

  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {

/*      for(jj = 0; jj < nc_index; jj++) {*/
/*        if((iptr[j + jj * nr] < 0) || (iptr[j + jj * nr] > nc))*/
/*          continue;*/
/*        XWXptr[jj + nc * jj] += pow(Xptr[j + (iptr[j + jj * nr] - 1) * nr], 2.0) * xweightsptr[j];*/
/*        mu0ptr[jj] += Xptr[j + (iptr[j + jj * nr] - 1) * nr] * xrresptr[j];*/
/*      }*/

      for(jj = 0; jj < nc; jj++) {
        XWXptr[jj + nc * jj] += pow(Xptr[j + nr * jj], 2.0) * xweightsptr[j];
        mu0ptr[jj] += Xptr[j + nr * jj] * xrresptr[j];
      }
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }
    k = orderptr[i] - 1;

    if(ISNA(weightsptr[k]))
      weightsptr[k] = 1.490116e-08;
    if(weightsptr[k] < 0.0)
      weightsptr[k] = -1.0 * weightsptr[k];
    if(weightsptr[k] > 1e+10)
      weightsptr[k] = 1e+10;
    if(nW > 1)
      weightsptr[k] *= Wptr[k];

    if(ISNA(scoreptr[k]))
      scoreptr[k] = 1.490116e-08;
    if(scoreptr[k] < -1e+10)
      scoreptr[k] = -1e+10;
    if(scoreptr[k] > 1e+10)
      scoreptr[k] = 1e+10;

    zptr[k] = etaptr[k] + scoreptr[k] / weightsptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];

    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWXptr[jj + nc * jj] += pow(Xptr[j + nr * jj], 2.0) * xweightsptr[j];
    mu0ptr[jj] += Xptr[j + nr * jj] * xrresptr[j];
  }

/*  for(jj = 0; jj < nc_index; jj++) {*/
/*    if((iptr[j + jj * nr] < 0) || (iptr[j + jj * nr] > nc))*/
/*      continue;*/
/*    XWXptr[jj + nc * jj] += pow(Xptr[j + (iptr[j + jj * nr] - 1) * nr], 2.0) * xweightsptr[j];*/
/*    mu0ptr[jj] += Xptr[j + (iptr[j + jj * nr] - 1) * nr] * xrresptr[j];*/
/*  }*/

  double edf0 = 0.0;
  double edf1 = 0.0;
  double qbetaprop = 0.0;
  GetRNGstate();
  if(ntau2 > 1) {
    for(j = 0; j < nc; j++) {
      edf1 = XWXptr[j + j * nc];
      XWXptr[j + nc * j] += dfptr[j] / (thetaptr[nc + j] * pow(pow(thetaptr[j], 2.0) + 0.000001, 0.5));
      XWXptr[j + nc * j] = 1.0 / XWXptr[j + nc * j];
      mu1ptr[j] = XWXptr[j + nc * j] * mu0ptr[j];
      gamma1ptr[j] = rnorm(mu1ptr[j], pow(XWXptr[j + j * nc], 0.5));
      edf0 = edf0 + edf1 * XWXptr[j + j * nc];
      qbetaprop += dnorm(gamma1ptr[j], mu1ptr[j], pow(XWXptr[j + nc * j], 0.5), 1);
    }
  } else {
    for(j = 0; j < nc; j++) {
      edf1 = XWXptr[j + j * nc];
      XWXptr[j + nc * j] += dfptr[j] / (thetaptr[nc] * pow(pow(thetaptr[j], 2.0) + 0.000001, 0.5));
      XWXptr[j + nc * j] = 1.0 / XWXptr[j + nc * j];
      mu1ptr[j] = XWXptr[j + nc * j] * mu0ptr[j];
      gamma1ptr[j] = rnorm(mu1ptr[j], pow(XWXptr[j + j * nc], 0.5));
      edf0 = edf0 + edf1 * XWXptr[j + j * nc];
      qbetaprop += dnorm(gamma1ptr[j], mu1ptr[j], pow(XWXptr[j + nc * j], 0.5), 1);
    }
  }
  PutRNGstate();
 
  SEXP edf;
  PROTECT(edf = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(edf)[0] = edf0;

  SEXP P1;
  PROTECT(P1 = eval_prior(getListElement(x, "prior"), theta2, rho));
  ++nProtected;
  double p1 = REAL(P1)[0];

  for(i = 0; i < nr; i++) {
    fitrptr[i] = 0.0;
    for(j = 0; j < nc_index; j++) {
      if((iptr[i + j * nr] < 0) || (iptr[i + j * nr] > nc))
        continue;
      fitrptr[i] += Xptr[i + (iptr[i + j * nr] - 1) * nr] * gamma1ptr[(iptr[i + j * nr] - 1)];
    }
  }

  for(i = 0; i < n; i++) {
    k = idptr[i] - 1;
    fitptr[i] = fitrptr[k];
    etaptr[i] += fitptr[i];
  }

  peta = map2par(getListElement(family, "map2par"), eta2, rho);
  double pibetaprop = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, id2, rho))[0];

  SEXP weights2;
  PROTECT(weights2 = iwls_eval(getListElement(getListElement(family, "hess"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *weights2ptr = REAL(weights2);
  ++nProtected;

  SEXP score2;
  PROTECT(score2 = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, id2, rho));
  double *score2ptr = REAL(score2);
  ++nProtected;

  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  for(jj = 0; jj < nc; jj++) {
    XWXptr[jj + nc * jj] = 0.0;
    mu0ptr[jj] = 0.0;
    mu1ptr[jj] = 0.0;
  }

  j = 0;
  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {

/*      for(jj = 0; jj < nc_index; jj++) {*/
/*        if((iptr[j + jj * nr] < 0) || (iptr[j + jj * nr] > nc))*/
/*          continue;*/
/*        XWXptr[jj + nc * jj] += pow(Xptr[j + (iptr[j + jj * nr] - 1) * nr], 2.0) * xweightsptr[j];*/
/*        mu0ptr[jj] += Xptr[j + (iptr[j + jj * nr] - 1) * nr] * xrresptr[j];*/
/*      }*/

      for(jj = 0; jj < nc; jj++) {
        XWXptr[jj + nc * jj] += pow(Xptr[j + nr * jj], 2.0) * xweightsptr[j];
        mu0ptr[jj] += Xptr[j + nr * jj] * xrresptr[j];
      }
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }
    k = orderptr[i] - 1;

    if(ISNA(weights2ptr[k]))
      weights2ptr[k] = 1.490116e-08;
    if(weights2ptr[k] < 0.0)
      weights2ptr[k] = -1.0 * weights2ptr[k];
    if(weights2ptr[k] > 1e+10)
      weights2ptr[k] = 1e+10;
    if(nW > 1)
      weightsptr[k] *= Wptr[k];

    if(ISNA(score2ptr[k]))
      score2ptr[k] = 1.490116e-08;
    if(score2ptr[k] < -1e+10)
      score2ptr[k] = -1e+10;
    if(score2ptr[k] > 1e+10)
      score2ptr[k] = 1e+10;

    zptr[k] = etaptr[k] + score2ptr[k] / weights2ptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];
    xweightsptr[j] += weights2ptr[k];
    xrresptr[j] += weights2ptr[k] * eptr[k];
  }

/*  for(jj = 0; jj < nc_index; jj++) {*/
/*    if((iptr[j + jj * nr] < 0) || (iptr[j + jj * nr] > nc))*/
/*      continue;*/
/*    XWXptr[jj + nc * jj] += pow(Xptr[j + (iptr[j + jj * nr] - 1) * nr], 2.0) * xweightsptr[j];*/
/*    mu0ptr[jj] += Xptr[j + (iptr[j + jj * nr] - 1) * nr] * xrresptr[j];*/
/*  }*/

  for(jj = 0; jj < nc; jj++) {
    XWXptr[jj + nc * jj] += pow(Xptr[j + nr * jj], 2.0) * xweightsptr[j];
    mu0ptr[jj] += Xptr[j + nr * jj] * xrresptr[j];
  }

  double qbeta = 0.0;
  if(ntau2 > 1) {
    for(j = 0; j < nc; j++) {
      XWXptr[j + nc * j] += dfptr[j] / (thetaptr[nc + j] * pow(pow(thetaptr[j], 2.0) + 0.000001, 0.5));
      XWXptr[j + nc * j] = 1.0 / XWXptr[j + nc * j];
      mu1ptr[j] = XWXptr[j + nc * j] * mu0ptr[j];
      qbeta += dnorm(thetaptr[j], mu1ptr[j], pow(XWXptr[j + nc * j], 0.5), 1);
      thetaptr[j] = gamma1ptr[j];
    }
  } else {
    for(j = 0; j < nc; j++) {
      XWXptr[j + nc * j] += dfptr[j] / (thetaptr[nc] * pow(pow(thetaptr[j], 2.0) + 0.000001, 0.5));
      XWXptr[j + nc * j] = 1.0 / XWXptr[j + nc * j];
      mu1ptr[j] = XWXptr[j + nc * j] * mu0ptr[j];
      qbeta += dnorm(thetaptr[j], mu1ptr[j], pow(XWXptr[j + nc * j], 0.5), 1);
      thetaptr[j] = gamma1ptr[j];
    }
  }

  SEXP P2;
  PROTECT(P2 = eval_prior(getListElement(x, "prior"), theta2, rho));
  ++nProtected;
  double p2 = REAL(P2)[0];

  SEXP alpha;
  PROTECT(alpha = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(alpha)[0] = (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1);

  SEXP loglik;
  PROTECT(loglik = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(loglik)[0] = pibetaprop;

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 4));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, theta2);
  SET_VECTOR_ELT(rval, 1, alpha);
  SET_VECTOR_ELT(rval, 2, edf);
  SET_VECTOR_ELT(rval, 3, loglik);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 4));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("parameters"));
  SET_STRING_ELT(nrval, 1, mkChar("alpha"));
  SET_STRING_ELT(nrval, 2, mkChar("edf"));
  SET_STRING_ELT(nrval, 3, mkChar("loglik"));
        
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(nProtected);
  return rval;
}


/* Censored normal: left = 0, right = Inf. */
SEXP cnorm_loglik(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, 1));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  int *checkptr = INTEGER(check);

  double ll = 0.0;
  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ll += pnorm5((-1.0 * muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 1);
    } else {
      ll += dnorm((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1) - log(sigmaptr[i]);
    }
  }

  REAL(rval)[0] = ll;
  UNPROTECT(1);
  return rval;
}


SEXP cnorm_score_mu(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      rvalptr[i] = -1 * mills / sigmaptr[i];
    } else {
      rvalptr[i] = (yptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0);
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cnorm_score_sigma(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      rvalptr[i] = mills * muptr[i] / sigmaptr[i];
    } else {
      rvalptr[i] = (yptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0) * (yptr[i] - muptr[i]) - 1.0;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cnorm_hess_mu(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills, d1;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      d1 = -muptr[i] / pow(sigmaptr[i], 2.0);
      rvalptr[i] = -1 * (-d1 / sigmaptr[i] * mills - pow(mills, 2.0) / pow(sigmaptr[i], 2.0));
    } else {
      rvalptr[i] = 1 / pow(sigmaptr[i], 2.0);
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cnorm_hess_sigma(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills, d1, d2;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      d1 = -muptr[i] / pow(sigmaptr[i], 2.0);
      d2 = d1 * -muptr[i];
      rvalptr[i] = -1 * ((-muptr[i] / sigmaptr[i] - (0.0 - muptr[i]) * d2) * mills - pow(muptr[i], 2.0) / pow(sigmaptr[i], 2.0) * pow(mills, 2.0));
    } else {
      rvalptr[i] = 2 / pow(sigmaptr[i], 2.0) * pow(yptr[i] - muptr[i], 2.0);
    }
  }

  UNPROTECT(1);
  return rval;
}


/* Censored normal, left = 0, right - Inf, with power parameter. */
SEXP cnorm_power_loglik(SEXP y, SEXP mu, SEXP sigma, SEXP lambda, SEXP check)
{
  int i, n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *lambdaptr = REAL(lambda);
  int *checkptr = INTEGER(check);

  double ll = 0.0;
  double onediv = 0.0;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ll += pnorm5(0.0 , muptr[i], sigmaptr[i], 1, 1);
    } else {
      onediv = 1.0 / lambdaptr[i];
      ll += dnorm(pow(yptr[i], onediv), muptr[i], sigmaptr[i], 1) -
        log(lambdaptr[i]) + (onediv - 1.0) * log(yptr[i]);
    }
  }

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, 1));
  REAL(rval)[0] = ll;

  UNPROTECT(1);

  return rval;
}


SEXP cnorm_power_score_lambda(SEXP y, SEXP mu, SEXP sigma, SEXP lambda, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  double *lambdaptr = REAL(lambda);
  int *checkptr = INTEGER(check);

  double tmp, tmp2, tmp3;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      rvalptr[i] = 0.0;
    } else {
      tmp = exp(-log(lambdaptr[i]));
      tmp2 = log(yptr[i]);
      tmp3 = pow(yptr[i], tmp);
      rvalptr[i] = (tmp * tmp3 * tmp2 * (tmp3 - muptr[i])) / pow(sigmaptr[i], 2.0) - tmp * tmp2 - 1.0;
    }
  }

  UNPROTECT(1);
  return rval;
}


/* Fast quantile computation */
void swapd(double *a, double *b)      
{ 
  double temp;
  temp = *a;
  *a = *b;
  *b = temp;
} 
    
void quicksort_body(double *x, int up, int down)
{ 
  int start, end;
  start = up;              
  end = down;                
  while(up < down) {            
    while(x[down] >= x[up] && up < down)      
      down--;              
    if(up != down) {                    
      swapd(&x[up], &x[down]);   
      up++;    
    }
    while(x[up] <= x[down] && up < down)        
      up++;                 
    if(up != down) {                    
      swapd(&x[up], &x[down]);   
      down--;   
    } 
  }       
  if(start < up)   
    quicksort_body(x, start, up - 1); 
  if(end > down)  
    quicksort_body(x, down + 1, end);  
}

void quicksort(int n, double *x)
{ 
  quicksort_body(x, 0, n - 1);    
}


SEXP quick_quantiles(SEXP X, SEXP samples)
{
  int i, j, ii;
  int iter, nr, nc, nProtected = 0;
  SEXP out, TMP, q1, q2, q3, names;
    
  nr = nrows(X);
  nc = ncols(X);
  iter = nrows(samples);
    
  PROTECT(names = allocVector(STRSXP, 3));
  ++nProtected;
        
  PROTECT(out = allocVector(VECSXP, 3));
  ++nProtected;

  PROTECT(TMP = allocVector(REALSXP, iter));
  ++nProtected;
    
  PROTECT(q1 = allocVector(REALSXP, nr));
  ++nProtected;
    
  PROTECT(q2 = allocVector(REALSXP, nr));
  ++nProtected;
    
  PROTECT(q3 = allocVector(REALSXP, nr));
  ++nProtected;
    
  double np11 = iter * 0.025;
  double np12 = iter * 0.5;
  double np13 = iter * 0.975;
    
  int np1 = iter - np11;
  int np2 = iter - np12;
  int np3 = iter - np13;
    
  double *Xptr, *sptr, *tptr, *q1ptr,*q2ptr, *q3ptr;
  Xptr = REAL(X);
  sptr = REAL(samples);
  tptr = REAL(TMP);
  q1ptr = REAL(q1);
  q2ptr = REAL(q2);
  q3ptr = REAL(q3);

  double tmp = 0.0;
    
  for(i = 0; i < nr; i++) {
    for(ii = 0; ii < iter; ii++) {
      tmp = 0.0;
      for(j = 0; j < nc; j++) {
        tmp += Xptr[i + j * nr] * sptr[ii + j * iter];
      }
      tptr[ii] = tmp;
    }

    quicksort(iter, tptr);
              
    if((np11 - floor(np11)) == 0.0) {
      q1ptr[i] = (tptr[np1 - 1] + tptr[np1]) / 2.0;
    } else {
      q1ptr[i] = tptr[np1 - 1];
    }
    if((np12 - floor(np12)) == 0.0) { 
      q2ptr[i] = (tptr[np2 - 1] + tptr[np2]) / 2.0;
    } else {
      q2ptr[i] = tptr[np2 - 1];
    }
    if((np13 - floor(np13)) == 0.0) { 
      q3ptr[i] = (tptr[np3 - 1] + tptr[np3]) / 2.0;
    } else {
      q3ptr[i] = tptr[np3 - 1];
    }
  }
        
  SET_VECTOR_ELT(out, 0, q1);
  SET_VECTOR_ELT(out, 1, q2);
  SET_VECTOR_ELT(out, 2, q3);
    
  SET_STRING_ELT(names, 0, mkChar("lo"));
  SET_STRING_ELT(names, 1, mkChar("med"));
  SET_STRING_ELT(names, 2, mkChar("up"));
    
  setAttrib(out, R_NamesSymbol, names);
    
  UNPROTECT(nProtected);
  return out;
}


SEXP fitted_matrix(SEXP X, SEXP samples)
{
  int i, j, ii;    
  int nr = nrows(X);
  int nc = ncols(X);
  int iter = nrows(samples);
    
  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, nr, iter));
  double *outptr = REAL(out);
    
  double *Xptr = REAL(X);
  double *sptr = REAL(samples);
  double tmp = 0.0;
    
  for(i = 0; i < nr; i++) {
    for(ii = 0; ii < iter; ii++) {
      tmp = 0.0;
      for(j = 0; j < nc; ++j) {
        tmp += Xptr[i + j * nr] * sptr[ii + j * iter];
      }
      outptr[i + ii * nr] = tmp;
    }
  }

  UNPROTECT(1);
  return out;
}


/* Survival integrals. */
SEXP survint(SEXP X, SEXP eta, SEXP width, SEXP gamma, SEXP eta2, SEXP check)
{
  double *Xptr = REAL(X);
  double *etaptr = REAL(eta);
  double *eta2ptr = REAL(eta2);
  double *gammaptr = REAL(gamma);
  double *widthptr = REAL(width);

  int nProtected = 0;

  int nr = nrows(X);
  int nc = ncols(X);

  int tnr = nrows(eta);
  int tnc = ncols(eta);

  int ok = INTEGER(check)[0];

  SEXP grad;
  PROTECT(grad = allocVector(REALSXP, nc));
  ++nProtected;
  double *gradptr = REAL(grad);

  SEXP hess;
  PROTECT(hess = allocMatrix(REALSXP, nc, nc));
  ++nProtected;
  double *hessptr = REAL(hess);

  int i, ii, j, jj, k, forward;
  double sum = 0.0;
  double tmp = 0.0;

  for(j = 0; j < nc; j++) {
    for(jj = 0; jj <= j; jj++) {
      hessptr[j + jj * nc] = 0.0;
      hessptr[jj + j * nc] = 0.0;
    }
  }

  SEXP tmat;
  PROTECT(tmat = duplicate(hess));
  ++nProtected;
  double *tmatptr = REAL(tmat);

  for(j = 0; j < nc; j++) {
    gradptr[j] = 0.0;
    for(i = 0; i < tnr; i++) {
      sum = 0.0;
      for(k = 1; k < (tnc - 1); k++) {
        sum += Xptr[k + i * tnc + nr * j] * etaptr[i + k * tnr];
      }
      sum += 0.5 * (Xptr[i * tnc + nr * j] * etaptr[i] +
        Xptr[(tnc - 1) + i * tnc + nr * j] * etaptr[i + (tnc - 1) * tnr]);
      sum *= widthptr[i] * gammaptr[i];
      gradptr[j] += sum;

      if(j < 1) {
        forward = tnc * i;
        for(k = 0; k < tnc; k++) {
          for(jj = 0; jj < nc; jj++) {
            for(ii = 0; ii <= jj; ii++) {
              tmp = Xptr[k + forward + jj * nr] * Xptr[k + forward + ii * nr];
              if(ok < 1) {
                tmp *= eta2ptr[i + k * tnr];
              } else {
                tmp *= etaptr[i + k * tnr];
              }
              if(k == 0 || k == (tnc - 1)) {
                tmatptr[jj + ii * nc] += tmp * 0.5;
              } else {
                tmatptr[jj + ii * nc] += tmp;
              }
            }
          }
        }
        for(jj = 0; jj < nc; jj++) {
          for(ii = 0; ii <= jj; ii++) {
            tmp = tmatptr[jj + ii * nc] * widthptr[i];
            hessptr[jj + ii * nc] += tmp * gammaptr[i];
            hessptr[ii + jj * nc] = hessptr[jj + ii * nc];
            tmatptr[jj + ii * nc] = 0.0;
          }
        }
      }
    }
  }

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));
  ++nProtected;

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, grad);
  SET_VECTOR_ELT(rval, 1, hess);

  SET_STRING_ELT(nrval, 0, mkChar("grad"));
  SET_STRING_ELT(nrval, 1, mkChar("hess"));
        
  setAttrib(rval, R_NamesSymbol, nrval); 

  UNPROTECT(nProtected);
  return rval;
}


/* Survival integrals with index matrix. */
SEXP survint_index(SEXP X, SEXP eta, SEXP width, SEXP gamma, SEXP eta2, SEXP check, SEXP index)
{
  double *Xptr = REAL(X);
  double *etaptr = REAL(eta);
  double *eta2ptr = REAL(eta2);
  double *gammaptr = REAL(gamma);
  double *widthptr = REAL(width);
  int *indexptr = INTEGER(index);

  int nProtected = 0;

  int nr = nrows(X);
  int nc = ncols(X);
  int nc_index = ncols(index);
  int tnr = nrows(eta);
  int tnc = ncols(eta);

  int ok = INTEGER(check)[0];

  SEXP grad;
  PROTECT(grad = allocVector(REALSXP, nc));
  ++nProtected;
  double *gradptr = REAL(grad);

  SEXP hess;
  PROTECT(hess = allocMatrix(REALSXP, nc, nc));
  ++nProtected;
  double *hessptr = REAL(hess);

  int i, ii, j, jj, k, forward;
  double sum = 0.0;
  double tmp = 0.0;

  for(j = 0; j < nc; j++) {
    gradptr[j] = 0.0;
    for(jj = 0; jj <= j; jj++) {
      hessptr[j + jj * nc] = 0.0;
      hessptr[jj + j * nc] = 0.0;
    }
  }

  SEXP tmat;
  PROTECT(tmat = duplicate(hess));
  ++nProtected;
  double *tmatptr = REAL(tmat);

  for(i = 0; i < tnr; i++) {
    forward = tnc * i;
    for(j = 0; j < nc_index; j++) {
      jj = indexptr[i + j * tnr] - 1;
      if(jj < 0) continue;
      sum = 0.0;
      for(k = 0; k < tnc; k++) {
        ii = indexptr[i] - 1;
        while(ii <= jj) {
          tmp = Xptr[k + forward + jj * nr] * Xptr[k + forward + ii * nr];
          if(ok < 1) {
            tmp *= eta2ptr[i + k * tnr];
          } else {
            tmp *= etaptr[i + k * tnr];
          }
          if(k == 0 || k == (tnc - 1)) {
            tmatptr[jj + ii * nc] += tmp * 0.5;
          } else {
            tmatptr[jj + ii * nc] += tmp;
          }
          ii++;
        }
        if((k > 0) & (k < (tnc - 1))) {
          sum += Xptr[k + i * tnc + nr * jj] * etaptr[i + k * tnr];
        }
      }
      sum += 0.5 * (Xptr[i * tnc + nr * jj] * etaptr[i] +
        Xptr[(tnc - 1) + i * tnc + nr * jj] * etaptr[i + (tnc - 1) * tnr]);
      sum *= widthptr[i] * gammaptr[i];
      gradptr[jj] += sum;
      ii = indexptr[i] - 1;
      while(ii <= jj) {
        tmp = tmatptr[jj + ii * nc] * widthptr[i];
        hessptr[jj + ii * nc] += tmp * gammaptr[i];
        hessptr[ii + jj * nc] = hessptr[jj + ii * nc];
        tmatptr[jj + ii * nc] = 0.0;
        ii++;
      }
    }
  }

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));
  ++nProtected;

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, grad);
  SET_VECTOR_ELT(rval, 1, hess);

  SET_STRING_ELT(nrval, 0, mkChar("grad"));
  SET_STRING_ELT(nrval, 1, mkChar("hess"));
        
  setAttrib(rval, R_NamesSymbol, nrval); 

  UNPROTECT(nProtected);
  return rval;
}


/* Survival integrals with derivatives. */
SEXP dsurvint(SEXP X, SEXP eta, SEXP width, SEXP gamma, SEXP eta2, SEXP check,
  SEXP dX, SEXP deta, SEXP deta2)
{
  double *Xptr = REAL(X);
  double *etaptr = REAL(eta);
  double *eta2ptr = REAL(eta2);
  double *gammaptr = REAL(gamma);
  double *widthptr = REAL(width);

  double *dXptr = REAL(dX);
  double *detaptr = REAL(deta);
  double *deta2ptr = REAL(deta2);

  int nProtected = 0;

  int nr = nrows(X);
  int nc = ncols(X);

  int tnr = nrows(eta);
  int tnc = ncols(eta);

  int ok = INTEGER(check)[0];

  SEXP grad;
  PROTECT(grad = allocVector(REALSXP, nc));
  ++nProtected;
  double *gradptr = REAL(grad);

  SEXP hess;
  PROTECT(hess = allocMatrix(REALSXP, nc, nc));
  ++nProtected;
  double *hessptr = REAL(hess);

  int i, ii, j, jj, k, forward;
  double sum = 0.0;
  double tmp = 0.0;
  double dtmp = 0.0;

  for(j = 0; j < nc; j++) {
    for(jj = 0; jj <= j; jj++) {
      hessptr[j + jj * nc] = 0.0;
      hessptr[jj + j * nc] = 0.0;
    }
  }

  SEXP tmat;
  PROTECT(tmat = duplicate(hess));
  ++nProtected;
  double *tmatptr = REAL(tmat);

  for(j = 0; j < nc; j++) {
    gradptr[j] = 0.0;
    for(i = 0; i < tnr; i++) {
      sum = 0.0;
      for(k = 1; k < (tnc - 1); k++) {
        sum += Xptr[k + i * tnc + nr * j] * etaptr[i + k * tnr] + dXptr[k + i * tnc + nr * j] * detaptr[i + k * tnr];
      }
      sum += 0.5 * (Xptr[i * tnc + nr * j] * etaptr[i] +
        Xptr[(tnc - 1) + i * tnc + nr * j] * etaptr[i + (tnc - 1) * tnr]) +
          0.5 * (dXptr[i * tnc + nr * j] * detaptr[i] +
        dXptr[(tnc - 1) + i * tnc + nr * j] * detaptr[i + (tnc - 1) * tnr]);
      sum *= widthptr[i] * gammaptr[i];
      gradptr[j] += sum;

      if(j < 1) {
        forward = tnc * i;
        for(k = 0; k < tnc; k++) {
          for(jj = 0; jj < nc; jj++) {
            for(ii = 0; ii <= jj; ii++) {
              tmp = Xptr[k + forward + jj * nr] * Xptr[k + forward + ii * nr];
              dtmp = dXptr[k + forward + jj * nr] * dXptr[k + forward + ii * nr];
              if(ok < 1) {
                tmp *= eta2ptr[i + k * tnr];
                dtmp *= deta2ptr[i + k * tnr];
              } else {
                tmp *= etaptr[i + k * tnr];
                dtmp *= detaptr[i + k * tnr];
              }
              if(k == 0 || k == (tnc - 1)) {
                tmatptr[jj + ii * nc] += (tmp * 0.5 + dtmp * 0.5);
              } else {
                tmatptr[jj + ii * nc] += (tmp + dtmp);
              }
            }
          }
        }
        for(jj = 0; jj < nc; jj++) {
          for(ii = 0; ii <= jj; ii++) {
            tmp = tmatptr[jj + ii * nc] * widthptr[i];
            hessptr[jj + ii * nc] += tmp * gammaptr[i];
            hessptr[ii + jj * nc] = hessptr[jj + ii * nc];
            tmatptr[jj + ii * nc] = 0.0;
          }
        }
      }
    }
  }

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));
  ++nProtected;

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, grad);
  SET_VECTOR_ELT(rval, 1, hess);

  SET_STRING_ELT(nrval, 0, mkChar("grad"));
  SET_STRING_ELT(nrval, 1, mkChar("hess"));
        
  setAttrib(rval, R_NamesSymbol, nrval); 

  UNPROTECT(nProtected);
  return rval;
}


/* Survival integrals with index matrix. */
SEXP dsurvint_index(SEXP X, SEXP eta, SEXP width, SEXP gamma, SEXP eta2, SEXP check, SEXP index,
  SEXP dX, SEXP deta, SEXP deta2)
{
  double *Xptr = REAL(X);
  double *etaptr = REAL(eta);
  double *eta2ptr = REAL(eta2);
  double *gammaptr = REAL(gamma);
  double *widthptr = REAL(width);

  double *dXptr = REAL(dX);
  double *detaptr = REAL(deta);
  double *deta2ptr = REAL(deta2);

  int *indexptr = INTEGER(index);

  int nProtected = 0;

  int nr = nrows(X);
  int nc = ncols(X);
  int nc_index = ncols(index);
  int tnr = nrows(eta);
  int tnc = ncols(eta);

  int ok = INTEGER(check)[0];

  SEXP grad;
  PROTECT(grad = allocVector(REALSXP, nc));
  ++nProtected;
  double *gradptr = REAL(grad);

  SEXP hess;
  PROTECT(hess = allocMatrix(REALSXP, nc, nc));
  ++nProtected;
  double *hessptr = REAL(hess);

  int i, ii, j, jj, k, forward;
  double sum = 0.0;
  double tmp = 0.0;
  double dtmp = 0.0;

  for(j = 0; j < nc; j++) {
    gradptr[j] = 0.0;
    for(jj = 0; jj <= j; jj++) {
      hessptr[j + jj * nc] = 0.0;
      hessptr[jj + j * nc] = 0.0;
    }
  }

  SEXP tmat;
  PROTECT(tmat = duplicate(hess));
  ++nProtected;
  double *tmatptr = REAL(tmat);

  for(i = 0; i < tnr; i++) {
    forward = tnc * i;
    for(j = 0; j < nc_index; j++) {
      jj = indexptr[i + j * tnr] - 1;
      if(jj < 0) continue;
      sum = 0.0;
      for(k = 0; k < tnc; k++) {
        ii = indexptr[i] - 1;
        while(ii <= jj) {
          tmp = Xptr[k + forward + jj * nr] * Xptr[k + forward + ii * nr];
          dtmp = dXptr[k + forward + jj * nr] * Xptr[k + forward + ii * nr];
          if(ok < 1) {
            tmp *= eta2ptr[i + k * tnr];
            dtmp *= deta2ptr[i + k * tnr];
          } else {
            tmp *= etaptr[i + k * tnr];
            dtmp *= detaptr[i + k * tnr];
          }
          if(k == 0 || k == (tnc - 1)) {
            tmatptr[jj + ii * nc] += 0.5 * (tmp + dtmp);
          } else {
            tmatptr[jj + ii * nc] += (tmp + dtmp);
          }
          ii++;
        }
        if((k > 0) & (k < (tnc - 1))) {
          sum += Xptr[k + i * tnc + nr * jj] * etaptr[i + k * tnr] + dXptr[k + i * tnc + nr * jj] * detaptr[i + k * tnr];
        }
      }
      sum += 0.5 * (Xptr[i * tnc + nr * jj] * etaptr[i] +
        Xptr[(tnc - 1) + i * tnc + nr * jj] * etaptr[i + (tnc - 1) * tnr]) +
          0.5 * (dXptr[i * tnc + nr * jj] * detaptr[i] +
        dXptr[(tnc - 1) + i * tnc + nr * jj] * detaptr[i + (tnc - 1) * tnr]);
      sum *= widthptr[i] * gammaptr[i];
      gradptr[jj] += sum;
      ii = indexptr[i] - 1;
      while(ii <= jj) {
        tmp = tmatptr[jj + ii * nc] * widthptr[i];
        hessptr[jj + ii * nc] += tmp * gammaptr[i];
        hessptr[ii + jj * nc] = hessptr[jj + ii * nc];
        tmatptr[jj + ii * nc] = 0.0;
        ii++;
      }
    }
  }

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));
  ++nProtected;

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, grad);
  SET_VECTOR_ELT(rval, 1, hess);

  SET_STRING_ELT(nrval, 0, mkChar("grad"));
  SET_STRING_ELT(nrval, 1, mkChar("hess"));
        
  setAttrib(rval, R_NamesSymbol, nrval); 

  UNPROTECT(nProtected);
  return rval;
}


/* Extract the XT matrix. */
SEXP extract_XT(SEXP X, SEXP TNR, SEXP TNC)
{
  int nr = nrows(X);
  int nc = ncols(X);
  int tnr = INTEGER(TNR)[0];
  int tnc = INTEGER(TNC)[0];
  int i, j;

  double *Xptr = REAL(X);

  SEXP XT;
  PROTECT(XT = allocMatrix(REALSXP, tnr, nc));
  double *XTptr = REAL(XT);

  for(i = 0; i < tnr; i++) {
    for(j = 0; j < nc; j++) {
      XTptr[i + j * tnr] = Xptr[(tnc - 1) + i * tnc + nr * j];
    }
  }

  UNPROTECT(1);
  return XT;
}


/* Fast block diagonal crossproduct with weights. */
SEXP do_XWX(SEXP x, SEXP w, SEXP index)
{
  int nr = nrows(x);
  int nc = ncols(x);
  int nc_index = ncols(index);
  int i, j, k;

  double *xptr = REAL(x);
  double *wptr = REAL(w);
  int *iptr = INTEGER(index);

  SEXP rval;
  PROTECT(rval = allocMatrix(REALSXP, nc, nc));
  double *rvalptr = REAL(rval);

  for(j = 0; j < nc; j++) {
    for(k = 0; k <= j; k++) {
      rvalptr[j + k * nc] = 0.0;
      rvalptr[k + j * nc] = 0.0;
    }
  }

  for(j = 0; j < nc_index; j++) {
    for(k = 0; k < nc_index; k++) {
      for(i = 0; i < nr; i++) {
        if((iptr[i + j * nr] < 0) || (iptr[i + k * nr] < 0))
          continue;
        rvalptr[iptr[i + j * nr] - 1 + (iptr[i + k * nr] - 1) * nc] += xptr[i + (iptr[i + j * nr] - 1) * nr] * (1.0 / wptr[i]) * xptr[i + (iptr[i + k * nr] - 1) * nr];
      }
    }
  }

  UNPROTECT(1);
  return rval;
}


/* Fast computation of fitted values with index matrix. */
SEXP sparse_matrix_fit_fun(SEXP x, SEXP b, SEXP index)
{
  double *xptr = REAL(x);
  double *bptr = REAL(b);
  int *iptr = INTEGER(index);

  int nr = nrows(x);
  int nc = ncols(x);
  int nc_index = ncols(index);
  int i, j;

  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, nr));
  double *rvalptr = REAL(rval);

  for(i = 0; i < nr; i++) {
    rvalptr[i] = 0.0;
    for(j = 0; j < nc_index; j++) {
      if((iptr[i + j * nr] < 0) || (iptr[i + j * nr] > nc))
        continue;
      rvalptr[i] += xptr[i + (iptr[i + j * nr] - 1) * nr] * bptr[(iptr[i + j * nr] - 1)];
    }
  }

  UNPROTECT(1);
  return rval;
}


/* Fast matrix scaling. */
SEXP scale_matrix(SEXP x, SEXP center, SEXP scale)
{
  int nr = nrows(x); 
  int nc = ncols(x);
  int i, j;

  double *xptr = REAL(x);
  double *centerptr = REAL(center);
  double *scaleptr = REAL(scale);

  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      xptr[i + nr * j] = (xptr[i + nr * j] - centerptr[j]) / scaleptr[j];
    }
  }

  return x;
}


/* Log-likelihood multivariate normal distribution. */
SEXP mvn_logdens(SEXP y1, SEXP y2, SEXP mu1, SEXP mu2, SEXP sigma1, SEXP sigma2, SEXP rho)
{
  int i, n = length(y1);

  double *y1ptr = REAL(y1);
  double *y2ptr = REAL(y2);
  double *mu1ptr = REAL(mu1);
  double *mu2ptr = REAL(mu2);
  double *sigma1ptr = REAL(sigma1);
  double *sigma2ptr = REAL(sigma2);
  double *rhoptr = REAL(rho);

  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

  double log2pi = -1.83787706640935;

  for(i = 0; i < n; i++) {
    dptr[i] = log2pi - log(sigma1ptr[i]) - log(sigma2ptr[i]) - 0.5 * log(1.0 - pow(rhoptr[i], 2.0)) -
      1.0 / (2.0 * (1.0 - pow(rhoptr[i], 2.0))) * (pow((y1ptr[i] - mu1ptr[i]) / sigma1ptr[i], 2.0) -
      2.0 * rhoptr[i] * ((y1ptr[i] - mu1ptr[i]) * (y2ptr[i] - mu2ptr[i])) / (sigma1ptr[i] * sigma2ptr[i]) +
      pow((y2ptr[i] - mu2ptr[i]) / sigma2ptr[i], 2.0));
  }

  UNPROTECT(1);

  return d;
}

SEXP bivnorm_loglik(SEXP y1, SEXP y2, SEXP mu1, SEXP mu2, SEXP sigma1, SEXP sigma2, SEXP rho)
{
  int i, n = length(y1);

  double *y1ptr = REAL(y1);
  double *y2ptr = REAL(y2);
  double *mu1ptr = REAL(mu1);
  double *mu2ptr = REAL(mu2);
  double *sigma1ptr = REAL(sigma1);
  double *sigma2ptr = REAL(sigma2);
  double *rhoptr = REAL(rho);
  double sum = 0.0;
  double log2pi = -1.83787706640935;

  for(i = 0; i < n; i++) {
    sum += log2pi - log(sigma1ptr[i]) - log(sigma2ptr[i]) - 0.5 * log(1.0 - pow(rhoptr[i], 2.0)) -
      1.0 / (2.0 * (1.0 - pow(rhoptr[i], 2.0))) * (pow((y1ptr[i] - mu1ptr[i]) / sigma1ptr[i], 2.0) -
      2.0 * rhoptr[i] * ((y1ptr[i] - mu1ptr[i]) * (y2ptr[i] - mu2ptr[i])) / (sigma1ptr[i] * sigma2ptr[i]) +
      pow((y2ptr[i] - mu2ptr[i]) / sigma2ptr[i], 2.0));
  }

  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  REAL(ll)[0] = sum;

  UNPROTECT(1);

  return ll;
}


SEXP log_dmvnorm(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP MJ, SEXP SJ, SEXP RJ)
{
  int i, j, jj, l;
  int n = INTEGER(N)[0];
  int k = INTEGER(K)[0];
  int info;

  SEXP Sigma;
  PROTECT(Sigma = allocMatrix(REALSXP, k, k));
  double *Sigmaptr = REAL(Sigma);
  double *Pptr = REAL(PAR);

  int *MJptr = INTEGER(MJ);
  int *SJptr = INTEGER(SJ);
  int rj = INTEGER(RJ)[0];

  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);
  double *Yptr = REAL(Y);

  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

  double lpi = - 0.5 * k * log(2.0 * 3.14159265358979323846);
  double det = 0.0;
  double sum = 0.0;

  for(j = 0; j < k; j++) {
    for(i = 0; i < k; i++) {
      Sigmaptr[i + k * j] = 0.0;
    }
  }

  for(i = 0; i < n; i++) {
    l = 0;
    for(j = 0; j < k; j++) {
      Sigmaptr[j + k * j] = pow(Pptr[i + n * (SJptr[j] - 1)], 2.0);
      for(jj = j + 1; jj < k; jj++) {
        Sigmaptr[jj + k * j] = Pptr[i + n * (rj + l - 1)] * Pptr[i + n * (SJptr[jj] - 1)] * Pptr[i + n * (SJptr[j] - 1)];
        Sigmaptr[j + k * jj] = Sigmaptr[jj + k * j];
        l = l + 1;
      }
      ymuptr[j] = Yptr[i + n * j] - Pptr[i + n * (MJptr[j] - 1)];
    }

    F77_CALL(dpotrf)("Upper", &k, Sigmaptr, &k, &info);

    det = 0.0;
    for(j = 0; j < k; j++)
      det += log(Sigmaptr[j + k * j]);
    det = det * 2.0;

    F77_CALL(dpotri)("Upper", &k, Sigmaptr, &k, &info);

    sum = 0.0;
    for(j = 0; j < k; j++) {
      for(jj = j + 1; jj < k; jj++)
        Sigmaptr[jj + j * k] = Sigmaptr[j + jj * k];
      for(jj = 0; jj < k; jj++)
        sum += Sigmaptr[j + jj * k] * ymuptr[j] * ymuptr[jj];
    }

    dptr[i] = lpi - 0.5 * det - 0.5 * sum;
  }

  UNPROTECT(3);
  return d;
}


SEXP mu_score_mvnorm(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP MJ, SEXP SJ, SEXP RJ, SEXP KJ)
{
  int i, j, jj, l;
  int n = INTEGER(N)[0];
  int k = INTEGER(K)[0];
  int kj = INTEGER(KJ)[0];
  int info;

  SEXP Sigma;
  PROTECT(Sigma = allocMatrix(REALSXP, k, k));
  double *Sigmaptr = REAL(Sigma);
  double *Pptr = REAL(PAR);

  int *MJptr = INTEGER(MJ);
  int *SJptr = INTEGER(SJ);
  int rj = INTEGER(RJ)[0];

  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);
  double *Yptr = REAL(Y);

  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

  double sum = 0.0;

  for(j = 0; j < k; j++) {
    for(i = 0; i < k; i++) {
      Sigmaptr[i + k * j] = 0.0;
    }
  }

  for(i = 0; i < n; i++) {
    l = 0;
    for(j = 0; j < k; j++) {
      Sigmaptr[j + k * j] = pow(Pptr[i + n * (SJptr[j] - 1)], 2.0);
      for(jj = j + 1; jj < k; jj++) {
        Sigmaptr[jj + k * j] = Pptr[i + n * (rj + l - 1)] * Pptr[i + n * (SJptr[jj] - 1)] * Pptr[i + n * (SJptr[j] - 1)];
        Sigmaptr[j + k * jj] = Sigmaptr[jj + k * j];
        l = l + 1;
      }
      ymuptr[j] = Yptr[i + n * j] - Pptr[i + n * (MJptr[j] - 1)];
    }

    F77_CALL(dpotrf)("Upper", &k, Sigmaptr, &k, &info);
    F77_CALL(dpotri)("Upper", &k, Sigmaptr, &k, &info);

    sum = 0.0;
    /* Fill lower part of row kj */
    for(jj = 0; jj < kj; jj++)
      Sigmaptr[kj + jj * k] = Sigmaptr[jj + kj * k];
    /* sum over row kj */
    for(jj = 0; jj < k; jj++)
      sum += Sigmaptr[kj + jj * k] * ymuptr[jj];

    dptr[i] = sum;
  }

  UNPROTECT(3);
  return d;
}

SEXP sigma_score_mvnorm(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP MJ, SEXP SJ, SEXP RJ, SEXP KJ)
{
  int i, j, jj, l;
  int n = INTEGER(N)[0];
  int k = INTEGER(K)[0];
  int kj = INTEGER(KJ)[0];
  int info;

  SEXP Sigma;
  PROTECT(Sigma = allocMatrix(REALSXP, k, k));
  double *Sigmaptr = REAL(Sigma);
  double *Pptr = REAL(PAR);

  int *MJptr = INTEGER(MJ);
  int *SJptr = INTEGER(SJ);
  int rj = INTEGER(RJ)[0];

  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);
  double *Yptr = REAL(Y);

  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

  double sum = 0.0;

  /* Here Sigma is the correlation matrix */
  for(j = 0; j < k; j++) {
    for(i = 0; i < k; i++) {
      Sigmaptr[i + k * j] = 0.0;
    }
  }

  for(i = 0; i < n; i++) {
    l = 0;
    for(j = 0; j < k; j++) {
      Sigmaptr[j + k * j] = 1.0;
      for(jj = j + 1; jj < k; jj++) {
        Sigmaptr[jj + k * j] = Pptr[i + n * (rj + l - 1)];
        Sigmaptr[j + k * jj] = Sigmaptr[jj + k * j];
        l = l + 1;
      }
      ymuptr[j] = (Yptr[i + n * j] - Pptr[i + n * (MJptr[j] - 1)]) / Pptr[i + n * (SJptr[j] - 1)];
    }

    F77_CALL(dpotrf)("Upper", &k, Sigmaptr, &k, &info);
    F77_CALL(dpotri)("Upper", &k, Sigmaptr, &k, &info);


    sum = 0.0;
    /* Fill lower part of row kj */
    for(jj = 0; jj < kj; jj++)
      Sigmaptr[kj + jj * k] = Sigmaptr[jj + kj * k];
    for(jj = 0; jj < k; jj++)
      sum += Sigmaptr[kj + jj * k] * ymuptr[jj];

    dptr[i] = ymuptr[kj] * sum - 1.0;
  }

  UNPROTECT(3);
  return d;
}

SEXP rho_score_mvnorm(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP MJ, SEXP SJ, SEXP RJ, SEXP KJ, SEXP LJ)
{
  int i, j, jj, l;
  int n = INTEGER(N)[0];
  int k = INTEGER(K)[0];
  int kj = INTEGER(KJ)[0];
  int lj = INTEGER(LJ)[0];
  int info;

  SEXP Sigma;
  PROTECT(Sigma = allocMatrix(REALSXP, k, k));
  double *Sigmaptr = REAL(Sigma);
  double *Pptr = REAL(PAR);

  int *MJptr = INTEGER(MJ);
  int *SJptr = INTEGER(SJ);
  int rj = INTEGER(RJ)[0];

  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);
  double *Yptr = REAL(Y);

  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

  double sumkj = 0.0;
  double sumlj = 0.0;
  double mu = 0.0;
  double eta = 0.0;
  double mu2 = 0.0;
  double eta2 = 0.0;
  double base = 0.0;
  double deriv = 0.0;

  /* Here Sigma is the correlation matrix */
  for(j = 0; j < k; j++) {
    for(i = 0; i < k; i++) {
      Sigmaptr[i + k * j] = 0.0;
    }
  }

  for(i = 0; i < n; i++) {
    l = 0;
    for(j = 0; j < k; j++) {
      Sigmaptr[j + k * j] = 1.0;
      for(jj = j + 1; jj < k; jj++) {
        Sigmaptr[jj + k * j] = Pptr[i + n * (rj + l - 1)];
        Sigmaptr[j + k * jj] = Sigmaptr[jj + k * j];
        l = l + 1;
      }
      ymuptr[j] = (Yptr[i + n * j] - Pptr[i + n * (MJptr[j] - 1)]) / Pptr[i + n * (SJptr[j] - 1)];
    }

    /*compute deriv*/
    mu = Sigmaptr[kj + k * lj];
    mu2 = pow(mu, 2.0);
    base = 1.0 - mu2;
    eta = mu / pow(base, 0.5);
    eta2 = pow(eta, 2.0);
    base = 1.0 + eta2;
    deriv = 1.0 / pow(base, 1.5);

    /* Inversion */
    F77_CALL(dpotrf)("Upper", &k, Sigmaptr, &k, &info);
    F77_CALL(dpotri)("Upper", &k, Sigmaptr, &k, &info);

/*
    Fill lower part of rows kj and lj
    !! This was not working correctly.
    for(jj = 0; jj < kj; jj++)
      Sigmaptr[kj + jj * k] = Sigmaptr[jj + kj * k];
    for(jj = 0; jj < lj; jj++)
      Sigmaptr[lj + jj * k] = Sigmaptr[jj + jj * k];
*/
    /* Fill lower part of InvSigma */
    for(j = 0; j < k; j++) {
     for(jj = 0; jj < j; jj++) {
      Sigmaptr[j + jj * k] = Sigmaptr[jj + j * k];
     }
    }

    sumkj = 0.0;
    sumlj = 0.0;
    for(jj = 0; jj < k; jj++) {
      sumkj += Sigmaptr[kj + jj * k] * ymuptr[jj];
      sumlj += Sigmaptr[lj + jj * k] * ymuptr[jj];
    }

    dptr[i] = (sumkj * sumlj - Sigmaptr[kj + lj * k]) * 0.5 * deriv;
  }

  UNPROTECT(3);
  return d;
}



/* Boosting updater. */
SEXP boost_fit(SEXP x, SEXP y, SEXP nu, SEXP rho)
{
  int i, j, k;
  int nProtected = 0;
  int n          = length(y);
  int fixed      = LOGICAL(getListElement(x, "fixed"))[0];

  SEXP state       = PROTECT(duplicate(getListElement(x, "state")));
  ++nProtected;
  double *thetaptr = REAL(getListElement(state, "parameters"));
  int *penFun = INTEGER(getListElement(x, "penaltyFunction"));

  /* Reto: changes */
  /* This was the main change: Reto 2017-01-24 */
  int ntau2, S_ind;
  int nc = length(getListElement(state, "parameters"));
  if(fixed) {
    /* ntau2 can be anything but 0 */
    S_ind = 0;
    ntau2 = 1;
  } else {
    S_ind = getListElement_index(x, "S");
    ntau2 = length(PROTECT(VECTOR_ELT(x, S_ind))); ++nProtected;
    nc   -= ntau2;
  }
  SEXP tau2       = PROTECT(allocVector(REALSXP,ntau2));
  ++nProtected;
  double *tau2ptr = REAL(tau2);

  /* Create weighted matrix */
  int X_ind = getListElement_index(x, "X");
  int nr = nrows(PROTECT(VECTOR_ELT(x, X_ind)));
  ++nProtected;

  /* More pointers needed. */
  double *eptr = REAL(y);
  double *xweightsptr = REAL(getListElement(x, "weights"));
  double *xrresptr    = REAL(getListElement(x, "rres"));
  double *XWptr       = REAL(getListElement(x, "XW"));
  double *XWXptr      = REAL(getListElement(x, "XWX"));
  double *Xptr        = REAL(PROTECT(VECTOR_ELT(x, X_ind))); ++nProtected;
  /* Reto: changes, setting empty pointer to 0 */
  double *Sptr = 0;
  int    *idptr       = INTEGER(getListElement(getListElement(x, "binning"), "match.index"));
  int    *indptr      = INTEGER(getListElement(getListElement(x, "binning"), "sorted.index"));
  int    *orderptr    = INTEGER(getListElement(getListElement(x, "binning"), "order"));

  /* Handling fitted.values. */
  double *fitrptr = REAL(getListElement(x, "fit.reduced"));
  double *fitptr  = REAL(getListElement(state, "fitted.values"));

  /* Start. */
  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  j = 0;
  int jj;

  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      for(jj = 0; jj < nc; jj++) {
        XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
      }
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j]    = 0.0;
    }
    k = orderptr[i] - 1;

    xweightsptr[j] += 1.0;
    xrresptr[j] += eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
  }

  /* Compute X'WX. */
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, &nc, &nc, &nr, &one,
    XWptr, &nc, Xptr, &nr, &zero, XWXptr, &nc);

  /* Add penalty matrix and variance parameter. */
  if(!fixed) {
    for(jj = 0; jj < ntau2; jj++) {
      if(penFun[jj] > 0) {
        Sptr = REAL(get_S_mat(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj), getListElement(state, "parameters"),
          getListElement(state, "fixed.hyper"), rho));
      } else {
        Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      }
      tau2ptr[jj] = thetaptr[nc + jj];
      for(i = 0; i < nc; i++) {
        for(j = 0; j < nc; j++) {
          XWXptr[i + nc * j] += 1 / tau2ptr[jj] * Sptr[i + nc * j];
        }
      }
    }
  }

  /* Cholesky decompostion of XWX. */
  SEXP L = PROTECT(duplicate(getListElement(x, "XWX")));
  ++nProtected;
  double *Lptr = REAL(L);

  for(j = 0; j < nc; j++) { 	/* Zero the lower triangle. */
    for(i = j + 1; i < nc; i++) {
      Lptr[i + nc * j] = 0.0;
    }
  }

  int info;
  F77_CALL(dpotrf)("Upper", &nc, Lptr, &nc, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV = PROTECT(duplicate(L));  ++nProtected;
  double *PINVptr = REAL(PINV);

  F77_CALL(dpotri)("Upper", &nc, PINVptr, &nc, &info);

  SEXP PINVL = PROTECT(duplicate(PINV)); ++nProtected;
  double *PINVLptr = REAL(PINVL);
  F77_CALL(dpotrf)("Upper", &nc, PINVLptr, &nc, &info);

  for(j = 0; j < nc; j++) {
    for(i = j + 1; i < nc; i++) {
      PINVptr[i + j * nc] = PINVptr[j + i * nc];
    }
  }

  /* Compute mu. */
  SEXP mu0 = PROTECT(allocVector(REALSXP, nc)); ++nProtected;
  double *mu0ptr = REAL(mu0);

  SEXP mu1 = PROTECT(allocVector(REALSXP, nc)); ++nProtected;
  double *mu1ptr = REAL(mu1);

  int k1 = 1;
  char *transa2 = "T";
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nr, &one,
    Xptr, &nr, xrresptr, &nr, &zero, mu0ptr, &nc);
  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,
    PINVptr, &nc, mu0ptr, &nc, &zero, mu1ptr, &nc);

  for(j = 0; j < nc; j++) {
    thetaptr[j] = REAL(nu)[0] * mu1ptr[j];
  }

  int *iptr = INTEGER(getListElement(getListElement(x, "sparse.setup"), "matrix"));
  int nc_index = ncols(getListElement(getListElement(x, "sparse.setup"), "matrix"));

  for(i = 0; i < nr; i++) {
    fitrptr[i] = 0.0;
    for(j = 0; j < nc_index; j++) {
      if((iptr[i + j * nr] < 0) || (iptr[i + j * nr] > nc))
        continue;
      fitrptr[i] += Xptr[i + (iptr[i + j * nr] - 1) * nr] * thetaptr[(iptr[i + j * nr] - 1)];
    }
  }

  double *yptr = REAL(y);
  double rss = 0.0;
  for(i = 0; i < n; i++) {
    k = idptr[i] - 1;
    fitptr[i] = fitrptr[k];
    rss += pow(fitptr[i] - yptr[i], 2.0);
  }
  REAL(getListElement(state, "rss"))[0] = rss;

  UNPROTECT(nProtected);

  return state;
}


/* Generalized Pareto. */
SEXP gpareto_score_xi(SEXP y, SEXP xi, SEXP sigma)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *xiptr = REAL(xi);
  double *sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);

  double ys = 0.0;
  double xi1 = 0.0;
  double xi1ys = 0.0;

  for(i = 0; i < n; i++) {
    ys = yptr[i] / sigmaptr[i];
    xi1 = 1.0 / xiptr[i];
    xi1ys = 1.0 + xiptr[i] * ys;
    rvalptr[i] = -((xi1 + 1.0) * (xiptr[i] * ys/xi1ys) - xi1 * log(xi1ys));
  }

  UNPROTECT(1);
  return rval;
}

SEXP gpareto_score_sigma(SEXP y, SEXP xi, SEXP sigma)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *xiptr = REAL(xi);
  double *sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);

  double ys = 0.0;

  for(i = 0; i < n; i++) {
    ys = yptr[i] / sigmaptr[i];
    rvalptr[i] = (1.0/xiptr[i] + 1.0) * (xiptr[i] * ys /(1.0 + xiptr[i] * ys)) - 1.0;
  }

  UNPROTECT(1);
  return rval;
}

SEXP gpareto_hess_xi(SEXP y, SEXP xi, SEXP sigma)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *xiptr = REAL(xi);
  double *sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);

  double ys = 0.0;
  double xi1 = 0.0;
  double xi2 = 0.0;
  double xiys = 0.0;
  double xi1ys = 0.0;
  double xiysxi1ys = 0.0;
  double xi1xiysxi1ys = 0.0;

  for(i = 0; i < n; i++) {
    ys = yptr[i] / sigmaptr[i];
    xi1 = 1.0 / xiptr[i];
    xi2 = pow(xiptr[i], 2.0);
    xiys = xiptr[i] * ys;
    xi1ys = 1.0 + xiys;
    xiysxi1ys = xiys / xi1ys;
    xi1xiysxi1ys = xi1 * xiysxi1ys;
    rvalptr[i] = (xi1 + 1.0) * (xiysxi1ys - pow(xiys, 2.0)/pow(xi1ys, 2.0)) - xi1xiysxi1ys - ((xi1 - xiptr[i] * 2.0/xi2) * log(xi1ys) + xi1xiysxi1ys);

    if(ISNA(rvalptr[i]))
      rvalptr[i] = 1.490116e-08;
  }

  UNPROTECT(1);
  return rval;
}

SEXP gpareto_hess_sigma(SEXP y, SEXP xi, SEXP sigma)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *xiptr = REAL(xi);
  double *sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);

  double s1 = 0.0;
  double ys = 0.0;
  double xiys1 = 0.0;
  double xiys = 0.0;
  double xi1ys = 0.0;

  for(i = 0; i < n; i++) {
    s1 = 1.0 / sigmaptr[i];
    ys = yptr[i] / sigmaptr[i];
    xiys1 = xiptr[i] * yptr[i] * s1;
    xiys = xiptr[i] * ys;
    xi1ys = 1.0 + xiys;

    rvalptr[i] = (1.0/xiptr[i] + 1.0) * ((xiys1 - xiptr[i] * yptr[i] * sigmaptr[i] * 2.0 * pow(s1, 2.0))/xi1ys + xiys1 * xiys1/pow(xi1ys, 2.0));
    rvalptr[i] *= -1.0;

    if(ISNA(rvalptr[i]))
      rvalptr[i] = 1.490116e-08;
  }

  UNPROTECT(1);
  return rval;
}

