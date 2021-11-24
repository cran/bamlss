#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define USE_FC_LEN_T

SEXP bamlss_glogis_score(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP bamlss_glogis_hesse(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP bamlss_glogis_density(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP bamlss_glogis_loglik(SEXP, SEXP, SEXP, SEXP);
SEXP bamlss_glogis_distr(SEXP, SEXP, SEXP, SEXP);
SEXP bamlss_glogis_quantile(SEXP, SEXP, SEXP, SEXP);
SEXP block_inverse(SEXP, SEXP, SEXP);
SEXP bivnorm_loglik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP boost_fit(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP cnorm_hess_mu(SEXP, SEXP, SEXP, SEXP);
SEXP cnorm_hess_sigma(SEXP, SEXP, SEXP, SEXP);
SEXP cnorm_loglik(SEXP, SEXP, SEXP, SEXP);
SEXP cnorm_power_loglik(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP cnorm_score_mu(SEXP, SEXP, SEXP, SEXP);
SEXP cnorm_score_sigma(SEXP, SEXP, SEXP, SEXP);
SEXP cpos(SEXP, SEXP);
SEXP do_XWX(SEXP, SEXP, SEXP);
SEXP dsurvint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dsurvint_index(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP extract_XT(SEXP, SEXP, SEXP);
SEXP fitted_matrix(SEXP, SEXP);
SEXP gmcmc_iwls(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP gmcmc_iwls_gp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP gmcmc_iwls_gp_diag_lasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP gpareto_hess_sigma(SEXP, SEXP, SEXP);
SEXP gpareto_hess_xi(SEXP, SEXP, SEXP);
SEXP gpareto_score_sigma(SEXP, SEXP, SEXP);
SEXP gpareto_score_xi(SEXP, SEXP, SEXP);
SEXP hatmat_trace(SEXP, SEXP);
SEXP hatmat_sumdiag(SEXP);
SEXP log_dmvnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP mu_score_mvnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP nnet_fitfun(SEXP, SEXP, SEXP);
SEXP process_derivs(SEXP, SEXP);
SEXP quick_quantiles(SEXP, SEXP);
SEXP rho_score_mvnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP scale_matrix(SEXP, SEXP, SEXP);
SEXP sigma_score_mvnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP sparse_matrix_fit_fun(SEXP, SEXP, SEXP);
SEXP sum_diag(SEXP, SEXP);
SEXP sum_diag2(SEXP, SEXP);
SEXP survint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP survint_index(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xbin_fun(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP log_dmvnormAR1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP mu_score_mvnormAR1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP sigma_score_mvnormAR1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rho_score_mvnormAR1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP ztnbinom_score_mu(SEXP, SEXP, SEXP);
SEXP ztnbinom_score_theta(SEXP, SEXP, SEXP);

SEXP xcenter(SEXP);

SEXP logNN_dens(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP logNN_score_mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP logNN_score_sigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP logNN_score_lambda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
  {"logNN_dens", (DL_FUNC) &logNN_dens, 6},
  {"logNN_score_mu", (DL_FUNC) &logNN_score_mu, 6},
  {"logNN_score_sigma", (DL_FUNC) &logNN_score_sigma, 6},
  {"logNN_score_lambda", (DL_FUNC) &logNN_score_lambda, 6},
  {"xcenter", (DL_FUNC) &xcenter, 1},
  {"bamlss_glogis_score", (DL_FUNC) &bamlss_glogis_score, 5},
  {"bamlss_glogis_hesse", (DL_FUNC) &bamlss_glogis_hesse, 5},
  {"bamlss_glogis_density", (DL_FUNC) &bamlss_glogis_density, 5},
  {"bamlss_glogis_loglik", (DL_FUNC) &bamlss_glogis_loglik, 4},
  {"bamlss_glogis_distr", (DL_FUNC) &bamlss_glogis_distr, 4},
  {"bamlss_glogis_quantile", (DL_FUNC) &bamlss_glogis_quantile, 4},
  {"block_inverse", (DL_FUNC) &block_inverse, 3},
  {"bivnorm_loglik", (DL_FUNC) &bivnorm_loglik, 7},
  {"boost_fit", (DL_FUNC) &boost_fit, 5},
  {"cnorm_hess_mu", (DL_FUNC) &cnorm_hess_mu, 4},
  {"cnorm_hess_sigma", (DL_FUNC) &cnorm_hess_sigma, 4},
  {"cnorm_loglik", (DL_FUNC) &cnorm_loglik, 4},
  {"cnorm_power_loglik", (DL_FUNC) &cnorm_power_loglik, 5},
  {"cnorm_score_mu", (DL_FUNC) &cnorm_score_mu, 4},
  {"cnorm_score_sigma", (DL_FUNC) &cnorm_score_sigma, 4},
  {"cpos", (DL_FUNC) &cpos, 2},
  {"do_XWX", (DL_FUNC) &do_XWX, 3},
  {"dsurvint", (DL_FUNC) &dsurvint, 9},
  {"dsurvint_index", (DL_FUNC) &dsurvint_index, 10},
  {"extract_XT", (DL_FUNC) &extract_XT, 3},
  {"fitted_matrix", (DL_FUNC) &fitted_matrix, 2},
  {"gmcmc_iwls", (DL_FUNC) &gmcmc_iwls, 11},
  {"gmcmc_iwls_gp", (DL_FUNC) &gmcmc_iwls_gp, 13},
  {"gmcmc_iwls_gp_diag_lasso", (DL_FUNC) &gmcmc_iwls_gp_diag_lasso, 11},
  {"gpareto_hess_sigma", (DL_FUNC) &gpareto_hess_sigma, 3},
  {"gpareto_hess_xi", (DL_FUNC) &gpareto_hess_xi, 3},
  {"gpareto_score_sigma", (DL_FUNC) &gpareto_score_sigma, 3},
  {"gpareto_score_xi", (DL_FUNC) &gpareto_score_xi, 3},
  {"hatmat_trace", (DL_FUNC) &hatmat_trace, 2},
  {"hatmat_sumdiag", (DL_FUNC) &hatmat_sumdiag, 1},
  {"log_dmvnorm", (DL_FUNC) &log_dmvnorm, 7},
  {"mu_score_mvnorm", (DL_FUNC) &mu_score_mvnorm, 8},
  {"nnet_fitfun", (DL_FUNC) &nnet_fitfun, 3},
  {"process_derivs", (DL_FUNC) &process_derivs, 2},
  {"quick_quantiles", (DL_FUNC) &quick_quantiles, 2},
  {"rho_score_mvnorm", (DL_FUNC) &rho_score_mvnorm, 9},
  {"scale_matrix", (DL_FUNC) &scale_matrix, 3},
  {"sigma_score_mvnorm", (DL_FUNC) &sigma_score_mvnorm, 8},
  {"sparse_matrix_fit_fun", (DL_FUNC) &sparse_matrix_fit_fun, 3},
  {"sum_diag", (DL_FUNC) &sum_diag, 2},
  {"sum_diag2", (DL_FUNC) &sum_diag2, 2},
  {"survint", (DL_FUNC) &survint, 6},
  {"survint_index", (DL_FUNC) &survint_index, 7},
  {"xbin_fun", (DL_FUNC) &xbin_fun, 6},
  {"log_dmvnormAR1", (DL_FUNC) &log_dmvnormAR1, 7},
  {"mu_score_mvnormAR1", (DL_FUNC) &mu_score_mvnormAR1, 8},
  {"sigma_score_mvnormAR1", (DL_FUNC) &sigma_score_mvnormAR1, 8},
  {"rho_score_mvnormAR1", (DL_FUNC) &rho_score_mvnormAR1, 7},
  {"ztnbinom_score_mu", (DL_FUNC) &ztnbinom_score_theta, 3},
  {"ztnbinom_score_theta", (DL_FUNC) &ztnbinom_score_theta, 3},
  {NULL, NULL, 0}
};

void R_init_sourcetools(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

