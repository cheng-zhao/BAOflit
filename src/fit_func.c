/*******************************************************************************
* fit_func.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "fit_func.h"
#include "multinest.h"
#include "save_res.h"
#include "cspline.h"
#include "linalg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*============================================================================*\
              Functions for model evaluation and parameter fitting
\*============================================================================*/

/******************************************************************************
Function `xi_template`:
  Compute the template 2-point correlation function.
Arguments:
  * `par`:      array of the Sigma_nl parameters;
  * `args`:     structure for storing fitting arguments.
******************************************************************************/
static inline void xi_template(const double *par, ARGS *args) {
  memset(args->xit[0], 0, sizeof(double) * args->ns * args->nxi);
  for (int k = 0; k < args->nxi; k++) {
    double Snl2 = par[k] * par[k];
    /* Compute the template power spectra. */
    if (args->Pnwt[k]) {
      for (size_t j = 0; j < args->nk; j++) {
        args->Pm[j] = (args->PBAO[j] * exp(-args->halfk2[j] * Snl2)
            + args->Pnw[j]) * args->Pnwt[k][j];
      }
    }
    else {
      for (size_t j = 0; j < args->nk; j++) {
        args->Pm[j] = args->PBAO[j] * exp(-args->halfk2[j] * Snl2)
            + args->Pnw[j];
      }
    }
    /* Integrate the template power spectra. */
    for (size_t i = 0; i < args->ns; i++) {
      size_t ioff = i * args->nk;
      for (size_t j = 0; j < args->nk; j++)
        args->xit[k][i] += args->Pm[j] * args->fac[ioff + j];
    }
    /* Compute the second derivative for interpolation. */
    cspline_ypp(args->st, args->xit[k], args->ns, args->xipp[k]);
  }
}

/******************************************************************************
Function `xi_residual`:
  Compute the difference between the data vector and the model 2PCFs.
Arguments:
  * `par`:      array of the free parameters;
  * `args`:     structure for storing fitting arguments.
******************************************************************************/
static inline void xi_residual(const double *par, ARGS *args) {
  for (int k = 0; k < args->nxi; k++) {
    /* Shift separation bins by alpha. */
    for (size_t i = args->idata[k]; i < args->edata[k]; i++)
      args->sm[i] = args->s[i] * par[0];
    /* Interpolate the template 2PCFs. */
    cspline_eval_array(args->st, args->xit[k], args->xipp[k], args->ns,
        args->sm + args->idata[k], args->xim + args->idata[k],
        args->edata[k] - args->idata[k]);

    /* Compute the residual */
    int i1 = args->idx_B[(k << 1)];
    int i2 = args->idx_B[(k << 1) + 1];
    double B1 = (i1 >= 0) ? par[i1 + 1] : BAOFLIT_DEFAULT_BIAS;
    double B2 = (i2 >= 0) ? par[i2 + 1] : BAOFLIT_DEFAULT_BIAS;
    for (size_t i = args->idata[k]; i < args->edata[k]; i++)
      args->xim[i] = args->data[i] - B1 * B2 * args->xim[i];
  }
}

/******************************************************************************
Function `least_square_fit`:
  Compute the nuisance parameters using the least-squared method.
Arguments:
  * `args`:     structure for storing fitting arguments.
******************************************************************************/
static inline void least_square_fit(ARGS *args) {
   size_t ntot = (size_t) args->npoly * (size_t) args->nxi;
  /* Compute the right-hand-side vector for least-squared fitting. */
  memset(args->apoly, 0, sizeof(double) * ntot);
  for (size_t i = 0; i < ntot; i++) {
    for (size_t j = 0; j < args->nbin; j++)
      args->apoly[i] += args->LS_Z[i * args->nbin + j] * args->xim[j];
  }
  /* Solve the linear equation by backward substitution. */
  bwd_subst(args->LS_U, args->apoly, ntot, args->apoly);

  /* Compute the residual with the nuisance parameters. */
  for (size_t j = 0; j < args->nbin; j++) {
    double poly = 0;
    for (size_t i = 0; i < ntot; i++)
      poly += args->apoly[i] * args->basis[j * args->npoly + i];
    args->xim[j] -= poly;
  }
}


/*============================================================================*\
                      Functions for likelihood evaluation
\*============================================================================*/

/******************************************************************************
Function `chi_squared`:
  Compute the chi-squared of the fit.
Arguments:
  * `par`:      array of free parameters;
  * `args`:     structure for storing fitting arguments.
Return:
  Address of the structure for 2PCF evaluation.
******************************************************************************/
static double chi_squared(const double *par, ARGS *args) {
  if (args->fit_Snl) xi_template(par + 1 + args->num_B, args);
  xi_residual(par, args);
  if (args->npoly) least_square_fit(args);

  /* Compute the chi-squred value using forward substitution. */
  fwd_subst(args->Rcov, args->xim, args->nbin, args->xim);
  double chi2 = 0;
  for (size_t i = 0; i < args->nbin; i++) chi2 += args->xim[i] * args->xim[i];
  return chi2;
}

/******************************************************************************
Function `best_fit`:
  Compute the best-fit multipoles given a set of parameters.
Arguments:
  * `conf`:     structure for storing configuration parameters;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int best_fit(const CONF *conf, ARGS *args) {
  if (args->fit_Snl) xi_template(args->pbest + 1 + args->num_B, args);
  xi_residual(args->pbest, args);
  if (args->npoly) least_square_fit(args);

  /* Allocate memory for the best-fit 2PCFs. */
  double *s = malloc(sizeof(double) * args->ns * args->nxi);
  if (!s) {
    P_ERR("failed to allocate memory for the best-fit 2PCFs\n");
    return BAOFLIT_ERR_MEMORY;
  }
  double *xi = malloc(sizeof(double) * args->ns * args->nxi);
  if (!xi) {
    P_ERR("failed to allocate memory for the best-fit 2PCFs\n");
    free(s);
    return BAOFLIT_ERR_MEMORY;
  }
  size_t *idx = malloc(sizeof(size_t) * args->nxi);
  if (!idx) {
    P_ERR("failed to allocate memory for the best-fit 2PCFs\n");
  }

  /* Recover the model 2PCFs. */
  size_t len = 0;
  for (int k = 0; k < args->nxi; k++) {
    idx[k] = len;
    double *sb = s + len;
    double *xib = xi + len;
    size_t ns = 0;
    /* Define the shifted abscissas. */
    for (size_t i = 0; i < args->ns; i++) {
      double ss = args->st[i] * args->pbest[0];
      if (ss < conf->fitmin[k]) continue;
      if (ss > conf->fitmax[k]) break;
      sb[ns++] = ss;
    }
    /* Interpolate the template 2PCFs. */
    cspline_eval_array(args->st, args->xit[k], args->xipp[k], args->ns,
        sb, xib, ns);
    /* Compute the model 2PCFs. */
    int i1 = args->idx_B[(k << 1)];
    int i2 = args->idx_B[(k << 1) + 1];
    double B1 = (i1 >= 0) ? args->pbest[i1 + 1] : BAOFLIT_DEFAULT_BIAS;
    double B2 = (i2 >= 0) ? args->pbest[i2 + 1] : BAOFLIT_DEFAULT_BIAS;
    for (size_t i = 0; i < ns; i++) {
      xib[i] *= B1 * B2;
      for (int j = 0; j < args->npoly; j++)
        xib[i] += args->apoly[j + k * args->npoly] * pow(sb[i], j - 2);
    }
    len += ns;
  }

  /* Save the best-fit model 2PCFs. */
  if (save_table(conf->fbest, s, xi, len, idx, args->nxi)) {
    free(s); free(xi); free(idx);
    return BAOFLIT_ERR_FILE;
  }

  free(s);
  free(xi);
  free(idx);
  return 0;
}

/******************************************************************************
Function `log_like`:
  Compute the log likelihood based on the chi-squared.
Arguments:
  * `par`:      array of input parameters;
  * `ndim`:     number of free parameters;
  * `npar`:     number of free plus derived parameters;
  * `lnlike`:   the evaluated log likelihood;
  * `context`:  arguments for computing the log likelihood.
******************************************************************************/
static void log_like(double *par, int *ndim, int *npar, double *lnlike,
    void *context) {
  (void) npar;
  ARGS *args = (ARGS *) context;

  /* Rescale the parameters according to the prior ranges. */
  for (int i = 0; i < *ndim; i++)
    par[i] = par[i] * (args->pmax[i] - args->pmin[i]) + args->pmin[i];

  /* Compute the log likelihood based on the chi-squared. */
  *lnlike = -0.5 * chi_squared(par, args);

  /* Add non-flat priors. */
  if (args->Btype == BIAS_PRIOR_GAUSS) {
    for (int i = 0; i < args->num_B; i++) {
      double d = (par[i + 1] - args->pcen[i]) / args->psig[i];
      *lnlike -= 0.5 * d * d;
    }
  }
}

/******************************************************************************
Function `dumper`:
  A function that is caled every updInt * 10 iterations.
  See `example_eggbox_C/eggbox.c` of the MultiNest library.
Arguments:
  * `nsample`:  number of samples in the posteior distribution;
  * `nlive`:    number of live points;
  * `npar`:     number of free plus derived parameters;
  * `physlive`: the last set of live points with their log likelihood values;
  * `poster`:   posterior distribution with nsample points;
  * `parinfo`:  information of the parameters;
  * `maxlike`:  maximum log likelihood value;
  * `logZ`:     log evidence value from the default (non-INS) mode;
  * `INSlogZ`:  log evidence value from the INS mode;
  * `logZerr`:  error on log evidence value;
  * `context`:  void pointer for any additional information.
******************************************************************************/
static void dumper(int *nsample, int *nlive, int *npar, double **physlive,
    double **poster, double **parinfo, double *maxlike, double *logZ,
    double *INSlogZ, double *logZerr, void *context) {
  (void) nlive;
  (void) npar;
  (void) physlive;
  (void) poster;
  (void) maxlike;
  (void) logZ;
  (void) INSlogZ;
  (void) logZerr;
  /* Record the best-fit parameters and maximum log likelihood. */
  ARGS *args = (ARGS *) context;
  memcpy(args->pbest, parinfo[0] + 3 * (*npar), sizeof(double) * (*npar));
  args->maxlnlike = *maxlike;
}

/******************************************************************************
Function `run_multinest`:
  Perform the MultiNest fit.
Arguments:
  * `conf`:     the structure for storing configurations;
  * `fit`:      the structure for storing information for the fit.
******************************************************************************/
void run_multinest(const CONF *conf, ARGS *args) {
  printf("Running BAO fit using MultiNest ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  /* Initialise MultiNest parameters */
  int IS = BAOFLIT_MN_IS;
  int mmodal = BAOFLIT_MN_MMODAL;
  int ceff = BAOFLIT_MN_CEFF;
  int nlive = conf->nlive;
  double efr = BAOFLIT_MN_EFR;
  double tol = conf->tol;
  int ndims = args->npar;
  int nPar = args->npar;
  int nClsPar = args->npar;
  int updInt = BAOFLIT_MN_UPD;
  double Ztol = BAOFLIT_MN_ZTOL;
  int maxModes = BAOFLIT_MN_MAXMODE;
  int seed = BAOFLIT_MN_SEED;
  int fb = BAOFLIT_MN_STDOUT;
  int resume = conf->resume ? 1 : 0;
  int outfile = 1;
  int initMPI = BAOFLIT_MN_INITMPI;
  double logZero = BAOFLIT_MN_LOGZERO;
  int maxiter = BAOFLIT_MN_MAXIT;
  int *pWrap = malloc(sizeof(int) * args->npar);
  for (int i = 0; i < args->npar; i++) pWrap[i] = BAOFLIT_MN_PWRAP;

  /* The template 2PCFs can be pre-computed if Sigma_nl values are fixed. */
  if (!args->fit_Snl) xi_template(conf->val_Snl, args);

  /* Run the MultiNest fit. */
  run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes,
        updInt, Ztol, conf->oroot, seed, pWrap, fb, resume, outfile, initMPI,
        logZero, maxiter, log_like, dumper, args);

  /* Evaluate and save the best-fit model. */
  if (best_fit(conf, args)) {
    P_WRN("failed to evaluate or save the best-fit model\n");
    return;
  }

  if (conf->verbose) {
    printf("  Best-fit parameters:\n    alpha: " OFMT_DBL "\n    bias:",
        args->pbest[0]);
    for (int i = 0; i < args->num_B; i++)
      printf(" " OFMT_DBL, args->pbest[i + 1]);

    if (args->fit_Snl) {
      printf("\n    Sigma_nl:");
      for (int i = 0; i < args->nxi; i++)
        printf(" " OFMT_DBL, args->pbest[1 + args->num_B + i]);
    }
    if (args->npoly) {
      printf("\n    Nuisance parameters:");
      size_t ntot = (size_t) args->npoly * (size_t) args->nxi;
      for (size_t i = 0; i < ntot; i++) printf(" " OFMT_DBL, args->apoly[i]);
    }
    printf("\n  Minimum chi-squared: " OFMT_DBL "\n", -2 * args->maxlnlike);
  }

  printf(FMT_DONE);
}
