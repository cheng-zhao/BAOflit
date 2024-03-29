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
  * `Snl`:      array of the Sigma_nl parameters;
  * `c`:        array of the c parameters;
  * `args`:     structure for storing fitting arguments.
******************************************************************************/
static inline void xi_template(const double *Snl,
#ifdef PARA_MODEL
    const double *c,
#endif
    ARGS *args) {
  memset(args->xit[0], 0, sizeof(double) * args->ns * args->nxi);
  for (int k = 0; k < args->nxi; k++) {
    double Snl2 = Snl[k] * Snl[k];
    /* Compute the template power spectra. */
#ifdef PARA_MODEL
    for (size_t j = 0; j < args->nk; j++) {
      args->Pm[j] = (args->PBAO[j] * exp(-args->halfk2[j] * Snl2)
          + args->Pnw[j]) * (1 + c[k] * args->k[j]);    /* args->k stores k^2 */
    }
#else
    if (args->has_nwt[k]) {
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
#endif
    /* Integrate the template power spectra. */
    for (size_t i = 0; i < args->ns; i++) {
      size_t ioff = i * args->nk;
      for (size_t j = 0; j < args->nk; j++)
        args->xit[k][i] += args->Pm[j] * args->fac[ioff + j];
      /* Compute (st^2 * xit) for interpolation. */
      args->xit[k][i] *= args->st2[i];
    }
    /* Compute the second derivative for interpolating (st^2 * xit). */
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
      args->xim[i] = args->data[i]
          - B1 * B2 * args->xim[i] / (args->sm[i] * args->sm[i]);
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
  for (int k = 0; k < args->nxi; k++) {
    double *apoly = args->apoly + k * args->npoly;
    for (size_t j = args->idata[k]; j < args->edata[k]; j++) {
      double poly = 0;
      for (int i = 0; i < args->npoly; i++)
        poly += apoly[i] * args->basis[j * args->npoly + i];
      args->xim[j] -= poly;
    }
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
#ifdef PARA_MODEL
  if (args->Snltype != BAOFLIT_PARAM_FIX) {
    const double *Snl = par + 1 + args->num_B;
    const double *c = (args->ctype == BAOFLIT_PARAM_FIX) ?
        args->cval : par + args->cidx;
    xi_template(Snl, c, args);
  }
  else if (args->ctype != BAOFLIT_PARAM_FIX)
    xi_template(args->Snlval, par + 1 + args->num_B, args);
#else
  if (args->Snltype != BAOFLIT_PARAM_FIX)
    xi_template(par + 1 + args->num_B, args);
#endif
  xi_residual(par, args);
  if (args->npoly) least_square_fit(args);

  /* Compute the chi-squred value using forward substitution. */
  fwd_subst(args->Rcov, args->xim, args->nbin, args->xim);
  double chi2 = 0;
  for (size_t i = 0; i < args->nbin; i++) chi2 += args->xim[i] * args->xim[i];
  return chi2;
}

/******************************************************************************
Function `eval_model`:
  Compute the model 2PCFs given the fitting parameters.
Arguments:
  * `conf`:     structure for storing configuration parameters;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int eval_model(CONF *conf, ARGS *args) {
  /* Allocate memory for the best-fit 2PCFs. */
  double *s = malloc(sizeof(double) * args->ns * args->nxi);
  if (!s) {
    P_ERR("failed to allocate memory for the model 2PCFs\n");
    return BAOFLIT_ERR_MEMORY;
  }
  double *xi = malloc(sizeof(double) * args->ns * args->nxi);
  if (!xi) {
    P_ERR("failed to allocate memory for the model 2PCFs\n");
    free(s);
    return BAOFLIT_ERR_MEMORY;
  }
  size_t *idx = malloc(sizeof(size_t) * args->nxi);
  if (!idx) {
    P_ERR("failed to allocate memory for the model 2PCFs\n");
    free(s); free(xi);
    return BAOFLIT_ERR_MEMORY;
  }
  const char *bname[3] =
      {"model_mean.dat", "model_maxlike.dat", "model_map.dat"};
  size_t ntot = (size_t) args->npoly * (size_t) args->nxi;

  /* Evaluate the mean/best-fit/maximum-a-posteriori parameters in order. */
  for (int mid = 0; mid < 3; mid++) {
    double *pbest = args->pmodel + mid * args->npar;
    /* Perform the fit. */
    double chi2 = chi_squared(pbest, args);
    memcpy(args->amodel + ntot * mid, args->apoly, sizeof(double) * ntot);

    /* Recover the model 2PCFs. */
    size_t len = 0;
    for (int k = 0; k < args->nxi; k++) {
      idx[k] = len;
      double *sb = s + len;
      double *xib = xi + len;
      size_t ns = 0;
      /* Define the shifted abscissas. */
      for (size_t i = 0; i < args->ns; i++) {
        double ss = args->st[i] * pbest[0];
        if (ss < conf->smin) continue;
        if (ss >= conf->smax) break;
        sb[ns++] = ss;
      }
      if (ns == 0) {
        P_WRN("the model is not available for `%s%s'\n", conf->oroot, bname[k]);
        continue;
      }
      /* Interpolate the template 2PCFs. */
      cspline_eval_array(args->st, args->xit[k], args->xipp[k], args->ns,
          sb, xib, ns);
      /* Compute the model 2PCFs. */
      int i1 = args->idx_B[(k << 1)];
      int i2 = args->idx_B[(k << 1) + 1];
      double B1 = (i1 >= 0) ? pbest[i1 + 1] : BAOFLIT_DEFAULT_BIAS;
      double B2 = (i2 >= 0) ? pbest[i2 + 1] : BAOFLIT_DEFAULT_BIAS;
      for (size_t i = 0; i < ns; i++) {
        xib[i] *= B1 * B2 / (sb[i] * sb[i]);
        for (int j = 0; j < args->npoly; j++)
          xib[i] += args->apoly[j + k * args->npoly] * pow(sb[i], j - 2);
      }
      /* Shift back the abscissas. */
      for (size_t i = 0; i < ns; i++) sb[i] /= pbest[0];
      len += ns;
    }

    /* Save the best-fit model 2PCFs. */
    if (save_model(conf->oroot, bname[mid], s, xi, len, idx, args->nxi, chi2)) {
      free(s); free(xi); free(idx);
      return BAOFLIT_ERR_FILE;
    }
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
  if (args->Btype == BAOFLIT_PRIOR_GAUSS) {
    for (int i = 0; i < args->num_B; i++) {
      double d = (par[i + 1] - args->Bcen[i]) / args->Bsig[i];
      *lnlike -= 0.5 * d * d;
    }
  }
  if (args->Snltype == BAOFLIT_PRIOR_GAUSS) {
    for (int i = 0; i < args->nxi; i++) {
      double d = (par[i + args->num_B + 1] - args->Snlcen[i]) / args->Snlsig[i];
      *lnlike -= 0.5 * d * d;
    }
  }
#ifdef PARA_MODEL
  if (args->ctype == BAOFLIT_PRIOR_GAUSS) {
    for (int i = 0; i < args->nxi; i++) {
      double d = (par[args->cidx + i] - args->ccen[i]) / args->csig[i];
      *lnlike -= 0.5 * d * d;
    }
  }
#endif
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
  (void) nsample;
  (void) nlive;
  (void) npar;
  (void) physlive;
  (void) poster;
  (void) maxlike;
  (void) logZ;
  (void) INSlogZ;
  (void) logZerr;
  /* Record the mean/best-fit/MAP parameters and maximum log likelihood. */
  ARGS *args = (ARGS *) context;
  memcpy(args->pmodel, parinfo[0], sizeof(double) * (*npar));
  memcpy(args->pmodel + (*npar), parinfo[0] + 2 * (*npar),
      sizeof(double) * 2 * (*npar));
  args->maxlnlike = *maxlike;
}

/******************************************************************************
Function `run_multinest`:
  Perform the MultiNest fit.
Arguments:
  * `conf`:     the structure for storing configurations;
  * `fit`:      the structure for storing information for the fit.
******************************************************************************/
void run_multinest(CONF *conf, ARGS *args) {
  printf("Running BAO fit using MultiNest ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  /* The template 2PCFs can be pre-computed if Sigma_nl & c values are fixed. */
#ifdef PARA_MODEL
  if (args->Snltype == BAOFLIT_PARAM_FIX && args->ctype == BAOFLIT_PARAM_FIX)
    xi_template(conf->val_Snl, conf->val_c, args);
#else
  if (args->Snltype == BAOFLIT_PARAM_FIX) xi_template(conf->val_Snl, args);
#endif

#ifdef DEBUG
  for (int i = 0; i < args->npar * 3; i++) args->pmodel[i] = 1;
  /* Parameter orders: alpha, B_1, B_2, ..., B_n, Snl_1, Snl_2, ..., Snl_n. */
  args->pmodel[0] = 1;
  args->pmodel[1] = 1;
  args->pmodel[2] = 3;
  args->pmodel[3] = 1;
  args->pmodel[4] = 1;
  args->pmodel[5] = 5;
  args->pmodel[6] = 1;
  args->pmodel[7] = 1;
  args->pmodel[8] = 7;
  args->maxlnlike = -0.5 * chi_squared(args->pmodel, args);
#else
  /* Pre-process the filename to be accessed by fortran. */
  size_t flen = strlen(conf->oroot);
  for (size_t i = flen; i < BAOFLIT_MN_FNAME_LEN; i++) conf->oroot[i] = ' ';

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

  /* Run the MultiNest fit. */
  run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes,
        updInt, Ztol, conf->oroot, seed, pWrap, fb, resume, outfile, initMPI,
        logZero, maxiter, log_like, dumper, args);

  free(pWrap);
  conf->oroot[flen] = '\0';
#endif

  /* Evaluate and save the best-fit model. */
  if (eval_model(conf, args)) {
    P_WRN("failed to evaluate or save the best-fit model\n");
    return;
  }

  if (conf->verbose) {
    const char *ptype[3] =
        {"Mean values of the", "Maximum-likelihood", "Maximum-a-posteriori"};
    for (int mid = 0; mid < 3; mid++) {
      double *pbest = args->pmodel + mid * args->npar;
      printf("  %s parameters:\n    alpha: " OFMT_DBL "\n    bias:",
          ptype[mid], pbest[0]);

      for (int i = 0; i < args->num_B; i++)
        printf(" " OFMT_DBL, pbest[i + 1]);

      if (args->Snltype != BAOFLIT_PARAM_FIX) {
        printf("\n    Sigma_nl:");
        for (int i = 0; i < args->nxi; i++)
          printf(" " OFMT_DBL, pbest[i + 1 + args->num_B]);
      }

#ifdef PARA_MODEL
      if (args->ctype != BAOFLIT_PARAM_FIX) {
        printf("\n    c:");
        for (int i = 0; i < args->nxi; i++)
          printf(" " OFMT_DBL, pbest[args->cidx + i]);
      }
#endif

      if (args->npoly) {
        printf("\n    Nuisance parameters:");
        size_t ntot = (size_t) args->npoly * (size_t) args->nxi;
        for (size_t i = 0; i < ntot; i++)
          printf(" " OFMT_DBL, args->amodel[ntot * mid + i]);
      }
      printf("\n");
    }
    printf("  Minimum chi-squared: " OFMT_DBL "\n", -2 * args->maxlnlike);
  }

  printf(FMT_DONE);
}
