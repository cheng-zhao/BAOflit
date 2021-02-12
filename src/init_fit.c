/*******************************************************************************
* init_fit.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "read_data.h"
#include "linalg.h"
#include "legauss.h"
#include "cspline.h"
#include "proc_pk.h"
#include "save_res.h"
#include "fit_args.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/*============================================================================*\
                Functions for intialising the fitting arguments
\*============================================================================*/

/******************************************************************************
Function `init_args`:
  Initialise the structure for arguments of the fit.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for arguments of the fit.
******************************************************************************/
static ARGS *init_args(const CONF *conf) {
  ARGS *args = calloc(1, sizeof *args);

  args->pmin = args->pmax = args->pbest = NULL;
  args->data = args->s = args->Rcov = NULL;
  args->idata = args->edata = NULL;
  args->k = args->fac = args->halfk2 = args->PBAO = args->Pnw = args->Pm = NULL;
  args->Pnwt = args->xit = args->xipp = NULL;
  args->st = args->sm = args->xim = NULL;
  args->apoly = args->basis = args->LS_U = args->LS_Z = NULL;
  args->idx_B = NULL;

  args->num_B = conf->num_B;
  args->nxi = conf->ninput;
  args->ns = (size_t) conf->ns;
  args->npoly = conf->npoly;
  args->pcen = conf->pcen_B;
  args->psig = conf->psig_B;
  args->fit_Snl = (conf->val_Snl) ? false : true;
  args->maxlnlike = -DBL_MAX;

  switch (conf->Btype) {
    case BIAS_PRIOR_FLAT:
    case BIAS_PRIOR_GAUSS:
      args->Btype = conf->Btype;
      break;
    default:
      P_ERR("invalid bias prior type: %d\n", conf->Btype);
      free(args);
      return NULL;
  }
  switch (conf->pkint) {
    case PK_INT_TRAPZ:
    case PK_INT_LEGAUSS:
      args->pkint = conf->pkint;
      break;
    default:
      P_ERR("invalid power spectra integration method: %d\n", conf->pkint);
      free(args);
      return NULL;
  }

  args->idata = calloc(conf->ninput, sizeof(size_t));
  args->edata = calloc(conf->ninput, sizeof(size_t));
  if (!args->idata || !args->edata) {
    P_ERR("failed to allocate memory for 2PCF positions in the data vector\n");
    args_destroy(args);
    return NULL;
  }

  return args;
}

/******************************************************************************
Function `init_param`:
  Initialise the fitting parameters.
Arguments:
  * `conf`:     structure for storing configurations;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int init_param(const CONF *conf, ARGS *args) {
  /* Compute the number of free parameters. */
  int npar = 1 + conf->num_B;
  if (args->fit_Snl) npar += conf->ninput;
  args->npar = npar;

  /* Allocate memory. */
  args->pmin = malloc(sizeof(double) * npar);
  args->pmax = malloc(sizeof(double) * npar);
  args->idx_B = malloc(sizeof(double) * conf->ninput * 2);
  args->pbest = malloc(sizeof(double) * npar);
  if (!args->pmin || !args->pmax || !args->idx_B || !args->pbest) {
    P_ERR("failed to allocate memory for the fitting parameters\n");
    return BAOFLIT_ERR_MEMORY;
  }

  /* Set prior limits. */
  args->pmin[0] = conf->pmin_a;
  args->pmax[0] = conf->pmax_a;
  for (int i = 0; i < args->num_B; i++) {
    args->pmin[i + 1] = conf->pmin_B[i];
    args->pmax[i + 1] = conf->pmax_B[i];
  }
  if (args->fit_Snl) {
    for (int i = 0; i < conf->ninput; i++) {
      args->pmin[i + conf->num_B + 1] = conf->pmin_Snl[i];
      args->pmax[i + conf->num_B + 1] = conf->pmax_Snl[i];
    }
  }

  /* Set prior indices. */
  for (int i = 0; i < conf->ninput; i++) {
    for (int k = 0; k < 2; k++) {
      int idx = -1;
      for (int j = 0; j < conf->num_B; j++) {
        if (conf->tracer[i][k] == conf->Bfit[j]) {
          idx = j;
          break;
        }
      }
      args->idx_B[i * 2 + k] = idx;
    }
  }

  if (save_param(conf)) return BAOFLIT_ERR_INIT;

  if (conf->verbose) {
    printf("  Number of fitting parameters: %d, names and priors saved\n",
        args->npar);
    printf("  Number of nuisance parameters: %d\n", args->npoly * args->nxi);
  }
  return 0;
}

/******************************************************************************
Function `read_data`:
  Read the data vector from files.
Arguments:
  * `conf`:     structure for storing configurations;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int read_data(const CONF *conf, ARGS *args) {
  double *s, *data;
  size_t n;
  if (read_table(conf->fdata[0], conf->comment, conf->dscol[0], conf->dxicol[0],
      conf->fitmin[0], conf->fitmax[0], &s, &data, &n)) return BAOFLIT_ERR_FILE;
  size_t max, size;
  args->edata[0] = max = size = n;

  for (int i = 1; i < conf->ninput; i++) {
    args->idata[i] = size;
    double *x, *y;
    if (read_table(conf->fdata[i], conf->comment, conf->dscol[i],
        conf->dxicol[i], conf->fitmin[i], conf->fitmax[i], &x, &y, &n))
      return BAOFLIT_ERR_FILE;
    if (size > SIZE_MAX - n) {
      P_ERR("too many samples in the data files\n");
      free(x); free(y); free(s); free(data);
      return BAOFLIT_ERR_FILE;
    }

    /* Enlarge memory if necessary. */
    if (size + n > max) {
      if (SIZE_MAX / 2 < max) {
        P_ERR("too many samples in the data files\n");
        free(x); free(y); free(s); free(data);
        return BAOFLIT_ERR_FILE;
      }
      max <<= 1;
      if (max < size + n) max = size + n;
      double *tmp = realloc(s, sizeof(double) * max);
      if (!tmp) {
        P_ERR("failed to allocate memory for storing separations\n");
        free(x); free(y); free(s); free(data);
        return BAOFLIT_ERR_MEMORY;
      }
      s = tmp;
      tmp = realloc(data, sizeof(double) * max);
      if (!tmp) {
        P_ERR("failed to allocate memory for the data vector\n");
        free(x); free(y); free(s); free(data);
        return BAOFLIT_ERR_MEMORY;
      }
      data = tmp;
    }

    /* Append records to the separation list and data vector. */
    memcpy(s + size, x, sizeof(double) * n);
    memcpy(data + size, y, sizeof(double) * n);
    free(x); free(y);
    size += n;
    args->edata[i] = args->idata[i] + n;
  }

  /* Reduce the memory cost if applicable. */
  double *tmp = realloc(s, sizeof(double) * size);
  if (tmp) s = tmp;
  tmp = realloc(data, sizeof(double) * size);
  if (tmp) data = tmp;

  args->nbin = size;
  args->s = s;
  args->data = data;

  if(conf->verbose) {
    printf("  Data vector constructed with %d %s, total length: %zu\n",
        args->nxi, (args->nxi > 1) ? "2PCFs" : "2PCF", args->nbin);
  }
  return 0;
}

/******************************************************************************
Function `get_cov`:
  Read or compute (upper triangular decomposition of) the covariance matrix.
Arguments:
  * `conf`:     structure for storing configurations;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int get_cov(const CONF *conf, ARGS * args) {
  if (conf->comp_cov) { /* compute covariance matrix using mocks */
    double *ximock = NULL;
    args->nmock = 0;
    for (int i = 0; i < conf->ninput; i++) {
      double *xi = NULL;
      size_t nmock = 0;
      size_t nbin = args->edata[i] - args->idata[i];
      if (read_mocks(conf->fmock[i], conf->comment, conf->mscol[i],
          conf->mxicol[i], conf->fitmin[i], conf->fitmax[i],
          args->s + args->idata[i], nbin, &nmock, &xi))
        return BAOFLIT_ERR_INIT;

      if (!args->nmock) {
        args->nmock = nmock;
        /* Allocate memory for the concatenated mock 2PCFs. */
        ximock = malloc(sizeof(double) * nmock * args->nbin);
        if (!ximock) {
          P_ERR("failed to allocate memory for storing all mock 2PCFs\n");
          free(xi);
          return BAOFLIT_ERR_MEMORY;
        }
      }
      else if (args->nmock != nmock) {
        P_ERR("inconsistent number of mocks in the lists\n");
        free(xi); free(ximock);
        return BAOFLIT_ERR_INIT;
      }

      /* Concatenate mock 2PCFs. */
      for (size_t j = 0; j < nmock; j++) {
        for (size_t k = 0; k < nbin; k++)
          ximock[j + (args->idata[i] + k) * nmock] = xi[j * nbin + k];
      }

      /* Concatenate mock 2PCFs and construct covariance generation matrix. */
      for (size_t j = 0; j < nbin; j++) {
        double sum = 0;
        size_t ioff = (args->idata[i] + j) * nmock;
        for (size_t k = 0; k < nmock; k++) {
          sum += xi[k * nbin + j];
          ximock[ioff + k] = xi[k * nbin + j];
        }
        sum /= nmock;
        for (size_t k = 0; k < nmock; k++) ximock[ioff + k] -= sum;
      }

      free(xi);
    }

    /* Check the dimension. */
    if (args->nmock - args->nbin < 3) {
      P_ERR("%zu mocks are not enough for %zu bins\n", args->nmock, args->nbin);
      return BAOFLIT_ERR_INIT;
    }

    /* QR decomposition: M = Q * R. */
    args->Rcov = qrdcmp(ximock, args->nmock, args->nbin);
    free(ximock);
    if (!args->Rcov) return BAOFLIT_ERR_INIT;

    /* Validate R. */
    for (size_t i = 0; i < args->nbin; i++) {
      double diag = args->Rcov[(((args->nbin << 1) + 1 - i) * i) >> 1];
      if (diag == 0 || !isfinite(diag)) {
        P_ERR("the covariance matrix is singular\n");
        return BAOFLIT_ERR_INIT;
      }
    }

    /* Renormalise R with the Hartlap (2007) correction. */
    double norm = 1 / sqrt(args->nmock - args->nbin - 2);
    size_t len = (args->nbin * (args->nbin + 1)) >> 1;
    for (size_t i = 0; i < len; i++) args->Rcov[i] *= norm;

    /* Save the upper triangular matrix if applicable. */
    if (conf->fcov && *conf->fcov) {
      if (save_Rcov(conf->fcov, args->Rcov, args->nbin, args->nmock))
        return BAOFLIT_ERR_SAVE;
    }
  }
  else {                /* read covariance matrix from file */
    if (read_Rcov(conf->fcov, &args->Rcov, args->nbin, &args->nmock))
      return BAOFLIT_ERR_FILE;
  }

  if (conf->verbose) {
    if (conf->comp_cov) {
      printf("  Covariance matrix constructed with %zu mocks", args->nmock);
      if (conf->fcov && *conf->fcov)
        printf(" and saved to file: `%s'", conf->fcov);
      printf("\n");
    }
    else
      printf("  Covariance matrix read from file successfully\n");
  }
  return 0;
}

/******************************************************************************
Function `init_pk`:
  Initialise the template power spectra.
Arguments:
  * `conf`:     structure for storing configurations;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int init_pk(const CONF *conf, ARGS *args) {
  /* Sample k values for the integration. */
  double *lnk, *fac;
  lnk = fac = NULL;
  double lnk_min = log(conf->kmin);
  double lnk_max = log(conf->kmax);
  switch (args->pkint) {
    case PK_INT_TRAPZ:
      /* Determine the number of k values. */
      args->nk = (size_t) conf->nlogk;
      if (!(lnk = malloc(sizeof(double) * args->nk))) {
        P_ERR("failed to allocate memory for k values\n");
        return BAOFLIT_ERR_MEMORY;
      }
      if (!(fac = malloc(sizeof(double) * args->nk * args->ns))) {
        P_ERR("failed to allocate memory for the integration\n");
        free(lnk);
        return BAOFLIT_ERR_MEMORY;
      }
      /* Evaluate ln(k) bins. */
      for (size_t i = 0; i < args->nk; i++)
        lnk[i] = lnk_min + i * (lnk_max - lnk_min) / (conf->nlogk - 1);
      /* The d ln(k) factor to be applied to the integrand */
      for (size_t i = 0; i < args->nk; i++) fac[i] = lnk[1] - lnk[0];
      /* Deal with the starting and ending points for the trapezoidal rule. */
      fac[0] *= 0.5;
      fac[args->nk - 1] *= 0.5;
      break;
    case PK_INT_LEGAUSS:
      /* Determine the number of k values. */
      args->nk = (size_t) conf->lgorder;
      if (!(lnk = malloc(sizeof(double) * args->nk))) {
        P_ERR("failed to allocate memory for k values\n");
        return BAOFLIT_ERR_MEMORY;
      }
      if (!(fac = malloc(sizeof(double) * args->nk * args->ns))) {
        P_ERR("failed to allocate memory for the integration\n");
        free(lnk);
        return BAOFLIT_ERR_MEMORY;
      }
      /* factors due to the change of inteval for Legendre-Gauss quadrature. */
      double fac1 = (lnk_max - lnk_min) * 0.5;
      double fac2 = (lnk_max + lnk_min) * 0.5;
      /* Evaluate ln(k) bins. */
      for (int i = 0; i < LEGAUSS_LEN_NONZERO(conf->lgorder); i++) {
        /* Look up the abscissas and weights. */
        int idx = i + LEGAUSS_IDX(conf->lgorder);
        double x = legauss_x[idx];
        double w = legauss_w[idx];

        /* Record both positive and negative abscissas. */
        lnk[i] = -x * fac1 + fac2;
        lnk[args->nk - i - 1] = x * fac1 + fac2;
        fac[i] = fac[args->nk - i - 1] = w * fac1;
      }
      /* For odd orders, there is also the abscissas x = 0. */
      if (conf->lgorder & 1) {
        int idx = LEGAUSS_IDX(conf->lgorder)
            + LEGAUSS_LEN_NONZERO(conf->lgorder);
        lnk[LEGAUSS_LEN_NONZERO(conf->lgorder)] = fac2;
        fac[LEGAUSS_LEN_NONZERO(conf->lgorder)] = legauss_w[idx] * fac1;
      }
      break;
  }

  /* Allocate memory for arrays related to pk. */
  args->halfk2 = malloc(sizeof(double) * args->nk);
  args->PBAO = malloc(sizeof(double) * args->nk);
  args->Pnw = malloc(sizeof(double) * args->nk);
  args->Pm = malloc(sizeof(double) * args->nk);
  if (!args->halfk2 || !args->PBAO || !args->Pnw || !args->Pm) {
    P_ERR("failed to allocate memory for pre-computed power spectra terms\n");
    free(lnk); free(fac);
    return BAOFLIT_ERR_MEMORY;
  }
  if (!(args->Pnwt = malloc(sizeof(double *) * args->nxi))) {
    P_ERR("failed to allocate memory for tracer non-wiggle power spectra\n");
    free(lnk); free(fac);
    return BAOFLIT_ERR_MEMORY;
  }
  for (int i = 0; i < args->nxi; i++) args->Pnwt[i] = NULL;
  if (conf->num_nwt) {
    if (!(args->Pnwt[0] = malloc(sizeof(double) * conf->num_nwt * args->nk))) {
      P_ERR("failed to allocate memory for tracer non-wiggle power spectra\n");
      free(lnk); free(fac);
      return BAOFLIT_ERR_MEMORY;
    }
    int j = 0;
    for (int i = 0; i < args->nxi; i++) {
      if (conf->fpnwt[i] && *(conf->fpnwt[i])) {
        args->Pnwt[i] = args->Pnwt[i] + j * args->nk;
        j += 1;
      }
    }
    if (j != conf->num_nwt) {
      P_ERR("unexpected number of non-wiggle tracer power spectrum: %d "
          "rather than %d\n", j, conf->num_nwt);
      return BAOFLIT_ERR_UNKNOWN;
    }
  }

  /* Read the linear matter power spectrum. */
  double *klin, *Plin;
  size_t nlin;
  if (read_table(conf->fplin, conf->comment, 1, 2, 0, DBL_MAX, &klin, &Plin,
      &nlin)) {
    free(lnk); free(fac);
    return BAOFLIT_ERR_INIT;
  }

  if (conf->fpnw) {     /* read the non-wiggle matter power spectrum */
    double *knw, *Pnw;
    size_t nnw;
    if (read_table(conf->fpnw, conf->comment, 1, 2, 0, DBL_MAX, &knw, &Pnw,
        &nnw)) {
      free(lnk); free(fac); free(klin); free(Plin);
      return BAOFLIT_ERR_INIT;
    }

    /* Normalise the non-wiggle matter power spectrum. */
    if (pk_norm(klin, Plin, nlin, knw, Pnw, nnw, conf->knorm)) {
      free(lnk); free(fac); free(klin); free(Plin); free(knw); free(Pnw);
      return BAOFLIT_ERR_INIT;
    }

    /* Interpolate the non-wiggle matter power spectrum. */
    if (pk_interp(knw, Pnw, nnw, lnk, args->Pnw, args->nk)) {
      free(lnk); free(fac); free(klin); free(Plin); free(knw); free(Pnw);
      return BAOFLIT_ERR_INIT;
    }
    free(knw);
    free(Pnw);
  }
  else {                /* compute the non-wiggle matter power spectrum */
    double *knw, *Pnw;
    knw = Pnw = NULL;
    knw = malloc(sizeof(double) * nlin);
    Pnw = malloc(sizeof(double) * nlin);
    if (!knw || !Pnw) {
      P_ERR("failed to allocate memory for the linear non-wiggle "
          "matter power spectrum\n");
      free(lnk); free(fac); free(klin); free(Plin);
      if (knw) free(knw);
      if (Pnw) free(Pnw);
      return BAOFLIT_ERR_MEMORY;
    }
    pk_nw_EH(klin, nlin, conf->hubble, conf->omega_m, conf->omega_b,
        conf->Tcmb, conf->pkns, Pnw);

    /* Normalise the non-wiggle matter power spectrum. */
    if (pk_norm(klin, Plin, nlin, klin, Pnw, nlin, conf->knorm)) {
      free(lnk); free(fac); free(klin); free(Plin); free(knw); free(Pnw);
      return BAOFLIT_ERR_INIT;
    }

    /* Interpolate the non-wiggle matter power spectrum. */
    memcpy(knw, klin, sizeof(double) * nlin);
    if (pk_interp(knw, Pnw, nlin, lnk, args->Pnw, args->nk)) {
      free(lnk); free(fac); free(klin); free(Plin); free(knw); free(Pnw);
      return BAOFLIT_ERR_INIT;
    }
    free(knw);
    free(Pnw);
  }

  /* Interpolate the linear matter power spectrum with wiggles. */
  if (pk_interp(klin, Plin, nlin, lnk, args->PBAO, args->nk)) {
    free(lnk); free(fac); free(klin); free(Plin);
    return BAOFLIT_ERR_INIT;
  }
  free(klin);
  free(Plin);
  for (size_t i = 0; i < args->nk; i++) args->PBAO[i] -= args->Pnw[i];

  /* Read linear non-wiggle tracer power spectrum. */
  if (conf->num_nwt) {
    for (int i = 0; i < args->nxi; i++) {
      if (!(conf->fpnwt[i]) || !(*conf->fpnwt[i])) continue;
      double *kt, *Pt;
      size_t nt;
      if (read_table(conf->fpnwt[i], conf->comment, 1, 2, 0, DBL_MAX,
          &kt, &Pt, &nt)) {
        free(lnk); free(fac);
        return BAOFLIT_ERR_INIT;
      }

      /* Interpolate the linear non-wiggle tracer power spectrum. */
      if (pk_interp(kt, Pt, nt, lnk, args->Pnwt[i], args->nk)) {
        free(lnk); free(fac); free(kt); free(Pt);
        return BAOFLIT_ERR_INIT;
      }
      free(kt);
      free(Pt);

      /* Compute the ratio between the tracer and matter power spectra. */
      for (size_t j = 0; j < args->nk; j++) args->Pnwt[i][j] /= args->Pnw[j];
    }
  }

  /* Pre-compute terms related to k. */
  args->k = lnk;
  args->fac = fac;
  for (size_t i = 0; i < args->nk; i++) {
    args->k[i] = exp(args->k[i]);       /* convert ln(k) to k */
    args->halfk2[i] = args->k[i] * args->k[i];
    args->fac[i] *= exp(-args->halfk2[i] * conf->damp_a * conf->damp_a)
        * args->halfk2[i] * args->k[i] / (2 * M_PI * M_PI);
    args->halfk2[i] *= 0.5;
  }
  for(size_t i = 1; i < args->ns; i++)
    memcpy(fac + i * args->nk, fac, sizeof(double) * args->nk);

  if (conf->verbose)
    printf("  Template power spectra initialised with %zu bins\n", args->nk);
  return 0;
}

/******************************************************************************
Function `sphbessel_j0`:
  Compute the 0-order spherical Bessel function of the first kind.
Arguments:
  * `x`:        the function argument.
Return:
  sinc(x).
******************************************************************************/
static inline double sphbessel_j0(double x) {
  if (x < 0.01) {
    double xx = x * x;
    return 1 - xx * 0x1.5555555555555p-3 + xx * xx * 0x1.1111111111111p-7;
  }
  else return sin(x) / x;
}

/******************************************************************************
Function `init_xi`:
  Initialise the template correlation function.
Arguments:
  * `conf`:     structure for storing configurations;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int init_xi(const CONF *conf, ARGS *args) {
  /* Allocate memory for the model 2PCFs. */
  args->st = malloc(sizeof(double) * args->ns);
  args->sm = malloc(sizeof(double) * args->nbin);
  args->xim = malloc(sizeof(double) * args->nbin);
  if (!args->st || !args->sm || !args->xim) {
    P_ERR("failed to allocate memory for the model 2PCFs\n");
    return BAOFLIT_ERR_MEMORY;
  }
  if (!(args->xit = malloc(sizeof(double *) * args->nxi)) ||
      !(args->xit[0] = malloc(sizeof (double) * args->ns * args->nxi))) {
    P_ERR("failed to alocate memory for the template 2PCFs\n");
    return BAOFLIT_ERR_MEMORY;
  }
  if (!(args->xipp = malloc(sizeof(double *) * args->nxi)) ||
      !(args->xipp[0] = malloc(sizeof(double) * args->ns * 2 * args->nxi))) {
    P_ERR("failed to allocate memory for template 2PCF interpolation\n");
    return BAOFLIT_ERR_MEMORY;
  }
  for (int i = 1; i < args->nxi; i++) {
    args->xit[i] = args->xit[0] + i * args->ns;
    args->xipp[i] = args->xipp[0] + i * args->ns * 2;
  }

  /* Pre-compute values if applicable. */
  for (size_t i = 0; i < args->ns; i++) {
    args->st[i] = conf->smin + conf->ds * i;
    for (size_t j = 0; j < args->nk; j++) {
      args->fac[i * args->nk + j] *= sphbessel_j0(args->st[i] * args->k[j]);
    }
  }

  if (conf->verbose)
    printf("  Model 2PCFs initialised with %zu bins\n", args->ns);
  return 0;
}

/******************************************************************************
Function `init_least_squared`:
  Initialise the matrices and vectors for least squared fitting.
Arguments:
  * `conf`:     structure for storing configurations;
  * `args`:     structure for storing fitting arguments.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int init_least_squared(const CONF *conf, ARGS *args) {
  if (!args->npoly) return 0;           /* no nuisance parameter */
  size_t ntot = (size_t) args->npoly * (size_t) args->nxi;
  /* Allocate memory. */
  args->apoly = malloc(sizeof(double) * ntot);
  args->basis = calloc(args->nbin * ntot, sizeof(double));
  args->LS_Z = malloc(sizeof(double) * args->nbin * ntot);
  if (!args->apoly || !args->basis || !args->LS_Z) {
    P_ERR("failed to allocate memory for least-squared fitting\n");
    return BAOFLIT_ERR_MEMORY;
  }

  /* Construct the basis matrix X for least-squared fitting. */
  for (int i = 0; i < args->nxi; i++) {
    for (int j = 0; j < args->npoly; j++) {
      size_t ioff = (i * args->npoly + j) * args->nbin;
      for (size_t k = args->idata[i]; k < args->edata[i]; k++)
        args->basis[ioff + k] = pow(args->s[k], j - 2);
    }
  }

  /* Construct Y = R^{-T} * X using forward substitution. */
  for (size_t i = 0; i < ntot; i++) {
    fwd_subst(args->Rcov, args->basis + i * args->nbin, args->nbin,
        args->LS_Z + i * args->nbin);
  }

  /* Construct U using QR decomposition: Y = Q1 * U. */
  if (!(args->LS_U = qrdcmp(args->LS_Z, args->nbin, ntot)))
    return BAOFLIT_ERR_INIT;

  /* Construct Y again since it was altered by the in-place QR decomposition. */
  for (size_t i = 0; i < ntot; i++) {
    fwd_subst(args->Rcov, args->basis + i * args->nbin, args->nbin,
        args->LS_Z + i * args->nbin);
  }

  /* Compute U^{-T} * Y using forward substitution. */
  double *col = malloc(sizeof(double) * ntot);
  if (!col) {
    P_ERR("failed to allocate memory for least-squared fit intialisation\n");
    return BAOFLIT_ERR_MEMORY;
  }
  for (size_t i = 0; i < args->nbin; i++) {
    for (size_t j = 0; j < ntot; j++) col[j] = args->LS_Z[j * args->nbin + i];
    fwd_subst(args->LS_U, col, ntot, col);
    for (size_t j = 0; j < ntot; j++) args->LS_Z[j * args->nbin + i] = col[j];
  }
  free(col);

  /* Compute Z = U^{-T} * Y * R^{-T} using backward substitution. */
  for (size_t i = 0; i < ntot; i++) {
    bwd_subst(args->Rcov, args->LS_Z + i * args->nbin, args->nbin,
        args->LS_Z + i * args->nbin);
  }

  /* Reduce the size of basis by verticle aligning non-zero blocks. */
  double *tmp = realloc(args->basis, sizeof(double) * args->nbin * args->npoly);
  if (!tmp) P_WRN("failed to reduce memory for the basis function\n");
  else args->basis = tmp;
  /* Save in the row-first order to reduce cache miss for the fit. */
  for (size_t j = 0; j < args->nbin; j++) {
    size_t ioff = j * args->npoly;
    for (int i = 0; i < args->npoly; i++)
      args->basis[ioff + i] = pow(args->s[j], i - 2);
  }

  if (conf->verbose)
    printf("  Matrices initialised for least-squared fitting\n");
  return 0;
}


/*============================================================================*\
             Interfaces for fitting initialisation and termination
\*============================================================================*/

/******************************************************************************
Function `init_fit`:
  Initialise the fitting process.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for storing fitting arguments.
******************************************************************************/
ARGS *init_fit(const CONF *conf) {
  printf("Initialising the fitting process ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  ARGS *args = init_args(conf);
  if (!args) return NULL;

  if (init_param(conf, args) || read_data(conf, args) || get_cov(conf, args) ||
      init_pk(conf, args) || init_xi(conf, args) ||
      init_least_squared(conf, args)) {
    args_destroy(args);
    return NULL;
  }

  printf(FMT_DONE);
  return args;
}

/******************************************************************************
Function `args_destroy`:
  Release memory allocated for the arguments of the fit.
Arguments:
  * `cf`:       structure for arguments of the fit.
******************************************************************************/
void args_destroy(ARGS *args) {
  if (!args) return;

  if (args->pmin) free(args->pmin);
  if (args->pmax) free(args->pmax);
  if (args->idx_B) free(args->idx_B);
  if (args->pbest) free(args->pbest);

  if (args->data) free(args->data);
  if (args->s) free(args->s);
  if (args->idata) free(args->idata);
  if (args->edata) free(args->edata);
  if (args->Rcov) free(args->Rcov);

  if (args->k) free(args->k);
  if (args->halfk2) free(args->halfk2);
  if (args->fac) free(args->fac);
  if (args->PBAO) free(args->PBAO);
  if (args->Pnw) free(args->Pnw);
  if (args->Pnwt) {
    if (args->Pnwt[0]) free(args->Pnwt[0]);
    free(args->Pnwt);
  }
  if (args->Pm) free(args->Pm);

  if (args->st) free(args->st);
  if (args->xit) {
    if (args->xit[0]) free(args->xit[0]);
    free(args->xit);
  }
  if (args->xipp) {
    if (args->xipp[0]) free (args->xipp[0]);
    free(args->xipp);
  }
  if (args->sm) free(args->sm);
  if (args->xim) free(args->xim);

  if (args->apoly) free(args->apoly);
  if (args->basis) free(args->basis);
  if (args->LS_U) free(args->LS_U);
  if (args->LS_Z) free(args->LS_Z);
  free(args);
}
