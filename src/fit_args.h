/*******************************************************************************
* fit_args.h: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __FIT_ARGS_H__
#define __FIT_ARGS_H__

#include "load_conf.h"
#include <stdio.h>
#include <stdbool.h>

/*============================================================================*\
               Structures for storing information during the fit
\*============================================================================*/

typedef enum {
  BAOFLIT_PARAM_FIX = -1,
  BAOFLIT_PRIOR_FLAT = 0,
  BAOFLIT_PRIOR_GAUSS = 1
} BAOFLIT_param_t;

typedef enum {
  PK_INT_TRAPZ = 0,
  PK_INT_LEGAUSS = 1
} pk_int_t;

typedef struct {
  int npar;             /* number of fitting parameters                      */
  double *pmin;         /* lower prior limit of the fitting parameters       */
  double *pmax;         /* upper prior limit of the fitting parameters       */
  int num_B;            /* number of free bias parameters                    */
  BAOFLIT_param_t Btype;        /* bias prior type                           */
  const double *Bcen;   /* center of Gaussian priors for the bias parameters */
  const double *Bsig;   /* sigma of Gaussian priors for the bias parameters  */
  int *idx_B;           /* indices of bias parameters for the fitted 2PCFs   */
  BAOFLIT_param_t Snltype;      /* Sigma_nl prior type                       */
  const double *Snlcen; /* center of Gaussian priors for Sigma_nl            */
  const double *Snlsig; /* sigma of Gaussian priors for Sigma_nl             */
#ifdef PARA_MODEL
  const double *Snlval; /* fixed Sigma_nl values                             */
  BAOFLIT_param_t ctype;        /* c prior type                              */
  const double *ccen;   /* center of Gaussian priors for c                   */
  const double *csig;   /* sigma of Gaussian priors for c                    */
  const double *cval;   /* fixed c values                                    */
  int cidx;             /* starting index of c parameters                    */
#endif
  double *pmodel;       /* mean/best-fit/MAP fitting parameters              */
  double *amodel;       /* mean/best-fit/MAP nuisance parameters             */
  double maxlnlike;     /* maximum log-likelihood                            */

  double *data;         /* the data vector                                   */
  double *s;            /* separations for the 2PCFs in the data vector      */
  size_t nbin;          /* total length of the data vector                   */
  int nxi;              /* number of partitions (2PCFs) of the data vector   */
  size_t *idata;        /* starting indices for segments of the data vector  */
  size_t *edata;        /* ending indices for segments of the data vector    */
  size_t nmock;         /* number of mocks for the covariance estimation     */
  double *Rcov;         /* triangular decomposition of the covariance matrix */

  pk_int_t pkint;       /* power spectra integration method                  */
  double *k;            /* abscissas of the power spectra to be integrated   */
  size_t nk;            /* number of sampling points for the power spectra   */
  double *fac;          /* factor multiplied to the power spectra integrand  */
  double *halfk2;       /* pre-computed k^2 / 2                              */
  double *PBAO;         /* BAO wiggles from the linear matter power spectra  */
  double *Pnw;          /* wiggle-free linear matter power spectrum          */
#ifndef PARA_MODEL
  double **Pnwt;        /* wiggle-free linear tracer power spectra           */
  const bool *has_nwt;  /* whether wiggle-free tracer power spectra exist    */
#endif
  double *Pm;           /* model power spectrum                              */

  size_t ns;            /* number of separation bins for the template 2PCF   */
  double *st;           /* separations for the template 2PCFs                */
  double *st2;          /* squared separations for the template 2PCFs        */
  double **xit;         /* template 2PCFs                                    */
  double **xipp;        /* second derivative of the 2PCFs for interpolation  */
  double *sm;           /* model separation bins shifted by alpha            */
  double *xim;          /* model 2PCF to be compared with the data vector    */

  int npoly;            /* number of nuisance parameter                      */
  double *apoly;        /* the nuisance parameters                           */
  double *basis;        /* the basis function for least-squared fitting      */
  double *LS_U;         /* triangular decomposition of the design matrix     */
  double *LS_Z;         /* matrix for least-squared fitting                  */
} ARGS;


/*============================================================================*\
                      Interfaces for the fitting arguments
\*============================================================================*/

/******************************************************************************
Function `init_fit`:
  Initialise the fitting process.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for storing fitting arguments.
******************************************************************************/
ARGS *init_fit(const CONF *conf);

/******************************************************************************
Function `args_destroy`:
  Release memory allocated for the arguments of the fit.
Arguments:
  * `cf`:       structure for arguments of the fit.
******************************************************************************/
void args_destroy(ARGS *args);

#endif
