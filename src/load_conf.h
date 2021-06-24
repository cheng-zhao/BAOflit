/*******************************************************************************
* load_conf.h: this file is part of the BAOflit program.

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

#ifndef __LOAD_CONF_H__
#define __LOAD_CONF_H__

#include <stdbool.h>

/*============================================================================*\
                   Data structure for storing configurations
\*============================================================================*/

typedef struct {
  char *fconf;          /* Name of the configuration file. */
  int ninput;           /* Number of input 2PCFs. */
  char **fdata;         /* DATA_FILE            */
  char **tracer;        /* TRACER               */
  int *dscol;           /* DATA_SEP_COL         */
  int *dxicol;          /* DATA_XI_COL          */
  double *fitmin;       /* FIT_SEP_MIN          */
  double *fitmax;       /* FIT_SEP_MAX          */
  char *fcov;           /* COV_FILE             */
  bool comp_cov;        /* Indicate whether to compute the covariance matrix. */
  char **fmock;         /* MOCK_LIST            */
  double cov_fac;       /* COV_RESCALE          */
  int *mscol;           /* MOCK_SEP_COL         */
  int *mxicol;          /* MOCK_XI_COL          */
  char comment;         /* FILE_COMMENT         */

  double pmin_a;        /* ALPHA_PRIOR_MIN      */
  double pmax_a;        /* ALPHA_PRIOR_MAX      */
  int num_B;            /* Number of free bias parameters. */
  char *Bfit;           /* TRACER_BIAS_FIT      */
  int Btype;            /* BIAS_PRIOR_TYPE      */
  double *pmin_B;       /* BIAS_PRIOR_MIN       */
  double *pmax_B;       /* BIAS_PRIOR_MAX       */
  double *pcen_B;       /* BIAS_PRIOR_CEN       */
  double *psig_B;       /* BIAS_PRIOR_SIG       */
  int Snltype;          /* SIGMA_TYPE           */
  double *val_Snl;      /* SIGMA_VALUE          */
  double *pmin_Snl;     /* SIGMA_PRIOR_MIN      */
  double *pmax_Snl;     /* SIGMA_PRIOR_MAX      */
  double *pcen_Snl;     /* SIGMA_PRIOR_CEN      */
  double *psig_Snl;     /* SIGMA_PRIOR_SIG      */
#ifdef PARA_MODEL
  int ctype;            /* C_TYPE               */
  double *val_c;        /* C_VALUE              */
  double *pmin_c;       /* C_PRIOR_MIN          */
  double *pmax_c;       /* C_PRIOR_MAX          */
  double *pcen_c;       /* C_PRIOR_CEN          */
  double *psig_c;       /* C_PRIOR_SIG          */
#endif
  int npoly;            /* NUM_NUISANCE         */

  char *fplin;          /* PK_LINEAR            */
  char *fpnw;           /* PK_NOBAO_MATTER      */
#ifndef PARA_MODEL
  char **fpnwt;         /* PK_NOBAO_TRACER      */
  int num_nwt;          /* Number of non-wiggle tracer power spectra. */
  bool *has_nwt;        /* Indicate if non-wiggle tracer power spectra exist. */
#endif
  double knorm;         /* K_NORM               */
  double kmin;          /* K_MIN                */
  double kmax;          /* K_MAX                */
  int pkint;            /* PK_INT_METHOD        */
  int nlogk;            /* NUM_LOG_K            */
  int lgorder;          /* LEGAUSS_ORDER        */
  double damp_a;        /* PK_INT_DAMP          */
  double smin;          /* S_MIN                */
  double smax;          /* S_MAX                */
  double ds;            /* S_BIN_SIZE           */
  int ns;               /* Number of separation bins. */

  double hubble;        /* HUBBLE               */
  double omega_m;       /* OMEGA_M              */
  double omega_b;       /* OMEGA_B              */
  double Tcmb;          /* CMB_TEMP             */
  double pkns;          /* PK_NS                */

  int nlive;            /* NUM_LIVE             */
  double tol;           /* TOLERANCE            */
  bool resume;          /* RESUME               */

  char *oroot;          /* OUTPUT_ROOT          */
  bool verbose;         /* VERBOSE              */
} CONF;


/*============================================================================*\
                      Interface for loading configurations
\*============================================================================*/

/******************************************************************************
Function `load_conf`:
  Read, check, and print configurations.
Arguments:
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  The structure for storing configurations.
******************************************************************************/
CONF *load_conf(const int argc, char *const *argv);

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf);

#endif
