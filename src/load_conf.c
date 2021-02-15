/*******************************************************************************
* load_conf.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "fit_args.h"
#include "legauss.h"
#include "libcfg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
/* Check existence of configuration parameters. */
#define CHECK_EXIST_PARAM(name, cfg, var)                       \
  if (!cfg_is_set((cfg), (var))) {                              \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return BAOFLIT_ERR_CFG;                                     \
  }
#define CHECK_EXIST_ARRAY(name, cfg, var, num)                  \
  if (!(num = cfg_get_size((cfg), (var)))) {                    \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return BAOFLIT_ERR_CFG;                                     \
  }

/* Check length of array. */
#define CHECK_ARRAY_LENGTH(name, cfg, var, fmt, num, nexp)      \
  if ((num) < (nexp)) {                                         \
    P_ERR("too few elements of " FMT_KEY(name) "\n");           \
    return BAOFLIT_ERR_CFG;                                     \
  }                                                             \
  if ((num) > (nexp)) {                                         \
    P_WRN("omitting the following " FMT_KEY(name) ":");         \
    for (int i = (nexp); i < (num); i++)                        \
      fprintf(stderr, " " fmt, (var)[i]);                       \
    fprintf(stderr, "\n");                                      \
  }

#define CHECK_STR_ARRAY_LENGTH(name, cfg, var, num, nexp)       \
  if ((num) < (nexp)) {                                         \
    P_ERR("too few elements of " FMT_KEY(name) "\n");           \
    return BAOFLIT_ERR_CFG;                                     \
  }                                                             \
  if ((num) > (nexp)) {                                         \
    P_WRN("omitting the following " FMT_KEY(name) ":\n");       \
    for (int i = (nexp); i < (num); i++)                        \
      fprintf(stderr, "  %s\n", (var)[i]);                      \
  }

/* Release memory for configuration parameters. */
#define FREE_ARRAY(x)           {if(x) free(x);}
#define FREE_STR_ARRAY(x)       {if(x) {if (*(x)) free(*(x)); free(x);}}

/* Print the warning and error messages. */
#define P_CFG_WRN(cfg)  cfg_pwarn(cfg, stderr, FMT_WARN);
#define P_CFG_ERR(cfg)  {                                       \
  cfg_perror(cfg, stderr, FMT_ERR);                             \
  cfg_destroy(cfg);                                             \
  return NULL;                                                  \
}


/*============================================================================*\
                    Functions called via command line flags
\*============================================================================*/

/******************************************************************************
Function `usage`:
  Print the usage of command line options.
******************************************************************************/
static void usage(void *args) {
  (void) args;
  printf("Usage: " BAOFLIT_CODE_NAME " [OPTION]\n\
Perform multi-tracer BAO fit with the MultiNest sampler.\n\
  -h, --help\n\
        Display this message and exit\n\
  -t, --template\n\
        Print a template configuration file to the standard output and exit\n\
  -c, --conf            " FMT_KEY(CONFIG_FILE) "     String\n\
        Specify the configuration file (default: `%s')\n\
  -d, --data            " FMT_KEY(DATA_FILE) "       String array\n\
        Set the input 2PCFs to be fitted\n\
  -T, --tracer          " FMT_KEY(TRACER) "          String array\n\
        Specify tracers of the 2PCFs, indicated by pairs of letters\n\
  -s, --data-s          " FMT_KEY(DATA_SEP_COL) "    Integer array\n\
        Column numbers (starting from 1) for separations in " FMT_KEY(DATA_FILE) "\n\
  -y, --data-xi         " FMT_KEY(DATA_XI_COL) "     Integer array\n\
        Column numbers (starting from 1) for 2PCFs in " FMT_KEY(DATA_FILE) "\n\
      --fit-min         " FMT_KEY(FIT_SEP_MIN) "     Double array\n\
        Minimum separations of 2PCFs for the fit\n\
      --fit-max         " FMT_KEY(FIT_SEP_MAX) "     Double array\n\
        Maximum separations of 2PCFs for the fit\n\
  -C, --cov             " FMT_KEY(COV_FILE) "        String\n\
        Set the covariance matrix (generated using this program) for the fit\n\
  -m, --mock            " FMT_KEY(MOCK_LIST) "       String array\n\
        Set the lists of mock 2PCFs used for covariance matrix estimations\n\
  -S, --mock-s          " FMT_KEY(MOCK_SEP_COL) "    Integer array\n\
        Column numbers (starting from 1) for separations in " FMT_KEY(MOCK_LIST) "\n\
  -Y, --mock-xi         " FMT_KEY(MOCK_XI_COL) "     Integer array\n\
        Column numbers (starting from 1) for 2PCFs in " FMT_KEY(MOCK_LIST) "\n\
      --comment         " FMT_KEY(FILE_COMMENT) "    Character\n\
        Specify comment symbols for input files\n\
      --alpha-min       " FMT_KEY(ALPHA_PRIOR_MIN) " Double\n\
        Set the lower prior limit of alpha\n\
      --alpha-max       " FMT_KEY(ALPHA_PRIOR_MAX) " Double\n\
        Set the upper prior limit of alpha\n\
  -B, --B-fit           " FMT_KEY(TRACER_BIAS_FIT) " Character array\n\
        Specify tracers with free bias parameters (B)\n\
      --B-prior-type    " FMT_KEY(BIAS_PRIOR_TYPE) " Integer\n\
        Prior type of the bias parameters\n\
      --B-min           " FMT_KEY(BIAS_PRIOR_MIN) "  Double array\n\
        Set the lower prior limits of B\n\
      --B-max           " FMT_KEY(BIAS_PRIOR_MAX) "  Double array\n\
        Set the upper prior limits of B\n\
      --B-center        " FMT_KEY(BIAS_PRIOR_CEN) "  Double array\n\
        Set the central values for Gaussian priors of B\n\
      --B-sigma         " FMT_KEY(BIAS_PRIOR_SIG) "  Double array\n\
        Set the standard deviations for Gaussian priors of B\n\
      --Snl-value       " FMT_KEY(SIGMA_VALUE) "     Double array\n\
        Set fixed values of the BAO dampling parameters (Sigma_nl)\n\
      --Snl-min         " FMT_KEY(SIGMA_PRIOR_MIN) " Double array\n\
        Set the lower prior limits of Sigma_nl\n\
      --Snl-max         " FMT_KEY(SIGMA_PRIOR_MAX) " Double array\n\
        Set the upper prior limits of Sigma_nl\n\
      --num-nuisance    " FMT_KEY(NUM_NUISANCE) "    Integer\n\
        Number of nuisance parameters for least-squared fits\n\
  -p, --pk-lin          " FMT_KEY(PK_LINEAR) "       String\n\
        Set the input linear matter power spectrum\n\
  -P, --pk-nw           " FMT_KEY(PK_NOBAO_MATTER) " String\n\
        Set the input linear non-wiggle matter power spectrum\n\
      --pk-tracer       " FMT_KEY(PK_NOBAO_TRACER) " String array\n\
        Set the input linear non-wiggle tracer power spectrum templates\n\
      --k-norm          " FMT_KEY(K_NORM) "          Double\n\
        The upper limit of k used for non-wiggle power spectra normalisations\n\
      --k-min           " FMT_KEY(K_MIN) "           Double\n\
        Minimum k value for sampling template power spectra\n\
      --k-max           " FMT_KEY(K_MAX) "           Double\n\
        Maximum k value for sampling template power spectra\n\
      --pk-int          " FMT_KEY(PK_INT_METHOD) "   Integer\n\
        Method for integrating the template power spectra\n\
      --num-log-k       " FMT_KEY(NUM_LOG_K) "       Integer\n\
        Number of log(k) bins for the trapezoidal integration\n\
      --lg-orderr       " FMT_KEY(LEGAUSS_ORDER) "   Integer\n\
        Order of Legendre-Gauss quadrature for power spectra integration\n\
      --pk-int-damp     " FMT_KEY(PK_INT_DAMP) "     Double\n\
        The damping parameter for power spectra integrations\n\
      --s-min           " FMT_KEY(S_MIN) "           Double\n\
        Minimum separation for sampling the model 2PCFs\n\
      --s-max           " FMT_KEY(S_MAX) "           Double\n\
        Maximum separation for sampling the model 2PCFs\n\
      --s-step          " FMT_KEY(S_BIN_SIZE) "      Double\n\
        Bin size of separations for the model 2PCFs\n\
      --hubble          " FMT_KEY(HUBBLE) "          Double\n\
        The dimensionless Hubble parameter for generating " FMT_KEY(PK_NOBAO_MATTER) "\n\
      --omega-m         " FMT_KEY(OMEGA_M) "         Double\n\
        Density parameter of matter at redshift 0\n\
      --omega-b         " FMT_KEY(OMEGA_B) "         Double\n\
        Density parameter of baryons at redshift 0\n\
      --CMB-temp        " FMT_KEY(CMB_TEMP) "        Double\n\
        Temperature of cosmic microwave background in Kelvin\n\
      --pk-ns           " FMT_KEY(PK_NS) "           Double\n\
        Scalar index of the primordial power spectrum\n\
  -n, --num-live        " FMT_KEY(NUM_LIVE) "        Integer\n\
        Number of live points for the MultiNest sampler\n\
  -e, --tolerance       " FMT_KEY(TOLERANCE) "       Double\n\
        Tolerance (stopping criteria) for the MultiNest sampler\n\
  -r, --resume          " FMT_KEY(RESUME) "          Boolean\n\
        Indicate whether to resume the previous MultiNest run\n\
  -o, --output          " FMT_KEY(OUTPUT_ROOT) "     String\n\
        Set the basename of outputs\n\
  -v, --verbose         " FMT_KEY(VERBOSE) "         Boolean\n\
        Indicate whether to display detailed standard outputs\n\
Consult the -t option for more information on the parameters\n\
Github repository: https://github.com/cheng-zhao/BAOflit\n\
Licence: MIT\n",
    DEFAULT_CONF_FILE);
  exit(0);
}

/******************************************************************************
Function `conf_template`:
  Print a template configuration file.
******************************************************************************/
static void conf_template(void *args) {
  (void) args;
  printf("# Configuration file for " BAOFLIT_CODE_NAME " (default: `"
DEFAULT_CONF_FILE "').\n\
# Format: keyword = value # comment\n\
#     or: keyword = [element1, element2]\n\
#    see: https://github.com/cheng-zhao/libcfg for details.\n\
# Some of the entries allow expressions, see\n\
#         https://github.com/cheng-zhao/libast for details.\n\
# NOTE that command line options have priority over this file.\n\
# Unnecessary entries can be left unset.\n\
\n\
#############################################################\n\
#  Specifications of the data vector and covariance matrix  #\n\
#############################################################\n\
\n\
DATA_FILE       = \n\
    # Filename for the 2-point correlation functions (2PCFs) to be fitted.\n\
    # String or string array. If multiple data files are supplied,\n\
    # the corresponding data vectors are joined for the simultaneous fit.\n\
TRACER          = \n\
    # Tracers of the 2PCFs, indicated by pairs of case-sensitive letters.\n\
    # String or string array, same dimension as `DATA_FILE`.\n\
    # If the two letters are different, the corresponding 2PCF is treated as\n\
    # cross correlations, while identical letters indicate auto 2PCFs.\n\
    # Biases of cross 2PCFs are computed as products of the two tracer biases.\n\
DATA_SEP_COL    = \n\
    # Column numbers (starting from 1) for separations in each `DATA_FILE`.\n\
    # Integer or integer array, same dimension as `DATA_FILE`.\n\
DATA_XI_COL     = \n\
    # Column numbers (starting from 1) for 2PCFs.\n\
    # Integer or integer array, same dimension as `DATA_FILE`.\n\
FIT_SEP_MIN     = \n\
FIT_SEP_MAX     = \n\
    # Minimum and maximum separations of each 2PCF for the fit.\n\
    # Double-precision numbers or arrays, same dimension as `DATA_FILE`.\n\
COV_FILE        = \n\
    # Filename for the covariance matrix.\n\
    # If it is unset, or set but the file does not exist, the covariance matrix\n\
    # will be computed using `MOCK_LIST`, and saved to `COV_FILE`.\n\
MOCK_LIST       = \n\
    # Filename for lists of 2PCFs from mocks.\n\
    # String or string array, same dimension as `DATA_FILE`.\n\
    # The order of mock lists should be consistent with that of `DATA_FILE`.\n\
    # For each mock list, the format of the 2PCF files should be identical.\n\
MOCK_SEP_COL    = \n\
    # Column numbers for separations of mock 2PCFs in each list.\n\
    # Integer or integer array, same dimension as `DATA_FILE`.\n\
MOCK_XI_COL     = \n\
    # Column numbers for 2PCFs of mocks in each list.\n\
    # Integer or integer array, same dimension as `DATA_FILE`.\n\
FILE_COMMENT    = \n\
    # Character indicating lines to be skipped for input files (unset: '%c%s.\n\
\n\
########################################\n\
#  Settings of the fitting parameters  #\n\
########################################\n\
\n\
ALPHA_PRIOR_MIN = \n\
ALPHA_PRIOR_MAX = \n\
    # Flat prior ranges of alpha.\n\
    # Double-precision numbers.\n\
TRACER_BIAS_FIT = \n\
    # Tracers with biases as free parameters, indicated by letters.\n\
    # Character or character array.\n\
    # The letters should be taken from the ones defined in `TRACER`.\n\
    # If the bias of a tracer is not supplied here, its value is set to 1.\n\
BIAS_PRIOR_TYPE = \n\
    # Prior type of the bias parameters (B, unset: %d).\n\
    # Integer, allowed values are:\n\
    # * 0: flat prior;\n\
    # * 1: Gaussian prior.\n\
BIAS_PRIOR_MIN  = \n\
BIAS_PRIOR_MAX  = \n\
    # Prior ranges of B, for both the flat and Gaussian priors.\n\
    # Double-precision number or array, same dimension as `TRACER_BIAS_FIT`.\n\
BIAS_PRIOR_CEN  = \n\
BIAS_PRIOR_SIG  = \n\
    # Mean and standard deviation for the Gaussian prior of B.\n\
    # Double-precision number or array, same dimension as `TRACER_BIAS_FIT`.\n\
SIGMA_VALUE     = \n\
    # The BAO damping parameter (Sigma_nl).\n\
    # Double-precision number or array, same dimension as `DATA_FILE`.\n\
    # If unset, Sigma_nl is fitted as free parameter.\n\
SIGMA_PRIOR_MIN = \n\
SIGMA_PRIOR_MAX = \n\
    # Flat prior ranges of Sigma_nl.\n\
    # Double-precision number or array, same dimension as `DATA_FILE`.\n\
    # They are used only if `SIGMA_VALUE` is unset.\n\
NUM_NUISANCE    = \n\
    # Number of noisance (polynomial) parameters for all 2PCFs (unset: %d).\n\
\n\
#################################\n\
#  Specifications of the model  #\n\
#################################\n\
\n\
PK_LINEAR       = \n\
    # Filename for the linear matter power spectrum.\n\
    # The first two columns must be k and P(k).\n\
PK_NOBAO_MATTER = \n\
    # Filename for the linear non-wiggle (no BAO) matter power spectrum.\n\
    # The first two columns must be k and P(k).\n\
    # If unset, the non-wiggle power spectrum will be computed following\n\
    # Eisenstein & Hu 1998 (arXiv:astro-ph/9709112).\n\
PK_NOBAO_TRACER = \n\
    # Filename for the linear non-wiggle tracer power spectrum.\n\
    # The first two columns must be k and P(k).\n\
    # String or string array, same dimension as `DATA_FILE`.\n\
    # It is used for modelling scale-dependent bias of tracers, see\n\
    # Zhao et al. in preparation.\n\
    # Leave it unset to disable the tracer bias model for all data, or set an\n\
    # empty string (\"\") to disable the model for a given 2PCF in `DATA_FILE`.\n\
    # See also Variu et al. in preparation.\n\
K_NORM          = \n\
    # The non-wiggle power spectra are normalised below this k value.\n\
    # Double-precision number.\n\
K_MIN           = \n\
K_MAX           = \n\
    # Minimum and maximum k values for integrating the power spectra.\n\
    # Double-precision numbers.\n\
PK_INT_METHOD   = \n\
    # Method for integrating the power spectra (unset: %d).\n\
    # Integer, allowed values are:\n\
    # * 0: trapezoidal integration with `NUM_LOG_K` sample points;\n\
    # * 1: Legendre-Gauss quadrature with the order of `LEGAUSS_ORDER`.\n\
NUM_LOG_K       = \n\
    # Integer, number of log(k) bins for integrating the power spectra.\n\
LEGAUSS_ORDER   = \n\
    # Order of Legendre-Gauss quadrature for integrating the power spectra.\n\
    # Integer between 4 and 32.\n\
PK_INT_DAMP     = \n\
    # The damping parameter for integrating power spectra.\n\
    # Double-precision number. See the a factor in Eq. (27) of\n\
    # Xu et al. 2012 (arXiv:1202.0091).\n\
S_MIN           = \n\
S_MAX           = \n\
S_BIN_SIZE      = \n\
    # Lower and upper limits, and width of separation bins for the model 2PCF.\n\
    # Double-precision numbers.\n\
\n\
###############################################################\n\
#  Cosmological parameters for the non-wiggle power spectrum  #\n\
###############################################################\n\
\n\
HUBBLE          = \n\
    # The non-deminsional Hubble parameter (h).\n\
OMEGA_M         = \n\
    # Density parameter of matter at redshift 0.\n\
OMEGA_B         = \n\
    # Density parameter of baryons at redshift 0.\n\
CMB_TEMP        = \n\
    # Temperature of cosmic microwave background in Kelvin.\n\
PK_NS           = \n\
    # Scalar index of the primordial power spectrum.\n\
\n\
############################################\n\
#  Configurations for parameter inference  #\n\
############################################\n\
\n\
NUM_LIVE        = \n\
    # Integer, number of live points for the MultiNest sampler.\n\
TOLERANCE       = \n\
    # Tolerance (stopping criteria) for the MultiNest sampler.\n\
    # Double-precision number.\n\
RESUME          = \n\
    # Boolean option, indicate whether to resume the MultiNest run (unset: %c).\n\
\n\
##############################\n\
#  Settings for the outputs  #\n\
##############################\n\
\n\
OUTPUT_ROOT     = \n\
    # String, basename of the outputs to be written by multinest.\n\
VERBOSE         = \n\
    # Boolean option, indicate whether to show detailed outputs (unset: %c).\n",
      DEFAULT_COMMENT ? DEFAULT_COMMENT : '\'', DEFAULT_COMMENT ? "')" : ")",
      DEFAULT_BIAS_PRIOR, DEFAULT_NUM_NUISANCE, DEFAULT_PK_INT_METHOD,
      DEFAULT_RESUME ? 'T' : 'F', DEFAULT_VERBOSE ? 'T' : 'F');
  exit(0);
}


/*============================================================================*\
                      Function for reading configurations
\*============================================================================*/

/******************************************************************************
Function `conf_init`:
  Initialise the structure for storing configurations.
Return:
  Address of the structure.
******************************************************************************/
static CONF *conf_init(void) {
  CONF *conf = calloc(1, sizeof *conf);
  if (!conf) return NULL;
  conf->fdata = conf->fmock = conf->fpnwt = conf->tracer = NULL;
  conf->dscol = conf->dxicol = conf->mscol = conf->mxicol = NULL;
  conf->fitmin = conf->fitmax = NULL;
  conf->pmin_B = conf->pmax_B = conf->pcen_B = conf->psig_B = NULL;
  conf->val_Snl = conf->pmin_Snl = conf->pmax_Snl = NULL;
  conf->fconf = conf->fcov = conf->fplin = conf->fpnw = conf->Bfit = NULL;
  conf->oroot = NULL;
  return conf;
}

/******************************************************************************
Function `conf_read`:
  Read configurations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  Interface of libcfg.
******************************************************************************/
static cfg_t *conf_read(CONF *conf, const int argc, char *const *argv) {
  if (!conf) {
    P_ERR("the structure for configurations is not initialised\n");
    return NULL;
  }
  cfg_t *cfg = cfg_init();
  if (!cfg) P_CFG_ERR(cfg);

  /* Functions to be called via command line flags. */
  const int nfunc = 2;
  const cfg_func_t funcs[] = {
    {   'h',        "help",             usage,          NULL},
    {   't',    "template",     conf_template,          NULL}
  };

  /* Configuration parameters. */
  const int npar = 47;
  const cfg_param_t params[] = {
    {'c', "conf"        , "CONFIG_FILE"    , CFG_DTYPE_STR , &conf->fconf   },
    {'d', "data"        , "DATA_FILE"      , CFG_ARRAY_STR , &conf->fdata   },
    {'T', "tracer"      , "TRACER"         , CFG_ARRAY_STR , &conf->tracer  },
    {'s', "data-s"      , "DATA_SEP_COL"   , CFG_ARRAY_INT , &conf->dscol   },
    {'y', "data-xi"     , "DATA_XI_COL"    , CFG_ARRAY_INT , &conf->dxicol  },
    { 0 , "fit-min"     , "FIT_SEP_MIN"    , CFG_ARRAY_DBL , &conf->fitmin  },
    { 0 , "fit-max"     , "FIT_SEP_MAX"    , CFG_ARRAY_DBL , &conf->fitmax  },
    {'C', "cov"         , "COV_FILE"       , CFG_DTYPE_STR , &conf->fcov    },
    {'m', "mock"        , "MOCK_LIST"      , CFG_ARRAY_STR , &conf->fmock   },
    {'S', "mock-s"      , "MOCK_SEP_COL"   , CFG_ARRAY_INT , &conf->mscol   },
    {'Y', "mock-xi"     , "MOCK_XI_COL"    , CFG_ARRAY_INT , &conf->mxicol  },
    { 0 , "comment"     , "FILE_COMMENT"   , CFG_DTYPE_CHAR, &conf->comment },
    { 0 , "alpha-min"   , "ALPHA_PRIOR_MIN", CFG_DTYPE_DBL , &conf->pmin_a  },
    { 0 , "alpha-max"   , "ALPHA_PRIOR_MAX", CFG_DTYPE_DBL , &conf->pmax_a  },
    {'B', "B-fit"       , "TRACER_BIAS_FIT", CFG_ARRAY_CHAR, &conf->Bfit    },
    { 0 , "B-prior-type", "BIAS_PRIOR_TYPE", CFG_DTYPE_INT , &conf->Btype   },
    { 0 , "B-min"       , "BIAS_PRIOR_MIN" , CFG_ARRAY_DBL , &conf->pmin_B  },
    { 0 , "B-max"       , "BIAS_PRIOR_MAX" , CFG_ARRAY_DBL , &conf->pmax_B  },
    { 0 , "B-center"    , "BIAS_PRIOR_CEN" , CFG_ARRAY_DBL , &conf->pcen_B  },
    { 0 , "B-sigma"     , "BIAS_PRIOR_SIG" , CFG_ARRAY_DBL , &conf->psig_B  },
    { 0 , "Snl-value"   , "SIGMA_VALUE"    , CFG_ARRAY_DBL , &conf->val_Snl },
    { 0 , "Snl-min"     , "SIGMA_PRIOR_MIN", CFG_ARRAY_DBL , &conf->pmin_Snl},
    { 0 , "Snl-max"     , "SIGMA_PRIOR_MAX", CFG_ARRAY_DBL , &conf->pmax_Snl},
    { 0 , "num-nuisance", "NUM_NUISANCE"   , CFG_DTYPE_INT , &conf->npoly   },
    {'p', "pk-lin"      , "PK_LINEAR"      , CFG_DTYPE_STR , &conf->fplin   },
    {'P', "pk-nw"       , "PK_NOBAO_MATTER", CFG_DTYPE_STR , &conf->fpnw    },
    { 0 , "pk-tracer"   , "PK_NOBAO_TRACER", CFG_ARRAY_STR , &conf->fpnwt   },
    { 0 , "k-norm"      , "K_NORM"         , CFG_DTYPE_DBL , &conf->knorm   },
    { 0 , "k-min"       , "K_MIN"          , CFG_DTYPE_DBL , &conf->kmin    },
    { 0 , "k-max"       , "K_MAX"          , CFG_DTYPE_DBL , &conf->kmax    },
    { 0 , "pk-int"      , "PK_INT_METHOD"  , CFG_DTYPE_INT , &conf->pkint   },
    { 0 , "num-log-k"   , "NUM_LOG_K"      , CFG_DTYPE_INT , &conf->nlogk   },
    { 0 , "lg-order"    , "LEGAUSS_ORDER"  , CFG_DTYPE_INT , &conf->lgorder },
    { 0 , "pk-int-damp" , "PK_INT_DAMP"    , CFG_DTYPE_DBL , &conf->damp_a  },
    { 0 , "s-min"       , "S_MIN"          , CFG_DTYPE_DBL , &conf->smin    },
    { 0 , "s-max"       , "S_MAX"          , CFG_DTYPE_DBL , &conf->smax    },
    { 0 , "s-step"      , "S_BIN_SIZE"     , CFG_DTYPE_DBL , &conf->ds      },
    { 0 , "hubble"      , "HUBBLE"         , CFG_DTYPE_DBL , &conf->hubble  },
    { 0 , "omega-m"     , "OMEGA_M"        , CFG_DTYPE_DBL , &conf->omega_m },
    { 0 , "omega-b"     , "OMEGA_B"        , CFG_DTYPE_DBL , &conf->omega_b },
    { 0 , "CMB-temp"    , "CMB_TEMP"       , CFG_DTYPE_DBL , &conf->Tcmb    },
    { 0 , "pk-ns"       , "PK_NS"          , CFG_DTYPE_DBL , &conf->pkns    },
    {'n', "num-live"    , "NUM_LIVE"       , CFG_DTYPE_INT , &conf->nlive   },
    {'e', "tolerance"   , "TOLERANCE"      , CFG_DTYPE_DBL , &conf->tol     },
    {'r', "resume"      , "RESUME"         , CFG_DTYPE_BOOL, &conf->resume  },
    {'o', "output"      , "OUTPUT_ROOT"    , CFG_DTYPE_STR , &conf->oroot   },
    {'v', "verbose"     , "VERBOSE"        , CFG_DTYPE_BOOL, &conf->verbose }
  };

  /* Register functions and parameters. */
  if (cfg_set_funcs(cfg, funcs, nfunc)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);
  if (cfg_set_params(cfg, params, npar)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read configurations from command line options. */
  int optidx;
  if (cfg_read_opts(cfg, argc, argv, BAOFLIT_PRIOR_CMD, &optidx))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read parameters from configuration file. */
  if (!cfg_is_set(cfg, &conf->fconf)) conf->fconf = DEFAULT_CONF_FILE;
  if (access(conf->fconf, R_OK))
    P_WRN("cannot access the configuration file: `%s'\n", conf->fconf);
  else if (cfg_read_file(cfg, conf->fconf, BAOFLIT_PRIOR_FILE)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  return cfg;
}


/*============================================================================*\
                      Functions for parameter verification
\*============================================================================*/

/******************************************************************************
Function `check_input`:
  Check whether an input file can be read.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int check_input(const char *fname, const char *key) {
  if (!fname || *fname == '\0') {
    P_ERR("the input " FMT_KEY(%s) " is not set\n", key);
    return BAOFLIT_ERR_CFG;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot access " FMT_KEY(%s) ": `%s'\n", key, fname);
    return BAOFLIT_ERR_FILE;
  }
  return 0;
}

/******************************************************************************
Function `check_outdir`:
  Check whether an output directory is accessible.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_outdir(char *fname, const char *key) {
  if (!fname || *fname == '\0') {
    P_ERR("the output " FMT_KEY(%s) " is not set\n", key);
    return BAOFLIT_ERR_CFG;
  }

  char *end;
  if ((end = strrchr(fname, BAOFLIT_PATH_SEP)) != NULL) {
    *end = '\0';
    if (access(fname, X_OK)) {
      P_ERR("cannot access the directory `%s'\n", fname);
      return BAOFLIT_ERR_FILE;
    }
    *end = BAOFLIT_PATH_SEP;
  }
  return 0;
}

/******************************************************************************
Function `conf_verify`:
  Verify configuration parameters.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int conf_verify(const cfg_t *cfg, CONF *conf) {
  int e, num;

  /* DATA_FILE */
  CHECK_EXIST_ARRAY(DATA_FILE, cfg, &conf->fdata, conf->ninput);
  if (conf->ninput > BAOFLIT_MAX_TRACER) {
    P_ERR("number of elements in " FMT_KEY(DATA_FILE) " cannot exceed %d\n",
        BAOFLIT_MAX_TRACER);
    return BAOFLIT_ERR_CFG;
  }
  for (int i = 0; i < conf->ninput; i++) {
    if ((e = check_input(conf->fdata[i], "DATA_FILE"))) return e;
  }

  /* TRACER */
  CHECK_EXIST_ARRAY(TRACER, cfg, &conf->tracer, num);
  CHECK_STR_ARRAY_LENGTH(TRACER, cfg, conf->tracer, num, conf->ninput);
  /* Simple validation. */
  for (int i = 0; i < conf->ninput; i++) {
    char *s = conf->tracer[i];
    if (!isalpha(s[0]) || !isalpha(s[1])) {
      P_ERR("invalid " FMT_KEY(TRACER) ": %s\n", s);
      return BAOFLIT_ERR_CFG;
    }
  }
  /* Check duplicates. */
  for (int i = 0; i < conf->ninput - 1; i++) {
    for (int j = i + 1; j < conf->ninput; j++) {
      if (conf->tracer[i][0] == conf->tracer[j][0] &&
          conf->tracer[i][1] == conf->tracer[j][1]) {
        P_ERR("duplicate " FMT_KEY(TRACER) ": %s\n", conf->tracer[i]);
        return BAOFLIT_ERR_CFG;
      }
    }
  }

  /* DATA_SEP_COL and DATA_XI_COL */
  CHECK_EXIST_ARRAY(DATA_SEP_COL, cfg, &conf->dscol, num);
  CHECK_ARRAY_LENGTH(DATA_SEP_COL, cfg, conf->dscol, "%d", num, conf->ninput);
  for (int i = 0; i < conf->ninput; i++) {
    if (conf->dscol[i] <= 0 || conf->dscol[i] >= BAOFLIT_MAX_FILE_COL) {
      P_ERR(FMT_KEY(DATA_SEP_COL) " must be positive and smaller than %d\n",
          BAOFLIT_MAX_FILE_COL);
      return BAOFLIT_ERR_CFG;
    }
  }
  CHECK_EXIST_ARRAY(DATA_XI_COL, cfg, &conf->dxicol, num);
  CHECK_ARRAY_LENGTH(DATA_XI_COL, cfg, conf->dxicol, "%d", num, conf->ninput);
  for (int i = 0; i < conf->ninput; i++) {
    if (conf->dxicol[i] <= 0 || conf->dxicol[i] >= BAOFLIT_MAX_FILE_COL) {
      P_ERR(FMT_KEY(DATA_XI_COL) " must be positive and smaller than %d\n",
          BAOFLIT_MAX_FILE_COL);
      return BAOFLIT_ERR_CFG;
    }
    if (conf->dxicol[i] == conf->dscol[i]) {
      P_ERR(FMT_KEY(DATA_XI_COL) " must be different from "
          FMT_KEY(DATA_SEP_COL) "\n");
      return BAOFLIT_ERR_CFG;
    }
  }

  /* FIT_SEP_MIN */
  CHECK_EXIST_ARRAY(FIT_SEP_MIN, cfg, &conf->fitmin, num);
  CHECK_ARRAY_LENGTH(FIT_SEP_MIN, cfg, conf->fitmin, OFMT_DBL, num,
      conf->ninput);
  double fitmin = DBL_MAX;
  for (int i = 0; i < conf->ninput; i++) {
    if (conf->fitmin[i] < 0) {
      P_ERR(FMT_KEY(FIT_SEP_MIN) " must be positive\n");
      return BAOFLIT_ERR_CFG;
    }
    if (fitmin > conf->fitmin[i]) fitmin = conf->fitmin[i];
  }

  /* FIT_SEP_MAX */
  CHECK_EXIST_ARRAY(FIT_SEP_MAX, cfg, &conf->fitmax, num);
  CHECK_ARRAY_LENGTH(FIT_SEP_MAX, cfg, conf->fitmax, OFMT_DBL, num,
      conf->ninput);
  double fitmax = -DBL_MAX;
  for (int i = 0; i < conf->ninput; i++) {
    if (conf->fitmin[i] >= conf->fitmax[i]) {
      P_ERR(FMT_KEY(FIT_SEP_MAX) " must be larger than " FMT_KEY(FIT_SEP_MIN)
          "\n");
      return BAOFLIT_ERR_CFG;
    }
    if (fitmax < conf->fitmax[i]) fitmax = conf->fitmax[i];
  }

  /* COV_FILE */
  conf->comp_cov = true;
  if (cfg_is_set(cfg, &conf->fcov)) {
    if (!access(conf->fcov, R_OK)) conf->comp_cov = false;
    else if ((e = check_outdir(conf->fcov, "COV_FILE"))) return e;
  }

  if (conf->comp_cov) {
    /* MOCK_LIST */
    CHECK_EXIST_ARRAY(MOCK_LIST, cfg, &conf->fmock, num);
    CHECK_STR_ARRAY_LENGTH(MOCK_LIST, cfg, conf->fmock, num, conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      if ((e = check_input(conf->fmock[i], "MOCK_LIST"))) return e;
    }

    /* MOCK_SEP_COL and MOCK_XI_COL*/
    CHECK_EXIST_ARRAY(MOCK_SEP_COL, cfg, &conf->mscol, num);
    CHECK_ARRAY_LENGTH(MOCK_SEP_COL, cfg, conf->mscol, "%d", num, conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      if (conf->mscol[i] <= 0 || conf->mscol[i] >= BAOFLIT_MAX_FILE_COL) {
        P_ERR(FMT_KEY(MOCK_SEP_COL) " must be positive and smaller than %d\n",
            BAOFLIT_MAX_FILE_COL);
        return BAOFLIT_ERR_CFG;
      }
    }
    CHECK_EXIST_ARRAY(MOCK_XI_COL, cfg, &conf->mxicol, num);
    CHECK_ARRAY_LENGTH(MOCK_XI_COL, cfg, conf->mxicol, "%d", num, conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      if (conf->mxicol[i] <= 0 || conf->mxicol[i] >= BAOFLIT_MAX_FILE_COL) {
        P_ERR(FMT_KEY(MOCK_XI_COL) " must be positive and smaller than %d\n",
            BAOFLIT_MAX_FILE_COL);
        return BAOFLIT_ERR_CFG;
      }
      if (conf->mxicol[i] == conf->mscol[i]) {
        P_ERR(FMT_KEY(MOCK_XI_COL) " must be different from "
            FMT_KEY(MOCK_SEP_COL) "\n");
        return BAOFLIT_ERR_CFG;
      }
    }
  }

  /* FILE_COMMENT */
  if (!cfg_is_set(cfg, &conf->comment)) conf->comment = DEFAULT_COMMENT;
  if (conf->comment && !isgraph(conf->comment)) {
    P_ERR("invalid " FMT_KEY(COMMENT) ", ASCII code: %d\n", conf->comment);
    return BAOFLIT_ERR_CFG;
  }

  /* ALPHA_PRIOR_MIN and ALPHA_PRIOR_MAX */
  CHECK_EXIST_PARAM(ALPHA_PRIOR_MIN, cfg, &conf->pmin_a);
  if (conf->pmin_a <= 0) {
    P_ERR(FMT_KEY(ALPHA_PRIOR_MIN) " must be positive\n");
    return BAOFLIT_ERR_CFG;
  }
  CHECK_EXIST_PARAM(ALPHA_PRIOR_MAX, cfg, &conf->pmax_a);
  if (conf->pmin_a >= conf->pmax_a) {
    P_ERR(FMT_KEY(ALPHA_PRIOR_MAX) " must be larger than "
        FMT_KEY(ALPHA_PRIOR_MIN) "\n");
    return BAOFLIT_ERR_CFG;
  }

  /* TRACER_BIAS_FIT */
  if ((conf->num_B = cfg_get_size(cfg, &conf->Bfit))) {
    /* Brute-force validation. */
    for (int i = 0; i < conf->num_B; i++) {
      char B = conf->Bfit[i];
      if (!isalpha(B)) {
        P_ERR(FMT_KEY(TRACER_BIAS_FIT) " must be capital letters\n");
        return BAOFLIT_ERR_CFG;
      }
      bool exist = false;
      for (int j = 0; j < conf->ninput; j++) {
        if (conf->tracer[j][0] == B || conf->tracer[j][1] == B) {
          exist = true;
          break;
        }
      }
      if (!exist) {
        P_ERR("tracer in " FMT_KEY(TRACER_BIAS_FIT) " not found in "
            FMT_KEY(TRACER) ": %c\n", B);
        return BAOFLIT_ERR_CFG;
      }
    }
    /* Check duplicates. */
    for (int i = 0; i < conf->num_B - 1; i++) {
      for (int j = i + 1; j < conf->num_B; j++) {
        if (conf->Bfit[i] == conf->Bfit[j]) {
          P_ERR("duplicate " FMT_KEY(TRACER_BIAS_FIT) ": %c\n", conf->Bfit[i]);
          return BAOFLIT_ERR_CFG;
        }
      }
    }

    /* BIAS_PRIOR_TYPE */
    if (!cfg_is_set(cfg, &conf->Btype)) conf->Btype = DEFAULT_BIAS_PRIOR;
    switch (conf->Btype) {
      case BIAS_PRIOR_FLAT:
        /* BIAS_PRIOR_MIN and BIAS_PRIOR_MAX */
        CHECK_EXIST_ARRAY(BIAS_PRIOR_MIN, cfg, &conf->pmin_B, num);
        CHECK_ARRAY_LENGTH(BIAS_PRIOR_MIN, cfg, conf->pmin_B, OFMT_DBL,
            num, conf->num_B);
        CHECK_EXIST_ARRAY(BIAS_PRIOR_MAX, cfg, &conf->pmax_B, num);
        CHECK_ARRAY_LENGTH(BIAS_PRIOR_MAX, cfg, conf->pmax_B, OFMT_DBL,
            num, conf->num_B);
        for (int i = 0; i < conf->num_B; i++) {
          if (conf->pmin_B[i] >= conf->pmax_B[i]) {
            P_ERR(FMT_KEY(BIAS_PRIOR_MAX) " must be larger than "
                FMT_KEY(BIAS_PRIOR_MIN) "\n");
            return BAOFLIT_ERR_CFG;
          }
        }
        break;
      case BIAS_PRIOR_GAUSS:
        /* BIAS_PRIOR_MIN and BIAS_PRIOR_MAX */
        CHECK_EXIST_ARRAY(BIAS_PRIOR_MIN, cfg, &conf->pmin_B, num);
        CHECK_ARRAY_LENGTH(BIAS_PRIOR_MIN, cfg, conf->pmin_B, OFMT_DBL,
            num, conf->num_B);
        CHECK_EXIST_ARRAY(BIAS_PRIOR_MAX, cfg, &conf->pmax_B, num);
        CHECK_ARRAY_LENGTH(BIAS_PRIOR_MAX, cfg, conf->pmax_B, OFMT_DBL,
            num, conf->num_B);
        for (int i = 0; i < conf->num_B; i++) {
          if (conf->pmin_B[i] >= conf->pmax_B[i]) {
            P_ERR(FMT_KEY(BIAS_PRIOR_MAX) " must be larger than "
                FMT_KEY(BIAS_PRIOR_MIN) "\n");
            return BAOFLIT_ERR_CFG;
          }
        }
        /* BIAS_PRIOR_CEN and BIAS_PRIOR_SIG */
        CHECK_EXIST_ARRAY(BIAS_PRIOR_CEN, cfg, &conf->pcen_B, num);
        CHECK_ARRAY_LENGTH(BIAS_PRIOR_CEN, cfg, conf->pcen_B, OFMT_DBL,
            num, conf->num_B);
        CHECK_EXIST_ARRAY(BIAS_PRIOR_SIG, cfg, &conf->psig_B, num);
        CHECK_ARRAY_LENGTH(BIAS_PRIOR_SIG, cfg, conf->psig_B, OFMT_DBL,
            num, conf->num_B);
        for (int i = 0; i < conf->num_B; i++) {
          if (conf->psig_B[i] <= 0) {
            P_ERR(FMT_KEY(BIAS_PRIOR_SIG) " must be positive\n");
            return BAOFLIT_ERR_CFG;
          }
          double low = conf->pcen_B[i] - conf->psig_B[i] * BAOFLIT_WARN_SIGMA;
          double high = conf->pcen_B[i] + conf->psig_B[i] * BAOFLIT_WARN_SIGMA;
          if (conf->pmin_B[i] > low || conf->pmax_B[i] < high) {
            P_WRN(FMT_KEY(BIAS_PRIOR_MIN) " (" OFMT_DBL ") or "
                FMT_KEY(BIAS_PRIOR_MAX) " (" OFMT_DBL ") is inside " OFMT_DBL
                " sigma range of the Gaussian prior: [" OFMT_DBL "," OFMT_DBL
                "]\n", conf->pmin_B[i], conf->pmax_B[i],
                (double) BAOFLIT_WARN_SIGMA, low, high);
          }
        }
        break;
      default:
        P_ERR("invalid " FMT_KEY(BIAS_PRIOR_TYPE) ": %d\n", conf->Btype);
        return BAOFLIT_ERR_CFG;
    }
  }

  /* SIGMA_VALUE */
  if ((num = cfg_get_size(cfg, &conf->val_Snl))) {
    CHECK_ARRAY_LENGTH(SIGMA_VALUE, cfg, conf->val_Snl, OFMT_DBL, num,
        conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      if (conf->val_Snl[i] < 0) {
        P_ERR(FMT_KEY(SIGMA_VALUE) " must be non-negative\n");
        return BAOFLIT_ERR_CFG;
      }
    }
  }
  else {
    /* SIGMA_PRIOR_MIN */
    CHECK_EXIST_ARRAY(SIGMA_PRIOR_MIN, cfg, &conf->pmin_Snl, num);
    CHECK_ARRAY_LENGTH(SIGMA_PRIOR_MIN, cfg, conf->pmin_Snl, OFMT_DBL,
        num, conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      if (conf->pmin_Snl[i] < 0) {
        P_ERR(FMT_KEY(SIGMA_PRIOR_MIN) " must be non-negative\n");
        return BAOFLIT_ERR_CFG;
      }
    }

    /* SIGMA_PRIOR_MAX */
    CHECK_EXIST_ARRAY(SIGMA_PRIOR_MAX, cfg, &conf->pmax_Snl, num);
    CHECK_ARRAY_LENGTH(SIGMA_PRIOR_MAX, cfg, conf->pmax_Snl, OFMT_DBL,
        num, conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      if (conf->pmin_Snl[i] >= conf->pmax_Snl[i]) {
        P_ERR(FMT_KEY(SIGMA_PRIOR_MAX) " must be larger than "
            FMT_KEY(SIGMA_PRIOR_MIN) "\n");
        return BAOFLIT_ERR_CFG;
      }
    }
  }

  /* NUM_NUISANCE */
  if (!cfg_is_set(cfg, &conf->npoly)) conf->npoly = DEFAULT_NUM_NUISANCE;
  if (conf->npoly < 0) {
    P_ERR(FMT_KEY(NUM_NUISANCE) " must be non-negative\n");
    return BAOFLIT_ERR_CFG;
  }

  /* PK_LINEAR */
  CHECK_EXIST_PARAM(PK_LINEAR, cfg, &conf->fplin);
  if ((e = check_input(conf->fplin, "PK_LINEAR"))) return e;

  /* PK_NOBAO_MATTER */
  if (cfg_is_set(cfg, &conf->fpnw)) {
    if ((e = check_input(conf->fpnw, "PK_NOBAO_MATTER"))) return e;
  }
  else {
    /* HUBBLE */
    CHECK_EXIST_PARAM(HUBBLE, cfg, &conf->hubble);
    if (conf->hubble <= 0) {
      P_ERR(FMT_KEY(HUBBLE) " must be positive\n");
      return BAOFLIT_ERR_CFG;
    }

    /* OMEGA_M */
    CHECK_EXIST_PARAM(OMEGA_M, cfg, &conf->omega_m);
    if (conf->omega_m <= 0 || conf->omega_m > 1) {
      P_ERR(FMT_KEY(OMEGA_M) " must be positive and not larger than 1\n");
      return BAOFLIT_ERR_CFG;
    }

    /* OMEGA_B */
    CHECK_EXIST_PARAM(OMEGA_B, cfg, &conf->omega_b);
    if (conf->omega_b <= 0 || conf->omega_b > conf->omega_m) {
      P_ERR(FMT_KEY(OMEGA_B) " must be positive and not larger than "
          FMT_KEY(OMEGA_M) "\n");
      return BAOFLIT_ERR_CFG;
    }

    /* CMB_TEMP */
    CHECK_EXIST_PARAM(CMB_TEMP, cfg, &conf->Tcmb);
    if (conf->Tcmb <= 0) {
      P_ERR(FMT_KEY(CMB_TEMP) " must be positive\n");
      return BAOFLIT_ERR_CFG;
    }

    /* PK_NS */
    CHECK_EXIST_PARAM(PK_NS, cfg, &conf->pkns);
    if (conf->pkns <= 0) {
      P_ERR(FMT_KEY(PK_NS) " must be positive\n");
      return BAOFLIT_ERR_CFG;
    }
  }

  /* PK_NOBAO_TRACER */
  conf->num_nwt = 0;
  if ((num = cfg_get_size(cfg, &conf->fpnwt))) {
    CHECK_STR_ARRAY_LENGTH(PK_NOBAO_TRACER, cfg, conf->fpnwt, num,
        conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      char *fname = conf->fpnwt[i];
      /* Omit empty strings. */
      if (((fname[0] == '\'' && fname[1] == '\'') ||
          (fname[0] == '"' && fname[1] == '"')) && fname[2] == '\0')
        *fname = '\0';
      else if (*fname) {
        if ((e = check_input(fname, "PK_NOBAO_TRACER"))) return e;
        conf->num_nwt += 1;
      }
    }
  }

  /* K_NORM */
  CHECK_EXIST_PARAM(K_NORM, cfg, &conf->knorm);
  if (conf->knorm <= 0) {
    P_ERR(FMT_KEY(K_NORM) " must be positive\n");
    return BAOFLIT_ERR_CFG;
  }

  /* K_MIN and K_MAX */
  CHECK_EXIST_PARAM(K_MIN, cfg, &conf->kmin);
  if (conf->kmin < 0) {
    P_ERR(FMT_KEY(K_MIN) " must be non-negative\n");
    return BAOFLIT_ERR_CFG;
  }
  CHECK_EXIST_PARAM(K_MAX, cfg, &conf->kmax);
  if (conf->kmin >= conf->kmax) {
    P_ERR(FMT_KEY(K_MAX) " must be larger than " FMT_KEY(K_MIN) "\n");
    return BAOFLIT_ERR_CFG;
  }

  /* PK_INT_METHOD */
  if (!cfg_is_set(cfg, &conf->pkint)) conf->pkint = DEFAULT_PK_INT_METHOD;
  switch (conf->pkint) {
    case PK_INT_TRAPZ:
      /* NUM_LOG_K */
      CHECK_EXIST_PARAM(NUM_LOG_K, cfg, &conf->nlogk);
      if (conf->nlogk <= 1) {
        P_ERR(FMT_KEY(NUM_LOG_K) " must be larger than 1\n");
        return BAOFLIT_ERR_CFG;
      }
      break;
    case PK_INT_LEGAUSS:
      P_ERR("Gaussian quadrature is currently not available\n"
          "Please set " FMT_KEY(PK_INT_METHOD) " to %d\n", PK_INT_TRAPZ);
      return BAOFLIT_ERR_CFG;
      /* LEGAUSS_ORDER */
      /*
      CHECK_EXIST_PARAM(LEGAUSS_ORDER, cfg, &conf->lgorder);
      if (conf->lgorder < LEGAUSS_MIN_ORDER ||
          conf->lgorder > LEGAUSS_MAX_ORDER) {
        P_ERR(FMT_KEY(LEGAUSS_ORDER) " must be between %d and %d\n",
            LEGAUSS_MIN_ORDER, LEGAUSS_MAX_ORDER);
        return BAOFLIT_ERR_CFG;
      }
      break;
      */
    default:
      P_ERR("invalid" FMT_KEY(PK_INT_METHOD) ": %d\n", conf->pkint);
      return BAOFLIT_ERR_CFG;
  }

  /* PK_INT_DAMP */
  CHECK_EXIST_PARAM(PK_INT_DAMP, cfg, &conf->damp_a);
  if (conf->damp_a < 0) {
    P_ERR(FMT_KEY(PK_INT_DAMP) " must be non-negative\n");
    return BAOFLIT_ERR_CFG;
  }

  /* S_MIN, S_MAX, and S_BIN_SIZE */
  CHECK_EXIST_PARAM(S_MIN, cfg, &conf->smin);
  if (conf->smin <= 0) {
    P_ERR(FMT_KEY(S_MIN) " must be positive\n");
    return BAOFLIT_ERR_CFG;
  }
  CHECK_EXIST_PARAM(S_MAX, cfg, &conf->smax);
  CHECK_EXIST_PARAM(S_BIN_SIZE, cfg, &conf->ds);
  if (conf->ds <= 0) {
    P_ERR(FMT_KEY(S_BIN_SIZE) " must be positive\n");
    return BAOFLIT_ERR_CFG;
  }
  if (conf->smin + conf->ds > conf->smax + DOUBLE_TOL) {
    P_ERR(FMT_KEY(S_MIN) " + " FMT_KEY(S_BIN_SIZE) " must not be larger than "
        FMT_KEY(S_MAX) "\n");
    return BAOFLIT_ERR_CFG;
  }
  double smax = conf->smin;
  conf->ns = 1;
  while (smax < conf->smax - DOUBLE_TOL) {
    smax += conf->ds;
    if (++conf->ns > BAOFLIT_MAX_SEP_BIN) {
      P_ERR("too many separations bins given " FMT_KEY(S_MIN) ", "
          FMT_KEY(S_MAX) ", and " FMT_KEY(S_BIN_SIZE) "\n");
      return BAOFLIT_ERR_CFG;
    }
  }
  if (smax > conf->smax + DOUBLE_TOL) {
    P_WRN("reduce " FMT_KEY(S_MAX) " to " OFMT_DBL " given "
        FMT_KEY(S_MIN) " and " FMT_KEY(S_BIN_SIZE) "\n", smax);
  }
  conf->smax = smax;

  /* Check if the separation range of the model covers the fitting range. */
  if (conf->smin * conf->pmax_a > fitmin ||
      conf->smax * conf->pmin_a < fitmax) {
    P_ERR("The separation range of model cannot cover the fitting range "
        "after shifting by alpha\n");
    return BAOFLIT_ERR_CFG;
  }

  /* NUM_LIVE */
  CHECK_EXIST_PARAM(NUM_LIVE, cfg, &conf->nlive);
  if (conf->nlive <= 0) {
    P_ERR(FMT_KEY(NUM_LIVE) " must be positive\n");
    return BAOFLIT_ERR_CFG;
  }

  /* TOLERANCE */
  CHECK_EXIST_PARAM(TOLERANCE, cfg, &conf->tol);
  if (conf->tol <= 0) {
    P_ERR(FMT_KEY(TOLERANCE) " must be positive\n");
    return BAOFLIT_ERR_CFG;
  }

  /* RESUME */
  if (!cfg_is_set(cfg, &conf->resume)) conf->resume = DEFAULT_RESUME;

  /* OUTPUT_ROOT */
  CHECK_EXIST_PARAM(OUTPUT_ROOT, cfg, &conf->oroot);
  if ((e = check_outdir(conf->oroot, "OUTPUT_ROOT"))) return e;
  /* Extend the length of OUTPUT_ROOT if necessary. */
  size_t len = strlen(conf->oroot);
  if (len + BAOFLIT_MN_SUFF_MAX >= BAOFLIT_MN_FNAME_LEN) {
    P_ERR(FMT_KEY(OUTPUT_ROOT) " is too long.\n");
  }
  else {
    char *tmp = realloc(conf->oroot, sizeof(char) * BAOFLIT_MN_FNAME_LEN);
    if (!tmp) {
      P_ERR("cannot extend " FMT_KEY(OUTPUT_ROOT) ".\n");
      return BAOFLIT_ERR_MEMORY;
    }
    conf->oroot = tmp;
  }

  /* VERBOSE */
  if (!cfg_is_set(cfg, &conf->verbose)) conf->verbose = DEFAULT_VERBOSE;

  return 0;
}


/*============================================================================*\
                      Function for printing configurations
\*============================================================================*/

/******************************************************************************
Function `conf_print`:
  Print configuration parameters.
Arguments:
  * `conf`:     structure for storing configurations.
******************************************************************************/
static void conf_print(const CONF *conf) {
  /* Configuration file */
  printf("\n  CONFIG_FILE     = %s", conf->fconf);

  /* Input data files. */
  printf("\n  DATA_FILE       = %s", conf->fdata[0]);
  for (int i = 1; i < conf->ninput; i++)
    printf("\n                    %s", conf->fdata[i]);
  printf("\n  TRACER          = %s", conf->tracer[0]);
  for (int i = 1; i < conf->ninput; i++) printf(" , %s", conf->tracer[i]);
  printf("\n  DATA_SEP_COL    = %d", conf->dscol[0]);
  for (int i = 1; i < conf->ninput; i++) printf(" , %d", conf->dscol[i]);
  printf("\n  DATA_XI_COL     = %d", conf->dxicol[0]);
  for (int i = 1; i < conf->ninput; i++) printf(" , %d", conf->dxicol[i]);

  printf("\n  FIT_SEP_MIN     = " OFMT_DBL, conf->fitmin[0]);
  for (int i = 1; i < conf->ninput; i++)
    printf(" , " OFMT_DBL, conf->fitmin[i]);
  printf("\n  FIT_SEP_MAX     = " OFMT_DBL, conf->fitmax[0]);
  for (int i = 1; i < conf->ninput; i++)
    printf(" , " OFMT_DBL, conf->fitmax[i]);

  if (conf->comp_cov) {
    if (conf->fcov) printf("\n  COV_FILE        = <W> %s", conf->fcov);
    printf("\n  MOCK_LIST       = %s", conf->fmock[0]);
    for (int i = 1; i < conf->ninput; i++)
      printf("\n                    %s", conf->fmock[i]);
    printf("\n  MOCK_SEP_COL    = %d", conf->mscol[0]);
    for (int i = 1; i < conf->ninput; i++) printf(" , %d", conf->mscol[i]);
    printf("\n  MOCK_XI_COL     = %d", conf->mxicol[0]);
    for (int i = 1; i < conf->ninput; i++) printf(" , %d", conf->mxicol[i]);
  }
  else printf("\n  COV_FILE        = <R> %s", conf->fcov);
  if (conf->comment) printf("\n  FILE_COMMENT    = '%c'", conf->comment);

  /* Fitting parameters. */
  printf("\n  ALPHA_PRIOR_MIN = " OFMT_DBL, conf->pmin_a);
  printf("\n  ALPHA_PRIOR_MAX = " OFMT_DBL, conf->pmax_a);

  if (conf->num_B) {
    printf("\n  TRACER_BIAS_FIT = '%c'", conf->Bfit[0]);
    for (int i = 1; i < conf->num_B; i++) printf(" , '%c'", conf->Bfit[i]);
    printf("\n  BIAS_PRIOR_TYPE = %d", conf->Btype);
    switch (conf->Btype) {
      case BIAS_PRIOR_FLAT:
        printf(" (flat)");
        printf("\n  BIAS_PRIOR_MIN  = " OFMT_DBL, conf->pmin_B[0]);
        for (int i = 1; i < conf->num_B; i++)
          printf(" , " OFMT_DBL, conf->pmin_B[i]);
        printf("\n  BIAS_PRIOR_MAX  = " OFMT_DBL, conf->pmax_B[0]);
        for (int i = 1; i < conf->num_B; i++)
          printf(" , " OFMT_DBL, conf->pmax_B[i]);
        break;
      case BIAS_PRIOR_GAUSS:
        printf(" (Gaussian)");
        printf("\n  BIAS_PRIOR_MIN  = " OFMT_DBL, conf->pmin_B[0]);
        for (int i = 1; i < conf->num_B; i++)
          printf(" , " OFMT_DBL, conf->pmin_B[i]);
        printf("\n  BIAS_PRIOR_MAX  = " OFMT_DBL, conf->pmax_B[0]);
        for (int i = 1; i < conf->num_B; i++)
          printf(" , " OFMT_DBL, conf->pmax_B[i]);
        printf("\n  BIAS_PRIOR_CEN  = " OFMT_DBL, conf->pcen_B[0]);
        for (int i = 1; i < conf->num_B; i++)
          printf(" , " OFMT_DBL, conf->pcen_B[i]);
        printf("\n  BIAS_PRIOR_SIG  = " OFMT_DBL, conf->psig_B[0]);
        for (int i = 1; i < conf->num_B; i++)
          printf(" , " OFMT_DBL, conf->psig_B[i]);
        break;
      default:
        P_ERR("unexpected " FMT_KEY(BIAS_PRIOR_TYPE) ": %d\n", conf->Btype);
        return;
    }
  }

  if (conf->val_Snl) {
    printf("\n  SIGMA_VALUE     = " OFMT_DBL, conf->val_Snl[0]);
    for (int i = 1; i < conf->ninput; i++)
      printf(" , " OFMT_DBL, conf->val_Snl[i]);
  }
  else {
    printf("\n  SIGMA_PRIOR_MIN = " OFMT_DBL, conf->pmin_Snl[0]);
    for (int i = 1; i < conf->ninput; i++)
      printf(" , " OFMT_DBL, conf->pmin_Snl[i]);
    printf("\n  SIGMA_PRIOR_MAX = " OFMT_DBL, conf->pmax_Snl[0]);
    for (int i = 1; i < conf->ninput; i++)
      printf(" , " OFMT_DBL, conf->pmax_Snl[i]);
  }
  printf("\n  NUM_NUISANCE    = %d", conf->npoly);

  /* Model evaluation. */
  printf("\n  PK_LINEAR       = %s", conf->fplin);
  if (conf->fpnw) printf("\n  PK_NOBAO_MATTER = %s", conf->fpnw);
  if (conf->fpnwt) {
    printf("\n  PK_NOBAO_TRACER = %s", conf->fpnwt[0]);
    for (int i = 1; i < conf->ninput; i++) {
      if (*conf->fpnwt[i]) printf("\n                    %s", conf->fpnwt[i]);
      else printf("\n                    \"\"");
    }
  }

  printf("\n  K_NORM          = " OFMT_DBL, conf->knorm);
  printf("\n  K_MIN           = " OFMT_DBL, conf->kmin);
  printf("\n  K_MAX           = " OFMT_DBL, conf->kmax);

  printf("\n  PK_INT_METHOD   = %d", conf->pkint);
  switch (conf->pkint) {
    case PK_INT_TRAPZ:
      printf(" (trapezoidal integration)");
      printf("\n  NUM_LOG_K       = %d", conf->nlogk);
      break;
    case PK_INT_LEGAUSS:
      printf(" (Legendre-Gauss quadrature)");
      printf("\n  LEGAUSS_ORDER   = %d", conf->lgorder);
      break;
    default:
      P_ERR("unexpected " FMT_KEY(BIAS_PRIOR_TYPE) ": %d\n", conf->Btype);
      return;
  }

  printf("\n  PK_INT_DAMP     = " OFMT_DBL, conf->damp_a);
  printf("\n  S_MIN           = " OFMT_DBL, conf->smin);
  printf("\n  S_MAX           = " OFMT_DBL, conf->smax);
  printf("\n  S_BIN_SIZE      = " OFMT_DBL, conf->ds);

  /* Cosmological parameters. */
  printf("\n  HUBBLE          = " OFMT_DBL, conf->hubble);
  printf("\n  OMEGA_M         = " OFMT_DBL, conf->omega_m);
  printf("\n  OMEGA_B         = " OFMT_DBL, conf->omega_b);
  printf("\n  CMB_TEMP        = " OFMT_DBL, conf->Tcmb);
  printf("\n  PK_NS           = " OFMT_DBL, conf->pkns);

  /* Parameter inference. */
  printf("\n  NUM_LIVE        = %d", conf->nlive);
  printf("\n  TOLERANCE       = " OFMT_DBL, conf->tol);
  printf("\n  RESUME          = %c", conf->resume ? 'T' : 'F');

  /* Outputs. */
  printf("\n  OUTPUT_ROOT     = %s", conf->oroot);
}


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
CONF *load_conf(const int argc, char *const *argv) {
  CONF *conf = conf_init();
  if (!conf) return NULL;

  cfg_t *cfg = conf_read(conf, argc, argv);
  if (!cfg) {
    conf_destroy(conf);
    return NULL;
  }

  printf("Loading configurations ...");
  fflush(stdout);

  if (conf_verify(cfg, conf)) {
    if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
    conf_destroy(conf);
    cfg_destroy(cfg);
    return NULL;
  }

  if (conf->verbose) conf_print(conf);

  if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
  cfg_destroy(cfg);

  printf(FMT_DONE);
  return conf;
}

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf) {
  if (!conf) return;
  FREE_STR_ARRAY(conf->fdata);
  FREE_STR_ARRAY(conf->tracer);
  FREE_ARRAY(conf->dscol);
  FREE_ARRAY(conf->dxicol);
  FREE_ARRAY(conf->fitmin);
  FREE_ARRAY(conf->fitmax);
  FREE_ARRAY(conf->fcov);
  FREE_STR_ARRAY(conf->fmock);
  FREE_ARRAY(conf->mscol);
  FREE_ARRAY(conf->mxicol);
  FREE_ARRAY(conf->Bfit);
  FREE_ARRAY(conf->pmin_B);
  FREE_ARRAY(conf->pmax_B);
  FREE_ARRAY(conf->pcen_B);
  FREE_ARRAY(conf->psig_B);
  FREE_ARRAY(conf->val_Snl);
  FREE_ARRAY(conf->pmin_Snl);
  FREE_ARRAY(conf->pmax_Snl);
  FREE_ARRAY(conf->fplin);
  FREE_ARRAY(conf->fpnw);
  FREE_STR_ARRAY(conf->fpnwt);
  FREE_ARRAY(conf->oroot);
  free(conf);
}
