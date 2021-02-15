/*******************************************************************************
* define.h: this file is part of the BAOflit program.

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

#ifndef __DEFINE_H__
#define __DEFINE_H__

/*============================================================================*\
                 Definitions of mathematical/physical constants
\*============================================================================*/
#ifndef M_PI
#define M_PI            0x1.921fb54442d18p+1    /* PI */
#endif
#ifndef M_E
#define M_E             0x1.5bf0a8b145769p+1    /* e */
#endif
#define DOUBLE_EPSILON  1e-16   /* ~ machine epsilon for double numbers */
#define DOUBLE_TOL      1e-8    /* tolerance for double number comparison */

/*============================================================================*\
                     Definitions for the MultiNest fitting
\*============================================================================*/
#define BAOFLIT_MN_IS           1       /* do Nested Importance Sampling?    */
#define BAOFLIT_MN_MMODAL       0       /* do mode separation?               */
#define BAOFLIT_MN_CEFF         0       /* run in constant efficiency mode?  */
#define BAOFLIT_MN_EFR          1       /* required efficiency               */
#define BAOFLIT_MN_UPD          2000    /* no. of iterations between updates */
#define BAOFLIT_MN_ZTOL         -1e90   /* threshold for ignoring logZ       */
#define BAOFLIT_MN_MAXMODE      100     /* maximum number of modes           */
#define BAOFLIT_MN_SEED         100     /* random seed for MultiNest         */
#define BAOFLIT_MN_STDOUT       0       /* feedback on standard output?      */
#define BAOFLIT_MN_INITMPI      0       /* initialise MPI routines?          */
#define BAOFLIT_MN_LOGZERO     -DBL_MAX /* threshold for ignoring loglike    */
#define BAOFLIT_MN_MAXIT        0       /* maximum number of iterations      */
#define BAOFLIT_MN_PWRAP        0       /* periodic boundary conditions?     */
#define BAOFLIT_MN_FNAME_LEN    1000    /* no. of characters of filename     */
#define BAOFLIT_MN_SUFF_MAX     32      /* maximum length of output suffices */

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
/* Default value for unset parameters. */
#define DEFAULT_CONF_FILE               "baoflit.conf"
#define DEFAULT_COMMENT                 '#'
#define DEFAULT_BIAS_PRIOR              0
#define DEFAULT_NUM_NUISANCE            3
#define DEFAULT_PK_INT_METHOD           0
#define DEFAULT_RESUME                  true
#define DEFAULT_VERBOSE                 true

/* Priority of parameters from different sources. */
#define BAOFLIT_PRIOR_CMD               5
#define BAOFLIT_PRIOR_FILE              1

/* Limits of configuration parameters. */
#define BAOFLIT_MAX_SEP_BIN     8192    /* maximum number of separation bins */
#define BAOFLIT_MAX_TRACER      100     /* maximum allowed number of tracers */
#define BAOFLIT_MAX_FILE_COL    1024    /* maximum allowed number of columns */

/* Warn if the Gaussian prior range is inside the following times the sigma. */
#define BAOFLIT_WARN_SIGMA      2.0

/*============================================================================*\
                            Definitions for file IO
\*============================================================================*/
#define BAOFLIT_PATH_SEP        '/'     /* separator for file paths     */
#define BAOFLIT_FILE_CHUNK      1048576 /* chunk size for ASCII file IO */
#define BAOFLIT_MAX_CHUNK       INT_MAX /* maximum allowed chunk size   */
/* Initial number of objects allocated for the catalogs.        */
#define BAOFLIT_DATA_INIT_NUM   128
/* Comment symbol for the output files.                         */
#define BAOFLIT_SAVE_COMMENT    '#'

/*============================================================================*\
                            Other runtime constants
\*============================================================================*/
#define BAOFLIT_CODE_NAME       "BAOflit"       /* name of the program */
#define BAOFLIT_DEFAULT_BIAS    1.0     /* value for fixed bias parameters */

/*============================================================================*\
                     Definitions for the format of outputs
\*============================================================================*/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m"          /* Magenta "Warning" */
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m"            /* Red "Error"       */
#define FMT_EXIT "\x1B[31;1mExit:\x1B[0m"               /* Red "Exit"        */
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"  /* Green "DONE"      */
#define FMT_FAIL "\r\x1B[70C[\x1B[31;1mFAIL\x1B[0m]\n"  /* Red "FAIL"        */
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"     /* Cyan keyword      */
#define OFMT_DBL "%.10lg"             /* Output format for double parameters */

/*============================================================================*\
                          Definitions for error codes
\*============================================================================*/
#define BAOFLIT_ERR_MEMORY      (-1)
#define BAOFLIT_ERR_ARG         (-2)
#define BAOFLIT_ERR_FILE        (-3)
#define BAOFLIT_ERR_CFG         (-4)
#define BAOFLIT_ERR_INIT        (-5)
#define BAOFLIT_ERR_SAVE        (-15)
#define BAOFLIT_ERR_UNKNOWN     (-99)

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/
#define P_ERR(...) fprintf(stderr, FMT_ERR " " __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN " " __VA_ARGS__)
#define P_EXT(...) fprintf(stderr, FMT_EXIT " " __VA_ARGS__)

/* k should vary the fastest to reduce cache miss. */
#define IDX(Ng,i,j,k)      (((size_t) (i) * (Ng) + (j)) * (Ng) + (k))

#endif

