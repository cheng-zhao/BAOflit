/*******************************************************************************
* save_res.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "save_res.h"
#include <stdlib.h>
#include <string.h>

/*============================================================================*\
                          Functions for saving results
\*============================================================================*/

/******************************************************************************
Function `save_Rcov`:
  Save the upper triangular decomposition of the covariance matrix to file.
Arguments:
  * `fname`:    filename for the output covariance matrix;
  * `Rcov`:     upper triangular decomposition of the covariance matrix;
  * `nbin`:     dimension of the covariance matrix;
  * `nmock`:    number of mocks for the covariance matrix estimation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_Rcov(const char *fname, const double *Rcov, const size_t nbin,
    const size_t nmock) {
  /* Check arguments. */
  if (!fname || !(*fname)) {
    P_ERR("the output filename is not set\n");
    return BAOFLIT_ERR_ARG;
  }
  if (!Rcov || !nbin || !nmock) {
    P_ERR("the covariance matrix is not initialised correctly\n");
    return BAOFLIT_ERR_ARG;
  }

  FILE *fp = fopen(fname, "w");
  if (!fp) {
    P_ERR("cannot open the file for writing: `%s'\n", fname);
    return BAOFLIT_ERR_FILE;
  }

  if (fwrite(&nmock, sizeof(size_t), 1, fp) != 1) {
    P_ERR("failed to write to file: `%s'\n", fname);
    fclose(fp);
    return BAOFLIT_ERR_SAVE;
  }

  if (fwrite(&nbin, sizeof(size_t), 1, fp) != 1) {
    P_ERR("failed to write to file: `%s'\n", fname);
    fclose(fp);
    return BAOFLIT_ERR_SAVE;
  }

  size_t lenR = (nbin * (nbin + 1)) / 2;
  if (fwrite(Rcov, sizeof(double) * lenR, 1, fp) != 1) {
    P_ERR("failed to write to file: `%s'\n", fname);
    fclose(fp);
    return BAOFLIT_ERR_SAVE;
  }

  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);
  return 0;
}

/******************************************************************************
Function `save_param`:
  Save the names and prior limits of fitting parameters.
Arguments:
  * `conf`:     structure for storing configuration parameters.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_param(const CONF *conf) {
  char *fname = calloc(BAOFLIT_MN_FNAME_LEN, sizeof(char));
  if (!fname) {
    P_ERR("failed to allocate memory for the filename\n");
    return BAOFLIT_ERR_MEMORY;
  }
  memcpy(fname, conf->oroot, BAOFLIT_MN_FNAME_LEN);

  /* Save the parameter names. */
  size_t idx = strlen(fname);
  const char pname[] = ".paramnames";
  if (idx + sizeof(pname) >= BAOFLIT_MN_FNAME_LEN) {
    P_ERR("the basename of output files is too long\n");
    return BAOFLIT_ERR_ARG;
  }
  memcpy(fname + idx, pname, sizeof(pname));

  FILE *fp = fopen(fname, "w");
  if (!fp) {
    P_ERR("cannot open file for writing: `%s'\n", fname);
    return BAOFLIT_ERR_FILE;
  }
  fprintf(fp, "alpha $\\alpha$\n");
  for (int i = 0; i < conf->num_B; i++)
    fprintf(fp, "B_%c $B_{\\rm %c}$\n", conf->Bfit[i], conf->Bfit[i]);
  if (!conf->val_Snl) {
    for (int i = 0; i < conf->ninput; i++) {
      fprintf(fp, "Snl_%c%c $\\Sigma_{\\rm nl, %c%c}$\n", conf->tracer[i][0],
          conf->tracer[i][1], conf->tracer[i][0], conf->tracer[i][1]);
    }
  }
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  /* Save the prior limits. */
  memset(fname + idx, 0, sizeof(pname));
  const char rname[] = ".ranges";
  if (idx + sizeof(rname) >= BAOFLIT_MN_FNAME_LEN) {
    P_ERR("the basename of the output files is too long\n");
    return BAOFLIT_ERR_ARG;
  }
  memcpy(fname + idx, rname, sizeof(rname));

  if (!(fp = fopen(fname, "w"))) {
    P_ERR("cannot open file for writing: `%s'\n", fname);
    return BAOFLIT_ERR_FILE;
  }
  fprintf(fp, "alpha " OFMT_DBL " " OFMT_DBL "\n", conf->pmin_a, conf->pmax_a);
  for (int i = 0; i < conf->num_B; i++) {
    fprintf(fp, "B_%c " OFMT_DBL " " OFMT_DBL "\n", conf->Bfit[i],
        conf->pmin_B[i], conf->pmax_B[i]);
  }
  if (!conf->val_Snl) {
    for (int i = 0; i < conf->ninput; i++) {
      fprintf(fp, "Snl_%c%c " OFMT_DBL " " OFMT_DBL "\n", conf->tracer[i][0],
          conf->tracer[i][1], conf->pmin_Snl[i], conf->pmax_Snl[i]);
    }
  }
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  free(fname);
  return 0;
}

/******************************************************************************
Function `save_table`:
  Save a table to a text file.
Arguments:
  * `bname`:    basename of the output file;
  * `suffix`:   suffix of the output filename;
  * `x`:        the first column of the table;
  * `y`:        the second column of the table;
  * `n` :       number of rows to be saved;
  * `idx`:      starting indices for segments (different 2PCFs);
  * `nidx`:     number of segements.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_table(char *bname, const char *suffix, const double *x,
    const double *y, const size_t n, const size_t *idx, const int nidx) {
  if (!bname || !(*bname) || !suffix || !(*suffix) || !x || !y || !n) {
    P_ERR("the table to be saved is not initialised\n");
    return BAOFLIT_ERR_ARG;
  }

  size_t slen = strlen(suffix);
  if (slen >= BAOFLIT_MN_SUFF_MAX) {
    P_ERR("suffix of the output file is too long: `%s'\n", suffix);
    return BAOFLIT_ERR_ARG;
  }
  size_t blen = strlen(bname);
  if (blen + slen >= BAOFLIT_MN_FNAME_LEN) {
    P_ERR("the output file name is too long: `%s%s'\n", bname, suffix);
    return BAOFLIT_ERR_ARG;
  }
  memcpy(bname + blen, suffix, slen);
  bname[blen + slen] = '\0';

  FILE *fp = fopen(bname, "w");
  if (!fp) {
    P_ERR("failed to open file for writing: `%s'\n", bname);
    return BAOFLIT_ERR_FILE;
  }

  if (nidx > 1) {
    fprintf(fp, "%c Starting indices for different 2PCFs:\n%c",
        BAOFLIT_SAVE_COMMENT, BAOFLIT_SAVE_COMMENT);
    for (int i = 0; i < nidx; i++) fprintf(fp, " %zu", idx[i]);
    fprintf(fp, "\n");
  }

  for (size_t i = 0; i < n; i++)
    fprintf(fp, OFMT_DBL " " OFMT_DBL "\n", x[i], y[i]);

  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", bname);

  bname[blen] = '\0';
  return 0;
}
