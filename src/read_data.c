/*******************************************************************************
* read_data.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "read_data.h"
#include "timsort.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

/*============================================================================*\
                      Functions for reading file by chunks
\*============================================================================*/

/******************************************************************************
Function `chunk_resize`:
  Enlarge the size of a chunk.
Arguments:
  * `chunk`:    address of the chunk;
  * `size`:     size of the chunk.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int chunk_resize(char **chunk, size_t *size) {
  /* Assume the arguments are not NULL. */
  size_t num;
  if (!(*chunk)) num = BAOFLIT_FILE_CHUNK;
  else {
    if (BAOFLIT_MAX_CHUNK / 2 < *size) return BAOFLIT_ERR_FILE;
    num = *size << 1;
  }

  char *tmp = realloc(*chunk, num * sizeof(char));
  if (!tmp) return BAOFLIT_ERR_MEMORY;

  *chunk = tmp;
  *size = num;
  return 0;
}

/******************************************************************************
Function `read_table`:
  Read two columns of an ASCII file as double arrays.
Arguments:
  * `fname`:    filename of the input catalog;
  * `comment`:  symbol indicating lines to be skipped;
  * `colx`:     number (starting from 1) of the first column to be read;
  * `coly`:     number (starting from 1) of the second column to be read;
  * `xmin`:     minimum value of the first column to be read;
  * `xmax`:     maximum value of the first column to be read;
  * `x`:        array for the first column;
  * `y`:        array for the second column;
  * `num`:      number of lines read successfully.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_table(const char *fname, const char comment, const int colx,
    const int coly, const double xmin, const double xmax,
    double **x, double **y, size_t *num) {
  /* Check arguments. */
  if (!fname || !(*fname)) {
    P_ERR("the input catalog is not set\n");
    return BAOFLIT_ERR_ARG;
  }
  if (colx <= 0 || coly <= 0 || colx == coly) {
    P_ERR("invalid numbers of columns: %d and %d\n", colx, coly);
    return BAOFLIT_ERR_ARG;
  }
  if (!x || !y || !num) {
    P_ERR("arrays for the table are not initialised\n");
    return BAOFLIT_ERR_ARG;
  }

  /* Construct the formatter string for reading lines. */
  int maxcol = (colx > coly) ? colx : coly;
  char *fmtr = calloc(maxcol * 4, sizeof(char));
  for (int i = 0; i < maxcol; i++) {
    fmtr[i * 4] = '%';
    if (colx - 1 == i || coly - 1 == i) {
      fmtr[i * 4 + 1] = 'l';
      fmtr[i * 4 + 2] = 'f';
    }
    else {
      fmtr[i * 4 + 1] = '*';
      fmtr[i * 4 + 2] = 's';
    }
    if (i == maxcol - 1) fmtr[i * 4 + 3] = '\0';
    else fmtr[i * 4 + 3] = ' ';
  }

  /* Open the file for reading. */
  FILE *fp;
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    free(fmtr);
    return BAOFLIT_ERR_FILE;
  }

  /* Prepare for the chunk. */
  char *chunk = NULL;
  size_t csize = 0;
  if (chunk_resize(&chunk, &csize)) {
    P_ERR("failed to allocate memory for reading the file by chunks\n");
    fclose(fp); free(fmtr);
    return BAOFLIT_ERR_MEMORY;
  }

  /* Allocate memory for the data. */
  size_t max = BAOFLIT_DATA_INIT_NUM;
  double *c1, *c2, *px;
  c1 = c2 = NULL;
  if (!(c1 = malloc(max * sizeof(double))) ||
      !(c2 = malloc(max * sizeof(double)))) {
    P_ERR("failed to allocate memory for the samples\n");
    fclose(fp); free(chunk); free(fmtr);
    if (c1) free(c1);
    if (c2) free(c2);
    return BAOFLIT_ERR_MEMORY;
  }
  px = (colx < coly) ? c1 : c2;         /* pointer to the x column */

  size_t n, nread, nrest;
  n = nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == comment || *p == '\0') {        /* comment or empty */
        p = endl + 1;
        continue;
      }

      /* Parse the line. */
      if (sscanf(p, fmtr, c1 + n, c2 + n) != 2) {
        P_ERR("failed to read line: %s\n", p);
        fclose(fp); free(chunk); free(fmtr); free(c1); free(c2);
        return BAOFLIT_ERR_FILE;
      }

      /* Check the x value to see if saving numbers in this line. */
      if (px[n] >= xmin && px[n] <= xmax) {
        /* Enlarge the memory for the data if necessary. */
        if (++n >= max) {
          if (SIZE_MAX / 2 < max) {
            P_ERR("too many samples in the file: `%s'\n", fname);
            fclose(fp); free(chunk); free(fmtr); free(c1); free(c2);
            return BAOFLIT_ERR_FILE;
          }
          max <<= 1;
          double *tmp = realloc(c1, sizeof(double) * max);
          if (!tmp) {
            P_ERR("failed to allocate memory for the samples\n");
            fclose(fp); free(chunk); free(fmtr); free(c1); free(c2);
            return BAOFLIT_ERR_MEMORY;
          }
          c1 = tmp;
          tmp = realloc(c2, sizeof(double) * max);
          if (!tmp) {
            P_ERR("failed to allocate memory for the samples\n");
            fclose(fp); free(chunk); free(fmtr); free(c1); free(c2);
            return BAOFLIT_ERR_MEMORY;
          }
          c2 = tmp;
          px = (colx < coly) ? c1 : c2;
        }
      }

      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        fclose(fp); free(chunk); free(fmtr); free(c1); free(c2);
        return BAOFLIT_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  free(chunk); free(fmtr);
  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    fclose(fp); free(c1); free(c2);
    return BAOFLIT_ERR_FILE;
  }
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  if (n <= 1) {
    P_ERR("too few elements read from file: `%s'\n", fname);
    free(c1); free(c2);
    return BAOFLIT_ERR_FILE;
  }

  /* Reduce the memory cost if applicable. */
  double *tmp = realloc(c1, sizeof(double) * n);
  if (tmp) c1 = tmp;
  tmp = realloc(c2, sizeof(double) * n);
  if (tmp) c2 = tmp;

  if (colx < coly) {
    *x = c1;
    *y = c2;
  }
  else {
    *x = c2;
    *y = c1;
  }
  /* Sort the two arrays by x. */
  tim_sort(*x, *y, n);

  *num = n;
  return 0;
}

/******************************************************************************
Function `read_mocks`:
  Read 2PCFs of mocks, from files in a list.
Arguments:
  * `fname`:    filename of the input catalog list;
  * `comment`:  symbol indicating lines to be skipped;
  * `scol`:     number (starting from 1) of the separations to be read;
  * `xicol`:    number (starting from 1) of the 2PCFs to be read;
  * `smin`:     minimum separation of interest;
  * `smax`:     maximum separation of interest;
  * `sref`:     separations of the data;
  * `nref`:     number of data points that should be read;
  * `nmock`:    number of mock 2PCFs that are read successfully;
  * `xi`:       array for all the mock 2PCFs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_mocks(const char *fname, const char comment, const int scol,
    const int xicol, const double smin, const double smax, const double *sref,
    const size_t nref, size_t *nmock, double **xi) {
  /* Check arguments. */
  if (!fname || !(*fname)) {
    P_ERR("the list for mock 2PCFs is not set\n");
    return BAOFLIT_ERR_ARG;
  }
  if (!sref || !nref) {
    P_ERR("the reference separations and number of data points are not set\n");
    return BAOFLIT_ERR_ARG;
  }
  if (!nmock || !xi) {
    P_ERR("arrays for the mock 2PCFs are not initialised\n");
    return BAOFLIT_ERR_ARG;
  }

  /* Prepare for the chunk for file reading. */
  char *chunk = NULL;
  size_t csize = 0;
  if (chunk_resize(&chunk, &csize)) {
    P_ERR("failed to allocate memory for reading the mock 2PCFs by chunks\n");
    return BAOFLIT_ERR_MEMORY;
  }

  /* Allocate memory for the mock 2PCFs. */
  size_t max = BAOFLIT_DATA_INIT_NUM;
  double *ximock = malloc(sizeof(double) * nref * max);
  if (!ximock) {
    P_ERR("failed to allocate memory for the mock 2PCFs\n");
    free(chunk);
    return BAOFLIT_ERR_MEMORY;
  }

  /* Open the file for reading. */
  FILE *fp;
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open the file for reading: `%s'\n", fname);
    free(chunk); free(ximock);
    return BAOFLIT_ERR_FILE;
  }

  size_t n, nread, nrest;
  n = nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */
    char *endl;

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == comment || *p == '\0') {        /* comment or empty */
        p = endl + 1;
        continue;
      }

      /* Read the mock 2PCF in this line. */
      double *x, *y;
      x = y = NULL;
      size_t num = 0;
      if (read_table(p, comment, scol, xicol, smin, smax, &x, &y, &num)) {
        fclose(fp); free(chunk); free(ximock);
        if (x) free(x);
        if (y) free(y);
        return BAOFLIT_ERR_FILE;
      }

      /* Check if the separations are consistent with those of data. */
      if (num != nref) {
        P_ERR("unexpected number of data points in file: `%s'\n", p);
        fclose(fp); free(chunk); free(ximock); free(x); free(y);
        return BAOFLIT_ERR_FILE;
      }
      for (size_t i = 0; i < num; i++) {
        if (sref[i] - x[i] < -DOUBLE_TOL || sref[i] - x[i] > DOUBLE_TOL) {
          fclose(fp); free(chunk); free(ximock); free(x); free(y);
          P_ERR("unexpected separation (" OFMT_DBL ") in file: `%s'\n",
              x[i], p);
          return BAOFLIT_ERR_FILE;
        }
      }

      /* Record the mock 2PCF. */
      memcpy(ximock + n * nref, y, sizeof(double) * nref);
      free(x); free(y);

      /* Enlarge the memory for storing all mock 2PCFs if necessary. */
      if (++n >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many samples in the file: `%s'\n", fname);
          fclose(fp); free(chunk); free(ximock);
          return BAOFLIT_ERR_FILE;
        }
        max <<= 1;
        double *tmp = realloc(ximock, sizeof(double) * nref * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the mock 2PCFs\n");
          fclose(fp); free(chunk); free(ximock);
          return BAOFLIT_ERR_MEMORY;
        }
        ximock = tmp;
      }

      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        fclose(fp); free(chunk); free(ximock);
        return BAOFLIT_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  free(chunk);
  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    fclose(fp); free(ximock);
    return BAOFLIT_ERR_FILE;
  }
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  /* Reduce the memory cost if applicable. */
  double *tmp = realloc(ximock, sizeof(double) * nref * n);
  if (tmp) ximock = tmp;

  *xi = ximock;
  *nmock = n;
  return 0;
}

/******************************************************************************
Function `read_Rcov`:
  Read the upper triangular decomposition of the covariance matrix from file.
Arguments:
  * `fname`:    filename for the input covariance matrix;
  * `Rcov`:     upper triangular decomposition of the covariance matrix;
  * `nbin`:     dimension of the covariance matrix;
  * `nmock`:    number of mocks for the covariance matrix estimation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_Rcov(const char *fname, double **Rcov, const size_t nbin,
    size_t *nmock) {
  /* Check arguments. */
  if (!fname || !(*fname)) {
    P_ERR("the input file for the covariance matrix is not set\n");
    return BAOFLIT_ERR_ARG;
  }
  if (!Rcov || !nbin || !nmock) {
    P_ERR("the fitting arguments are not initialised correctly\n");
    return BAOFLIT_ERR_ARG;
  }

  FILE *fp = fopen(fname, "r");
  if (!fp) {
    P_ERR("cannot open the file for reading: `%s'\n", fname);
    return BAOFLIT_ERR_FILE;
  }

  size_t num = 0;
  if (fread(&num, sizeof(size_t), 1, fp) != 1) {
    P_ERR("failed to read the number of mocks from file: `%s'\n", fname);
    fclose(fp);
    return BAOFLIT_ERR_FILE;
  }
  if (!num) {
    P_ERR("invalid number of mocks read from file: `%s'\n", fname);
    fclose(fp);
    return BAOFLIT_ERR_FILE;
  }
  *nmock = num;

  if (fread(&num, sizeof(size_t), 1, fp) != 1) {
    P_ERR("failed to read the dimension of covariance matrix from file: `%s'\n",
        fname);
    fclose(fp);
    return BAOFLIT_ERR_FILE;
  }
  if (num != nbin) {
    P_ERR("unexpected dimension of covariance matrix: %zu rather than %zu\n",
        num, nbin);
    fclose(fp);
    return BAOFLIT_ERR_FILE;
  }

  size_t lenR = (num * (num + 1)) / 2;
  double *cov = malloc(sizeof(double) * lenR);
  if (!cov) {
    P_ERR("failed to allocate memory for the covariance matrix\n");
    fclose(fp);
    return BAOFLIT_ERR_MEMORY;
  }
  if (fread(cov, sizeof(double) * lenR, 1, fp) != 1) {
    P_ERR("failed to read the covariance matrix from file: `%s'\n", fname);
    fclose(fp); free(cov);
    return BAOFLIT_ERR_FILE;
  }

  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  *Rcov = cov;
  return 0;
}
