/*******************************************************************************
* read_data.h: this file is part of the BAOflit program.

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

#ifndef __READ_DATA_H__
#define __READ_DATA_H__

#include <stdio.h>

/*============================================================================*\
                          Functions for reading files
\*============================================================================*/

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
    double **x, double **y, size_t *num);

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
    const size_t nref, size_t *nmock, double **xi);

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
    size_t *nmock);

#endif
