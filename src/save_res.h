/*******************************************************************************
* save_res.h: this file is part of the BAOflit program.

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

#ifndef __SAVE_RES_H__
#define __SAVE_RES_H__

#include "load_conf.h"
#include <stdio.h>

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
    const size_t nmock);

/******************************************************************************
Function `save_param`:
  Save the names and prior limits of fitting parameters.
Arguments:
  * `conf`:     structure for storing configuration parameters.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_param(const CONF *conf);

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
    const double *y, const size_t n, const size_t *idx, const int nidx);

#endif
