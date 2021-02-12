/*******************************************************************************
* cspline.h: this file is part of the BAOflit program.

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

#ifndef __CSPLINE_H__
#define __CSPLINE_H__

#include <stdio.h>

/*******************************************************************************
  Implementation of the "natural" cubic spline interpolation algorithm.
  ref: https://doi.org/10.5281/zenodo.3611922
  see also: https://arxiv.org/abs/2001.09253

  The original source codes are released under a CC0 license by Haysn Hornbeck.
*******************************************************************************/

/*============================================================================*\
                    Interface for cubic spline interpolation
\*============================================================================*/

/******************************************************************************
Function `cspline_ypp`:
  Compute the second derivative of sample points.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `y`:        y coordinates of the sample points;
  * `n`:        number of sample points;
  * `ypp`:      array containing (2 * n) elements, with the first n elements
                being the second derivative of `y`.
******************************************************************************/
void cspline_ypp(const double *x, const double *y, const size_t n,
    double *ypp);

/******************************************************************************
Function `cspline_eval_array`:
  Evaluate the cubic spline interpolation for an array with ascending abscissas.
Arguments:
  * `x`:        x coordinates of the data to be interpolated;
  * `y`:        y coordinates of the data to be interpolated;
  * `ypp`:      second derivative of `y`;
  * `n`:        number of data points to be interpolated;
  * `xv`:       x coordinates of the output array;
  * `yv`:       y coordinates of the output array;
  * `nv`:       number of output data points.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cspline_eval_array(const double *x, const double *y, const double *ypp,
    const size_t n, const double *xv, double *yv, const size_t nv);

#endif

