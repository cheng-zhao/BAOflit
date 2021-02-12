/*******************************************************************************
* linalg.h: this file is part of the BAOflit program.

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

#ifndef __LINALG_H__
#define __LINALG_H__

#include <stdio.h>

/*******************************************************************************
  Implementation of QR decomposition with the Householder algorithm.
  ref: https://people.inf.ethz.ch/gander/papers/qrneu.pdf
*******************************************************************************/

/******************************************************************************
Function `qrdcmp`:
  Perform QR decomposition of a matrix, and return the upper triangular R.
Arguments:
  * `A`:        the matrix to be decomposed;
  * `m`:        number of rows of the matrix;
  * `n`:        number of columns of the matrix.
Return:
  The decomposed upper triangular matrix R.
******************************************************************************/
double *qrdcmp(double *A, const size_t m, const size_t n);

/******************************************************************************
Function `fwd_subst`:
  Solve Ax = b with forward substitution, given the lower triangular matrix `A`
  and a vector `b`.
Arguments:
  * `A`:        the lower triangular matrix, represented as array;
  * `b`:        the vector on the right hand side;
  * `n`:        length of the vector `b`;
  * `x`:        the solution.
******************************************************************************/
void fwd_subst(const double *A, const double *b, const size_t n, double *x);

/******************************************************************************
Function `bwd_subst`:
  Solve Ax = b with backward substitution, given the upper triangular matrix `A`
  and a vector `b`.
Arguments:
  * `A`:        the lower triangular matrix, represented as array;
  * `b`:        the vector on the right hand side;
  * `n`:        length of the vector `b`;
  * `x`:        the solution.
******************************************************************************/
void bwd_subst(const double *A, const double *b, const size_t n, double *x);

#endif
