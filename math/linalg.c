/*******************************************************************************
* linalg.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "linalg.h"
#include <stdlib.h>
#include <math.h>

/* Shortcut for matrix indexing. */
#define MATRIX(A,m,i,j)   (A)[(i) + (j) * (m)]

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
double *qrdcmp(double *A, const size_t m, const size_t n) {
  if (!m || !n || m < n) {
    P_ERR("invalid dimension of the matrix: %zu x %zu\n", m, n);
    return NULL;
  }
  double *R = malloc(sizeof(double) * ((n * (n + 1)) >> 1));
  if (!R) {
    P_ERR("failed to allocate memory for the upper triangular matrix\n");
    return NULL;
  }

  for (size_t j = 0; j < n - 1; j++) {
    double sum = 0;
    for (size_t i = j; i < m; i++) sum += MATRIX(A,m,i,j) * MATRIX(A,m,i,j);
    double d = (MATRIX(A,m,j,j) > 0) ? sqrt(sum) : -sqrt(sum);
    MATRIX(A,m,j,j) += d;
    /* The sign of `d` is different from that of the reference, and the sqrt
       of `fak` is not evaluated, for efficiency consideration. */
    double fak = d * MATRIX(A,m,j,j);
    for (size_t i = j + 1; i < n; i++) {
      sum = 0;
      for (size_t k = j; k < m; k++) sum += MATRIX(A,m,k,j) * MATRIX(A,m,k,i);
      /* In the reference `fak` is divided twice through A_{kj}.
         Here it is divided only once as sqrt is omitted. */
      sum /= fak;
      for (size_t k = j; k < m; k++) MATRIX(A,m,k,i) -= MATRIX(A,m,k,j) * sum;
    }
    MATRIX(A,m,j,j) = -d;     /* correct for A_{jj} */
  }
  if (m != n) {       /* deal with the last element */
    double sum = 0;
    for (size_t k = n - 1; k < m; k++)
      sum += MATRIX(A,m,k,n-1) * MATRIX(A,m,k,n-1);
    MATRIX(A,m,n-1,n-1) = (MATRIX(A,m,n-1,n-1) > 0) ? -sqrt(sum) : sqrt(sum);
  }

  size_t l = (n << 1) - 1;
  for (size_t i = 0; i < n; i++) {
    size_t offset = ((l - i) * i) >> 1;
    for (size_t j = i; j < n; j++) R[offset + j] = MATRIX(A,m,i,j);
  }
  return R;
}

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
void fwd_subst(const double *A, const double *b, const size_t n, double *x) {
  for (size_t i = 0; i < n; i++) {
    x[i] = b[i];
    for (size_t j = 0; j < i; j++) {
      size_t offset = (((n << 1) - 1 - j) * j) >> 1;
      x[i] -= A[offset + i] * x[j];
    }
    x[i] /= A[(((n << 1) + 1 - i) * i) >> 1];
  }
}

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
void bwd_subst(const double *A, const double *b, const size_t n, double *x) {
  for (size_t ii = n; ii != 0; ii--) {
    size_t i = ii - 1;
    x[i] = b[i];
    size_t offset = (((n << 1) - 1 - i) * i) >> 1;
    for (size_t j = ii; j < n; j++) x[i] -= A[offset + j] * x[j];
    x[i] /= A[(((n << 1) + 1 - i) * i) >> 1];
  }
}

