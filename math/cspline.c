/*******************************************************************************
* cspline.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "cspline.h"
#include <stdlib.h>

/*******************************************************************************
  Implementation of the "natural" cubic spline interpolation algorithm.
  ref: https://doi.org/10.5281/zenodo.3611922
  see also: https://arxiv.org/abs/2001.09253

  The original source codes are released under a CC0 license by Haysn Hornbeck.
*******************************************************************************/

/* Error code */
#define CSPLINE_FAILURE         (-1)

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
    double *ypp) {
  if (n < 2 || !ypp) return;
  double *cp = ypp + n;

  double newx = x[1];
  double newy = y[1];
  double c = x[1] - x[0];
  double newd = (y[1] - y[0]) / c;

  /* natural condition: second derivative = 0 */
  cp[0] = cp[n - 1] = ypp[0] = ypp[n - 1] = 0;

  /* forward substitution */
  size_t j = 1;
  while (j < n - 1) {
    double oldx = newx;
    double oldy = newy;
    double a = c;
    double oldd = newd;

    newx = x[j + 1];
    newy = y[j + 1];
    c = newx - oldx;
    newd = (newy - oldy) / c;

    double b = (c + a) * 2;
    double invd = 1 / (b - a * cp[j - 1]);
    double d = (newd - oldd) * 6;

    ypp[j] = (d - a * ypp[j - 1]) * invd;
    cp[j] = c * invd;

    j += 1;
  }

  /* backward substitution */
  while (j) {
    j -= 1;
    ypp[j] -= cp[j] * ypp[j + 1];
  }
}

/******************************************************************************
Function `cspline_eval`:
  Evaluate the cubic spline interpolation.
  It can be further optimised for uniformly spaced sample points.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `y`:        y coordinates of the sample points;
  * `ypp`:      second derivative of `y`;
  * `xv`:       x coordinate of the value to be evaluated;
  * `i`:        index of the value to be evaluated.
Return:
  The value for the given x coordinate.
******************************************************************************/
static double cspline_eval(const double *x, const double *y, const double *ypp,
    const double xv, const size_t i) {
  size_t j = i + 1;
  double ba = x[j] - x[i];
  double xa = xv - x[i];
  double bx = x[j] - xv;
  double ba2 = ba * ba;

  double lower = xa * y[j] + bx * y[i];
  double c = (xa * xa - ba2) * xa * ypp[j];
  double d = (bx * bx - ba2) * bx * ypp[i];

  /* 1/6 = 0x1.5555555555555p-3 */
  return (lower + 0x1.5555555555555p-3 * (c + d)) / ba;
}

/******************************************************************************
Function `bin_search`:
  Binary search the x coordinate for interpolation.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `n`:        number of the sample points;
  * `xv`:       x coordinate of the value to be evaluated;
  * `istart`:   starting index for the search;
  * `iend`:     ending index for the search.
Return:
  Index of the value in the sample to be evaluated.
******************************************************************************/
static inline size_t bin_search(const double *x, const double xv,
    const size_t istart, const size_t iend) {
  size_t l = istart;
  size_t u = iend;
  while (l <= u) {
    size_t i = (l + u) >> 1;
    if (x[i + 1] <= xv) l = i + 1;
    else if (x[i] > xv) u = i - 1;
    else return i;
  }
  return SIZE_MAX;
}

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
    const size_t n, const double *xv, double *yv, const size_t nv) {
  if (!nv) return 0;
  /* Process the last point first, to reduce the searching range. */
  size_t end = bin_search(x, xv[nv - 1], 0, n - 1);
  if (end >= n - 1) {
    if (xv[nv - 1] == x[n - 1]) yv[nv - 1] = y[n - 1];
    else return CSPLINE_FAILURE;
  }
  else yv[nv - 1] = cspline_eval(x, y, ypp, xv[nv - 1], end);
  if (nv == 1) return 0;

  /* Process the rest of the points. */
  if (end < n - 1) end += 1;
  size_t pos = 0;
  for (size_t i = 0; i < nv - 1; i++) {
    pos = bin_search(x, xv[i], pos, end);
    if (pos >= n - 1) {
      if (xv[i] == x[n - 1]) yv[i] = y[n - 1];
      else return CSPLINE_FAILURE;
    }
    else yv[i] = cspline_eval(x, y, ypp, xv[i], pos);
  }
  return 0;
}
