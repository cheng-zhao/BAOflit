/*******************************************************************************
* proc_pk.h: this file is part of the BAOflit program.

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

#ifndef __PROC_PK_H__
#define __PROC_PK_H__

#include <stdio.h>

/*============================================================================*\
                 Functions for pre-processing the power spectra
\*============================================================================*/

/******************************************************************************
Function `pk_interp`:
  Interpolate the power spectrum at given abscissas (k values).
  Note that the original power spectrum is modified.
Arguments:
  * `k`:        abscissas of the input power spectrum in ascending order;
  * `P`:        power of the input power spectrum;
  * `num`:      number of data points of the input power spectrum;
  * `lnkout`:   ln(k) for the power spectrum to be interpolated at;
  * `Pout`:     array for storing the interpolated power spectrum;
  * `nout`:     number of data points for the interpolated power spectrum.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int pk_interp(double *k, double *P, const size_t num, const double *lnkout,
    double *Pout, const size_t nout);

/******************************************************************************
Function `pk_norm`:
  Normalise the linear non-wiggle power spectrum based on the one with wiggles.
Arguments:
  * `k`:        k of the linear power spectrum with BAO wiggles;
  * `P`:        linear power spectrum with BAO wiggles;
  * `n`:        number of points for the linear power spectrum with wiggles;
  * `knw`:      k of the linear non-wiggle power spectrum;
  * `Pnw`:      linear non-wiggle power spectrum;
  * `nnw`:      number of points for the linear non-wiggle power spectrum;
  * `knorm`:    maximum k value used for the normalisation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int pk_norm(const double *k, const double *P, const size_t n, const double *knw,
    double *Pnw, const size_t nnw, const double knorm);

/******************************************************************************
Function `pk_nw_EH`:
  Compute the linear non-wiggle power spectrum based on the fitting formulae of
  Eisenstein & Hu (1998).
Arguments:
  * `k`:        k values at which the power spectrum is evaluated;
  * `n`:        number of data points for the power spectrum;
  * `h`:        the dimensionless hubble parameter;
  * `Om`:       density parameter of matter at redshift 0;
  * `Ob`:       density parameter of baryons at redshift 0;
  * `Tcmb`:     the temperature of cosmic microwave background;
  * `ns`:       scalar index of the initial power spectrum;
  * `Pnw`:      the resulting linear non-wiggle power spectrum.
******************************************************************************/
void pk_nw_EH(const double *k, const size_t n, const double h, const double Om,
    const double Ob, const double Tcmb, const double ns, double *Pnw);

#endif
