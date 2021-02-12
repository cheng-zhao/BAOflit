/*******************************************************************************
* proc_pk.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "cspline.h"
#include "proc_pk.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
    double *Pout, const size_t nout) {
  if (!k || !P || !lnkout || !Pout || !num || !nout) {
    P_ERR("the power spectra are not initialised correctly\n");
    return BAOFLIT_ERR_ARG;
  }

  /* Interpolate in log scale. */
  for (size_t i = 0; i < num; i++) {
    k[i] = log(k[i]);
    P[i] = log(P[i]);
  }

  if (k[0] > lnkout[0] || k[num - 1] < lnkout[nout - 1]) {
    P_ERR("interpolation range (in log scale) [" OFMT_DBL "," OFMT_DBL "] "
        "is outside the k range of the power spectrum: [" OFMT_DBL "," OFMT_DBL
        "]\n", lnkout[0], lnkout[nout - 1], k[0], k[num - 1]);
    return BAOFLIT_ERR_INIT;
  }

  double *ypp = malloc(sizeof(double) * num * 2);
  if (!ypp) {
    P_ERR("failed to allocate memory for interpolation\n");
    return BAOFLIT_ERR_MEMORY;
  }

  cspline_ypp(k, P, num, ypp);
  if (cspline_eval_array(k, P, ypp, num, lnkout, Pout, nout)) {
    P_ERR("failed to interpolate the power spectrum\n");
    return BAOFLIT_ERR_INIT;
  }

  /* Output P rather than ln(P). */
  for (size_t i = 0; i < nout; i++) Pout[i] = exp(Pout[i]);

  free(ypp);
  return 0;
}

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
    double *Pnw, const size_t nnw, const double knorm) {
  /* Check arguments. */
  if (!k || !P || !knw || !Pnw || !n || !nnw) {
    P_ERR("the fitting arguments are not initialised correctly\n");
    return BAOFLIT_ERR_ARG;
  }
  if (k[0] > knorm || knw[0] > knorm) {
    P_ERR("not enough k values below " OFMT_DBL " for normalisation\n", knorm);
    return BAOFLIT_ERR_INIT;
  }

  size_t m, mnw;        /* number of points used for normalisation */
  for (m = 0; m < n; m++) {
    if (k[m] > knorm) break;
  }
  if (m == n - 1) {
    P_ERR("k range of the linear matter power spectrum is too small for "
        "normalisation: [" OFMT_DBL "," OFMT_DBL "]\n", k[0], k[m]);
    return BAOFLIT_ERR_INIT;
  }
  for (mnw = 0; mnw < nnw; mnw++) {
    if (knw[mnw] > knorm) break;
  }
  if (mnw == nnw - 1) {
    P_ERR("k range of the linear non-wiggle matter power spectrum is too "
        "small for normalisation: [" OFMT_DBL "," OFMT_DBL "]\n",
        knw[0], knw[mnw]);
    return BAOFLIT_ERR_INIT;
  }

  /* Check if the k values for the two power spectra are identical. */
  if (m == mnw) {
    bool samex = true;
    for (size_t i = 0; i < m; i++) {
      if (k[i] - knw[i] > DOUBLE_TOL || k[i] - knw[i] < -DOUBLE_TOL) {
        samex = false;
        break;
      }
    }
    if (samex) {
      double sum1, sum2;
      sum1 = sum2 = 0;
      for (size_t i = 0; i < m; i++) sum1 += P[i];
      for (size_t i = 0; i < m; i++) sum2 += Pnw[i];
      sum1 /= sum2;
      for (size_t i = 0; i < nnw; i++) Pnw[i] *= sum1;
      return 0;
    }
  }

  /* Determine the power spectrum to be interpolated. */
  size_t nsp, nv;
  const double *kref, *Pref, *xv;
  if (k[0] < knw[0]) {  /* interpolate the power spectrum with wiggles */
    nsp = n;
    nv = mnw;
    kref = k;
    Pref = P;
    xv = knw;
  }
  else {                /* interpolate the non-wiggle power spectrum */
    nsp = nnw;
    nv = m;
    kref = knw;
    Pref = Pnw;
    xv = k;
  }

  /* Allocate memory for interpolation. */
  double *ksp, *Psp, *kv, *Pv;
  ksp = Psp = kv = Pv = NULL;
  ksp = malloc(sizeof(double) * nsp);
  Psp = malloc(sizeof(double) * nsp);
  kv = malloc(sizeof(double) * nv);
  Pv = malloc(sizeof(double) * nv);
  if (!ksp || !Psp || !kv || !Pv) {
    P_ERR("failed to allocate memory for power spectrum normalisation\n");
    if (ksp) free(ksp);
    if (Psp) free(Psp);
    if (kv) free(kv);
    if (Pv) free(Pv);
    return BAOFLIT_ERR_MEMORY;
  }

  /* Copy the arrays as they are modified during interpolation. */
  memcpy(ksp, kref, sizeof(double) * nsp);
  memcpy(Psp, Pref, sizeof(double) * nsp);
  for (size_t i = 0; i < nv; i++) kv[i] = log(xv[i]);

  if (pk_interp(ksp, Psp, nsp, kv, Pv, nv)) {
    free(ksp); free(Psp); free(kv); free(Pv);
    return BAOFLIT_ERR_INIT;
  }

  double sum = 0;
  for (size_t i = 0; i < nv; i++) sum += Pv[i];
  if (k[0] < knw[0]) {  /* interpolate the power spectrum with wiggles */
    double sum2 = 0;
    for (size_t i = 0; i < nv; i++) sum2 += Pnw[i];
    sum /= sum2;
  }
  else {                /* interpolate the non-wiggle power spectrum */
    double sum2 = 0;
    for (size_t i = 0; i < nv; i++) sum2 += P[i];
    sum = sum2 / sum;
  }
  for (size_t i = 0; i < nnw; i++) Pnw[i] *= sum;

  free(ksp); free(Psp); free(kv); free(Pv);
  return 0;
}

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
    const double Ob, const double Tcmb, const double ns, double *Pnw) {
  double Omh2 = Om * h * h;
  double Obh2 = Ob * h * h;
  double Ofac = Ob / Om;
  double s = 44.5 * log(9.83 / Omh2) / sqrt(1 + 10 * pow(Obh2, 0.75));
  double alpha = 1 - 0.328 * log(431 * Omh2) * Ofac
      + 0.38 * log(22.3 * Omh2) * Ofac * Ofac;
  double T2 = Tcmb / 2.7;
  T2 *= T2;

  for (size_t i = 0; i < n; i++) {
    Pnw[i] = Om * h * (alpha + (1 - alpha) / (1 + pow(0.43 * k[i] * s, 4)));
    Pnw[i] = k[i] * T2 / Pnw[i];
    Pnw[i] = log(2 * M_E + 1.8 * Pnw[i]) / (log(2 * M_E + 1.8 * Pnw[i])
        + (14.2 + 731 / (1 + 62.5 * Pnw[i])) * Pnw[i] * Pnw[i]);
    Pnw[i] *= Pnw[i] * pow(k[i], ns);
  }
}

