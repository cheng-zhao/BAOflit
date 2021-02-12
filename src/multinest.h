/*******************************************************************************
* multinest.h: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* This file is taken from the MultiNest library:
    https://github.com/farhanferoz/MultiNest
    Copyright (c) 2010 Farhan Feroz & Mike Hobson
  but with minor modifications.

  It was distributed following the MultiNest License Agreement.
  See https://github.com/farhanferoz/MultiNest for details.

*******************************************************************************/

#ifndef __MULTINEST_H__
#define __MULTINEST_H__

#ifdef __INTEL_COMPILER	/* if MultiNest was built with ifort */
  #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__  /* if MultiNest was built with gfortran */
  #define NESTRUN __nested_MOD_nestrun
#else
  #error Check (nm libnest3.a | grep nestrun) and edit multinest.h
#endif

#include "define.h"
#include <string.h>

/*============================================================================*\
                            C Interface to MultiNest
\*============================================================================*/

extern void NESTRUN(int *, int *, int *, int *, double *, double *, int *,
    int *, int *, int *, int *, double *, char *, int *, int *, int *, int *,
    int *, int *, double *, int *, void (*Loglike)(double *, int *, int *,
    double *, void *), void (*dumper)(int *, int *, int *, double **, double **,
    double **, double *, double *, double *, double *, void *), void *context);

void run(int IS, int mmodal, int ceff, int nlive, double tol, double efr,
    int ndims, int nPar, int nClsPar, int maxModes, int updInt, double Ztol,
    char root[], int seed, int *pWrap, int fb, int resume, int outfile,
    int initMPI, double logZero, int maxiter,
    void (*LogLike)(double *, int *, int *, double *, void *),
    void (*dumper)(int *, int *, int *, double **, double **, double **,
        double *, double *, double *, double *, void *),
    void *context) {
  for (size_t i = strlen(root); i < BAOFLIT_MN_FNAME_LEN; i++) root[i] = ' ';
  NESTRUN(&IS, &mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar,
      &maxModes, &updInt, &Ztol, root, &seed, pWrap, &fb, &resume, &outfile,
      &initMPI, &logZero, &maxiter, LogLike, dumper, context);
}

#endif
