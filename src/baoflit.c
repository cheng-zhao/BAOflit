/*******************************************************************************
* baoflit.c: this file is part of the BAOflit program.

* BAOflit: Baryon Acoustic Oscillation Fitter for muLtI-Tracers.

* Github repository:
        https://github.com/cheng-zhao/BAOflit

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "fit_args.h"
#include "fit_func.h"

int main(int argc, char *argv[]) {
  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return BAOFLIT_ERR_CFG;
  }

  ARGS *args;
  if (!(args = init_fit(conf))) {
    printf(FMT_FAIL);
    P_EXT("failed to initialise the fitting process\n");
    return BAOFLIT_ERR_INIT;
  }

  /* Run the MultiNest fit. */
  run_multinest(conf, args);

  conf_destroy(conf);
  args_destroy(args);
  return 0;
}
