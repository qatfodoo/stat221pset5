source("kuatefodouop_functions.R")

## Fit EM with the refined model

load("./dat/1rout_y.dat")

theta.est <- smoothed_EM(y)
