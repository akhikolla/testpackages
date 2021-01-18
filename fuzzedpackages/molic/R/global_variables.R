## R CMD check fails if we dont make these globals (due to NSE)
## https://www.r-bloggers.com/no-visible-binding-for-global-variable/
utils::globalVariables(c(
  "z",            # .sim_internal
  ".region",      # plot.outlier
  ".dev",         # plot.outlier
  "Deviance",     # plot.multiple_models
  "response",     # plot.multiple_models
  "..quantile..", # plot.multiple_models
  "x1",           # plot.multiple_models
  "x2",           # plot.multiple_models
  "y1",           # plot.multiple_models
  "y2",           # plot.multiple_models
  "y")            # plot.multiple_models
)
