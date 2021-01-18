library(Countr)

library(rbenchmark)

source("mcShaneCode.R")
data(fertility)

## config: choose parameters (selection process not shown here)
children <- fertility$children
shape <- 1.116
scale <- rep(2.635, length(children))
rep <- 1000
nstepsConv <- c(132, 24, 132, 24, 132, 36)
ntermsSeries <- c(20, 17)
conv_series_acc <- 1e-7

## performance model
perf <- benchmark(direct0 =
                      dWeibullCount_loglik(children, shape, scale, "conv_direct",
                                           1, TRUE, nstepsConv[1],
                                           conv_extrap = FALSE),
                  direct1 =
                      dWeibullCount_loglik(children, shape, scale, "conv_direct",
                                           1, TRUE, nstepsConv[2],
                                           conv_extrap = TRUE),
                  naive0 =
                      dWeibullCount_loglik(children, shape, scale, "conv_naive",
                                           1, TRUE, nstepsConv[3],
                                           conv_extrap = FALSE),
                  naive1 = dWeibullCount_loglik(children, shape, scale,
                                                "conv_naive",
                                                1, TRUE, nstepsConv[4],
                                                conv_extrap = TRUE),
                  dePril0 = dWeibullCount_loglik(children, shape, scale,
                                                 "conv_dePril",
                                                 1, TRUE, nstepsConv[5],
                                                 conv_extrap = FALSE),
                  dePril1 = dWeibullCount_loglik(children, shape, scale,
                                                 "conv_dePril",
                                                 1, TRUE, nstepsConv[6],
                                                 conv_extrap = TRUE),
                  series_mat =
                      dWeibullCount_loglik(children, shape, scale,
                                           "series_mat", 1, TRUE,
                                           series_terms = ntermsSeries[1]),
                  series_acc =
                      dWeibullCount_loglik(children, shape, scale,
                                           "series_acc", 1, TRUE,
                                           series_terms = ntermsSeries[2],
                                           series_acc_eps = conv_series_acc),
                  mcShane = dWeibullCount_McShane(scale, shape,
                                                  children, jmax = 150),
                  replications = rep, order = "relative",
                  columns = c("test", "replications", "relative", "elapsed")
                  )

print(perf)

save.image()
