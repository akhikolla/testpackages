## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  error = TRUE,
  comment = "#>"
)

## ----bm-fun, eval=TRUE--------------------------------------------------------
library("fddm")
library("rtdists")
library("RWiener")
source(system.file("extdata", "Gondan_et_al_density.R", package = "fddm", mustWork = TRUE))
library("microbenchmark")

rt_benchmark_vec <- function(RT, resp, V, A, t0 = 1e-4, W = 0.5, SV = 0.0,
                             err_tol = 1e-6, times = 100, unit = "ns") {

  fnames <- c("fs_SWSE_17", "fs_SWSE_14", "fb_SWSE_17", "fb_SWSE_14",
              "fs_Gon_17", "fs_Gon_14", "fb_Gon_17", "fb_Gon_14",
              "fs_Nav_17", "fs_Nav_14", "fb_Nav_17", "fb_Nav_14",
              "fl_Nav_09", "RWiener", "Gondan", "rtdists")
  nf <- length(fnames) # number of functions being benchmarked
  nV <- length(V)
  nA <- length(A)
  nW <- length(W)
  nSV <- length(SV)
  resp <- rep(resp, length(RT)) # for RWiener

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol = 4+nf, nrow = nV*nA*nW*nSV))
  colnames(mbm_res) <- c('V', 'A', 'W', 'SV', fnames)
  row_idx <- 1

  # Loop through each combination of parameters and record microbenchmark results
  for (v in 1:nV) {
    for (a in 1:nA) {
      for (w in 1:nW) {
        for (sv in 1:nSV) {
          mbm <- microbenchmark(
          fs_SWSE_17 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2017", scale = "small",
                             err_tol = err_tol),
          fs_SWSE_14 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2014", scale = "small",
                             err_tol = err_tol),
          fb_SWSE_17 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2017", scale = "both",
                             err_tol = err_tol),
          fb_SWSE_14 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2014", scale = "both",
                             err_tol = err_tol),
          fs_Gon_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2017", scale = "small",
                            err_tol = err_tol),
          fs_Gon_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2014", scale = "small",
                            err_tol = err_tol),
          fb_Gon_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2017", scale = "both",
                            err_tol = err_tol),
          fb_Gon_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2014", scale = "both",
                            err_tol = err_tol),
          fs_Nav_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2017", scale = "small",
                            err_tol = err_tol),
          fs_Nav_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2014", scale = "small",
                            err_tol = err_tol),
          fb_Nav_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2017", scale = "both",
                            err_tol = err_tol),
          fb_Nav_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2014", scale = "both",
                            err_tol = err_tol),
          fl_Nav_09 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            scale = "large", err_tol = err_tol),
          RWiener = dwiener(RT, resp = resp, alpha = A[a],
                            delta = V[v], tau = t0, beta = W[w],
                            give_log = FALSE),
          Gondan = fs(t = RT-t0, a = A[a], v = V[v],
                      w = W[w], eps = err_tol), # only "lower" resp
          rtdists = ddiffusion(RT, resp, a = A[a], v = V[v],
                               t0 = t0, z = W[w]*A[a]),
          times = times, unit = unit)
        # add the v, a, and w values to the dataframe
        mbm_res[row_idx, 1] <- V[v]
        mbm_res[row_idx, 2] <- A[a]
        mbm_res[row_idx, 3] <- W[w]
        mbm_res[row_idx, 4] <- SV[sv]
        # add the median microbenchmark results to the dataframe
        for (i in 1:nf) {
          mbm_res[row_idx, 4+i] <- median(mbm[mbm[,1] == fnames[i],2])
        }
        # iterate start value
        row_idx = row_idx + 1
        }
      }
    }
  }
  return(mbm_res)
}

rt_benchmark_ind <- function(RT, resp, V, A, t0 = 1e-4, W = 0.5, SV = 0.0,
                             err_tol = 1e-6, times = 100, unit = "ns") {
  fnames <- c("fs_SWSE_17", "fs_SWSE_14", "fb_SWSE_17", "fb_SWSE_14",
              "fs_Gon_17", "fs_Gon_14", "fb_Gon_17", "fb_Gon_14",
              "fs_Nav_17", "fs_Nav_14", "fb_Nav_17", "fb_Nav_14",
              "fl_Nav_09", "RWiener", "Gondan", "rtdists")
  nf <- length(fnames) # number of functions being benchmarked
  nRT <- length(RT)
  nV <- length(V)
  nA <- length(A)
  nW <- length(W)
  nSV <- length(SV)

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol = 5+nf, nrow = nRT*nV*nA*nW*nSV))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'SV', fnames)
  row_idx <- 1

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          for (sv in 1:nSV) {
            mbm <- microbenchmark(
            fs_SWSE_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_SWSE_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fb_SWSE_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_SWSE_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            fs_Gon_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_Gon_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fb_Gon_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_Gon_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            fs_Nav_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_Nav_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fb_Nav_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_Nav_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            fl_Nav_09 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              scale = "large", err_tol = err_tol),
            RWiener = dwiener(RT[rt], resp = resp, alpha = A[a],
                              delta = V[v], tau = t0, beta = W[w],
                              give_log = FALSE),
            Gondan = fs(t = RT[rt]-t0, a = A[a], v = V[v],
                        w = W[w], eps = err_tol), # only "lower" resp
            rtdists = ddiffusion(RT[rt], resp, a = A[a], v = V[v],
                                 t0 = t0, z = W[w]*A[a]),
            times = times, unit = unit)
          # add the v, a, and w values to the dataframe
          mbm_res[row_idx, 1] <- RT[rt]
          mbm_res[row_idx, 2] <- V[v]
          mbm_res[row_idx, 3] <- A[a]
          mbm_res[row_idx, 4] <- W[w]
          mbm_res[row_idx, 5] <- SV[sv]
          # add the median microbenchmark results to the dataframe
          for (i in 1:nf) {
            mbm_res[row_idx, 5+i] <- median(mbm[mbm[,1] == fnames[i],2])
          }
          # iterate start value
          row_idx = row_idx + 1
          }
        }
      }
    }
  }
  return(mbm_res)
}

## ----bm-run, eval=FALSE-------------------------------------------------------
#  # Define parameter space
#  RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
#  A <- c(0.25, 0.5, 1, 2.5, 5)
#  V <- c(-5, -2, 0, 2, 5)
#  t0 <- 1e-4 # must be nonzero for RWiener
#  W <- c(0.2, 0.5, 0.8)
#  SV <- c(0, 0.5, 1, 1.5)
#  err_tol <- 1e-6 # this is the setting from rtdists
#  
#  # Run benchmark tests
#  bm_vec <- rt_benchmark_vec(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
#                             W = W, SV = SV, err_tol = err_tol,
#                             times = 1000, unit = "ns")
#  bm_ind <- rt_benchmark_ind(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
#                             W = W, SV = SV, err_tol = err_tol,
#                             times = 100, unit = "ns")

## ----bm-run-and-save, eval=FALSE, include=FALSE-------------------------------
#  # Define parameter space
#  RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
#  A <- c(0.25, 0.5, 1, 2.5, 5)
#  V <- c(-5, -2, 0, 2, 5)
#  t0 <- 1e-4 # must be nonzero for RWiener
#  W <- c(0.2, 0.5, 0.8)
#  SV <- c(0, 0.5, 1, 1.5)
#  err_tol <- 1e-6 # this is the setting from rtdists
#  
#  # Run benchmark tests
#  bm_vec <- rt_benchmark_vec(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
#                             W = W, SV = SV, err_tol = err_tol,
#                             times = 1000, unit = "ns")
#  save(bm_vec, compress = "xz", compression_level = 9,
#       file = "inst/extdata/bm_vec.Rds")
#  # load(system.file("extdata", "bm_vec.Rds", package = "fddm", mustWork = TRUE))
#  bm_ind <- rt_benchmark_ind(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
#                             W = W, SV = SV, err_tol = err_tol,
#                             times = 100, unit = "ns")
#  save(bm_ind, compress = "xz", compression_level = 9,
#       file = "inst/extdata/bm_ind.Rds")
#  # load(system.file("extdata", "bm_ind.Rds", package = "fddm", mustWork = TRUE))

## ----bm-violin, eval=TRUE, fig.height=5---------------------------------------
library("reshape2")
library("ggplot2")

# load data, will be in the variable 'bm_vec'
load(system.file("extdata", "bm_vec.Rds", package = "fddm", mustWork = TRUE))

t_idx <- match("SV", colnames(bm_vec))
bm_vec[, -seq_len(t_idx)] <- bm_vec[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_vec <- melt(bm_vec, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_vec <- c("fb_SWSE_17", "fb_SWSE_14", "fb_Gon_17", "fb_Gon_14",
               "fb_Nav_17", "fb_Nav_14", "fs_SWSE_17", "fs_SWSE_14",
               "fs_Gon_17", "fs_Gon_14", "fs_Nav_17", "fs_Nav_14",
               "fl_Nav_09", "RWiener", "Gondan", "rtdists")
Color_vec <- c("#e000b4", "#ff99eb", "#e68a00", "#ffb366",
               "#006699", "#66ccff", "#9900cc", "#cc99ff",
               "#c2a500", "#d7db42", "#336600", "#33cc33",
               "#996633", "#ff9999", "#ff5050", "#990000")

mi <- min(bm_vec[, -seq_len(t_idx)])
ma <- max(bm_vec[, (t_idx+1):(ncol(bm_vec)-4)])

ggplot(mbm_vec, aes(x = factor(FuncName, levels = Names_vec), y = time,
                    color = factor(FuncName, levels = Names_vec),
                    fill = factor(FuncName, levels = Names_vec))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = Color_vec, guide = FALSE) +
       scale_fill_manual(values = Color_vec, guide = FALSE) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       coord_cartesian(ylim = c(mi, ma)) +
       labs(title = "Distribution of median benchmark times",
            subtitle = "Dashed lines represent mean benchmark times",
            x = "Method", y = "Time (microseconds)",
            color = "Method") +
       theme_bw() +
       theme(panel.border = element_blank(),
             plot.title = element_text(size = 23),
             plot.subtitle = element_text(size = 16),
             axis.text.x = element_text(size = 16, angle = 90,
                                        vjust = 0.5, hjust = 1),
             axis.text.y = element_text(size = 16),
             axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             legend.position = "none")

## ----bm-meq-prep, eval=TRUE---------------------------------------------------
# load data, will be in the variable 'bm_ind'
load(system.file("extdata", "bm_ind.Rds", package = "fddm", mustWork = TRUE))
bm_ind[["RTAA"]] <- bm_ind[["RT"]] / bm_ind[["A"]] / bm_ind[["A"]]
bm_ind <- bm_ind[, c(1, 2, ncol(bm_ind), 3:(ncol(bm_ind)-1)) ]

t_idx <- match("SV", colnames(bm_ind))
bm_ind[,-seq_len(t_idx)] <- bm_ind[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_ind <- melt(bm_ind, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_meq <- c("fb_SWSE_17", "fs_SWSE_14", "fl_Nav_09",
               "RWiener", "Gondan", "rtdists")
Color_meq <- c("#e000b4", "#cc99ff", "#996633",
               "#ff9999", "#ff5050", "#990000")
mbm_meq <- subset(mbm_ind, FuncName %in% Names_meq)

## ----bm-meq-rtaa, eval=TRUE---------------------------------------------------
ggplot(mbm_meq, aes(x = RTAA, y = time,
                    color = factor(FuncName, levels = Names_meq),
                    fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_log10() +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = bquote(frac(rt, a^2) ~ ", effective response time, " ~ log[10]),
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")

## ----bm-meq-w, eval=TRUE------------------------------------------------------
ggplot(mbm_meq, aes(x = W, y = time,
                    color = factor(FuncName, levels = Names_meq),
                    fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "w, relative starting point (a priori bias)",
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")

## ----bm-meq-v, eval=TRUE------------------------------------------------------
ggplot(mbm_meq, aes(x = V, y = time,
                    color = factor(FuncName, levels = Names_meq),
                    fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "v, drift rate",
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")

## ----bm-meq-sv, eval=TRUE-----------------------------------------------------
ggplot(mbm_meq, aes(x = SV, y = time,
                    color = factor(FuncName, levels = Names_meq),
                    fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "sv, inter-trial variability in the drift rate",
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")

## ----fit-pkg, eval=FALSE------------------------------------------------------
#  library("fddm")
#  library("rtdists")
#  library("microbenchmark")

## ----fit-loglik-fun, eval=TRUE------------------------------------------------
ll_fb_SWSE_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_Gon_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "Gondan", summation_small = "2017",
                scale = "both", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fs_SWSE_14 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "small", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fs_Gon_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "Gondan", summation_small = "2017",
                scale = "small", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fl_Nav_09 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "Navarro", scale = "large", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_RTDists <- function(pars, rt, resp, truth) {
  rtu <- rt[truth == "upper"]
  rtl <- rt[truth == "lower"]
  respu <- resp[truth == "upper"]
  respl <- resp[truth == "lower"]

  densu <- ddiffusion(rtu, respu, a = pars[[3]], v = pars[[1]],
                      z = pars[[5]]*pars[[3]], t0 = pars[[4]], sv = pars[[6]])
  densl <- ddiffusion(rtl, respl, a = pars[[3]], v = pars[[2]],
                      z = pars[[5]]*pars[[3]], t0 = pars[[4]], sv = pars[[6]])

  densities <- c(densu, densl)
  if (any(densities <= 0)) return(1e6)
  return(-sum(log(densities)))
}

## ----fit-fun, eval=TRUE-------------------------------------------------------
rt_fit <- function(data, id_idx = NULL, rt_idx = NULL, response_idx = NULL,
                   truth_idx = NULL, response_upper = NULL, err_tol = 1e-6,
                   times = 100, unit = "ns") {

  # Format data for fitting
  if (all(is.null(id_idx), is.null(rt_idx), is.null(response_idx),
      is.null(truth_idx), is.null(response_upper))) {
    df <- data # assume input data is already formatted
  } else {
    if(any(data[,rt_idx] < 0)) {
      stop("Input data contains negative response times; fit will not be run.")
    }
    if(any(is.na(data[,response_idx]))) {
      stop("Input data contains invalid responses (NA); fit will not be run.")
    }

    nr <- nrow(data)
    df <- data.frame(id = character(nr),
                     rt = double(nr),
                     response = character(nr),
                     truth = character(nr),
                     stringsAsFactors = FALSE)

    if (!is.null(id_idx)) { # relabel identification tags
      for (i in 1:length(id_idx)) {
        idi <- unique(data[,id_idx[i]])
        for (j in 1:length(idi)) {
          df[["id"]][data[,id_idx[i]] == idi[j]] <- paste(
            df[["id"]][data[,id_idx[i]] == idi[j]], idi[j], sep = " ")
        }
      }
      df[["id"]] <- trimws(df[["id"]], which = "left")
    }

    df[["rt"]] <- as.double(data[,rt_idx])

    df[["response"]] <- "lower"
    df[["response"]][data[,response_idx] == response_upper] <- "upper"

    df[["truth"]] <- "lower"
    df[["truth"]][data[,truth_idx] == response_upper] <- "upper"
  }

  # Preliminaries
  ids <- unique(df[["id"]])
  nids <- max(length(ids), 1) # if inds is null, there is only one individual

  init_vals <- data.frame(v1 = c( 0,  10, -.5,  0,  0,  0,  0,  0,  0,   0,  0),
                          v0 = c( 0, -10,  .5,  0,  0,  0,  0,  0,  0,   0,  0),
                          a  = c( 1,   1,   1, .5,  5,  1,  1,  1,  1,   1,  1),
                          t0 = c( 0,   0,   0,  0,  0,  0,  0,  0,  0,   0,  0),
                          w  = c(.5,  .5,  .5, .5, .5, .5, .5, .2, .8,  .5, .5),
                          sv = c( 1,   1,   1,  1,  1,  1,  1,  1,  1, .05,  5))
  ninit_vals <- nrow(init_vals)

  algo_names <- c("fb_SWSE_17", "fb_Gon_17", "fs_SWSE_14",
                  "fs_Gon_17", "fl_Nav_09", "rtdists")
  nalgos <- length(algo_names)
  ni <- nalgos*ninit_vals

  # Initilize the result dataframe
  cnames <- c("ID", "Algorithm", "Convergence", "Objective", "Iterations",
              "FuncEvals", "BmTime")
  res <- data.frame(matrix(ncol = length(cnames), nrow = nids*ninit_vals*nalgos))
  colnames(res) <- cnames

  # label the result dataframe
  res[["ID"]] <- rep(ids, each = ni) # label individuals
  res[["Algorithm"]] <- rep(algo_names, each = ninit_vals) # label algorithms

  # Loop through each individual
  for (i in 1:nids) {
    # extract data for id i
    dfi <- df[df[["id"]] == ids[i], ]
    rti <- dfi[["rt"]]
    respi <- dfi[["response"]]
    truthi <- dfi[["truth"]]

    # starting value for t0 must be smaller than the smallest rt
    min_rti <- min(rti)
    t0_lo <- 0.01*min_rti
    t0_me <- 0.50*min_rti
    t0_hi <- 0.99*min_rti
    init_vals[["t0"]] <- c(rep(t0_me, 5), t0_lo, t0_hi, rep(t0_me, 4))

    # loop through all of the starting values
    for (j in 1:ninit_vals) {
      # get number of evaluations
      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+0*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+0*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+0*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+0*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_Gon_17,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+1*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+1*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+1*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+1*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fs_SWSE_14,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+2*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+2*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+2*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+2*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fs_Gon_17,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+3*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+3*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+3*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+3*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fl_Nav_09,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+4*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+4*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+4*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+4*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_RTDists,
                     rt = rti, resp = respi, truth = truthi,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+5*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+5*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+5*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+5*ninit_vals+j] <- temp[["evaluations"]][[1]]

      # microbenchmark
      mbm <- microbenchmark(
        fb_SWSE_17 = nlminb(init_vals[j,], ll_fb_SWSE_17, err_tol = err_tol,
                            rt = rti, resp = respi, truth = truthi,
                            # limits:   vu,   vl,   a,      t0, w,  sv
                            lower = c(-Inf, -Inf, .01,       0, 0,   0),
                            upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_Gon_17 = nlminb(init_vals[j,], ll_fb_Gon_17, err_tol = err_tol,
                           rt = rti, resp = respi, truth = truthi,
                           # limits:   vu,   vl,   a,      t0, w,  sv
                           lower = c(-Inf, -Inf, .01,       0, 0,   0),
                           upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fs_SWSE_14 = nlminb(init_vals[j,], ll_fs_SWSE_14, err_tol = err_tol,
                            rt = rti, resp = respi, truth = truthi,
                            # limits:   vu,   vl,   a,      t0, w,  sv
                            lower = c(-Inf, -Inf, .01,       0, 0,   0),
                            upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fs_Gon_17 = nlminb(init_vals[j,], ll_fs_Gon_17, err_tol = err_tol,
                           rt = rti, resp = respi, truth = truthi,
                           # limits:   vu,   vl,   a,      t0, w,  sv
                           lower = c(-Inf, -Inf, .01,       0, 0,   0),
                           upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fl_Nav_09 = nlminb(init_vals[j,], ll_fl_Nav_09, err_tol = err_tol,
                           rt = rti, resp = respi, truth = truthi,
                           # limits:   vu,   vl,   a,      t0, w,  sv
                           lower = c(-Inf, -Inf, .01,       0, 0,   0),
                           upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        rtdists = nlminb(init_vals[j,], ll_RTDists,
                         rt = rti, resp = respi, truth = truthi,
                         # limits:   vu,   vl,   a,      t0, w,  sv
                         lower = c(-Inf, -Inf, .01,       0, 0,   0),
                         upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        times = times, unit = unit
      )
      for (k in 1:nalgos) {
        res[["BmTime"]][(i-1)*ni+(k-1)*ninit_vals+j] <- median(
          mbm[mbm[["expr"]] == algo_names[k], 2])
      }
    }
  }
  return(res)
}

## ----fit-run, eval=FALSE------------------------------------------------------
#  data(med_dec, package = "fddm")
#  med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]
#  fit <- rt_fit(med_dec, id_idx = c(2,1), rt_idx = 8, response_idx = 7,
#                truth_idx = 5, response_upper = "blast", err_tol = 1e-6,
#                times = 25, unit = "ns")

## ----fit-run-save, eval=FALSE, include=FALSE----------------------------------
#  save(fit, compress = "xz", compression_level = 9,
#       file = "inst/extdata/bm_fit.Rds")

## ----fit-packages, eval=TRUE--------------------------------------------------
library("reshape2")
library("ggplot2")

## ----fitting-prep, eval=TRUE--------------------------------------------------
fit_prep <- function(fit, eps = 1e-4) {
  nr <- nrow(fit)
  fit[["Obj_diff"]] <- rep(0, nr)

  ids <- unique(fit[["ID"]])
  nids <- length(ids)
  algos <- unique(fit[["Algorithm"]])
  nalgos <- length(algos)

  ninit <- nrow(fit[fit[["ID"]] == ids[1] & fit[["Algorithm"]] == algos[1], ])
  for (i in 1:nids) {
    for (j in 1:ninit) {
      idx <- which(fit[["ID"]] == ids[i])[ninit*(0:(nalgos-1)) + j]
      objs <- fit[idx, "Objective"]
      min_obj <- min(objs)
      abs_min_obj <- abs(min_obj)
      obj_diffs <- objs - min(objs)
      fit[idx, "Obj_diff"] <- ifelse(obj_diffs <= eps*abs_min_obj, 0,
        ifelse(obj_diffs > eps*abs_min_obj & obj_diffs <= 2*abs_min_obj, 1, 3))
    }
  }

  fit[["BmTime"]] <- fit[["BmTime"]]*1e-6 # convert to milliseconds
  fit[["Convergence"]] <- ifelse(fit[["Convergence"]] < 1, 0, 1)

  return(fit)
}

# load data, will be in the variable 'fit'
load(system.file("extdata", "bm_fit.Rds", package = "fddm", mustWork = TRUE))
fit <- fit_prep(fit)

Names <- c("fb_SWSE_17", "fb_Gon_17", "fs_SWSE_14",
           "fs_Gon_17", "fl_Nav_09", "rtdists")
Color <- c("#e000b4", "#e68a00", "#cc99ff",
          "#c2a500", "#996633", "#990000")
Shape <- c(21, 25)
Sizes <- c(0, 3, 3)
Stroke <- c(0, 1, 1)
Fills <- c("#ffffff00", "#ffffff00", "#80808099")

## ----fit-mbm, eval=TRUE, fig.height=6-----------------------------------------
fit_mbm <- melt(fit, id.vars = c("Algorithm", "Convergence", "Obj_diff"),
                measure.vars = "BmTime", value.name = "BmTime")[,-4]

mibm <- min(fit[fit[["Algorithm"]] != "rtdists", "BmTime"])
mabm <- max(fit[fit[["Algorithm"]] != "rtdists", "BmTime"])

ggplot(fit_mbm, aes(x = factor(Algorithm, levels = Names),
                    y = BmTime)) +
  geom_violin(trim = TRUE, alpha = 0.5,
              aes(color = factor(Algorithm, levels = Names),
                  fill = factor(Algorithm, levels = Names))) +
  geom_boxplot(width = 0.2, outlier.shape = NA,
               fill = "white", alpha = 0.4,
               aes(color = factor(Algorithm, levels = Names))) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .5, linetype = "dashed",
               color = Color) +
  scale_color_manual(values = Color, guide = FALSE) +
  scale_shape_manual(values = Shape,
                     name = "Convergence Code",
                     breaks = c(0, 1),
                     labels = c("Success", "Failure")) +
  scale_size_manual(values = Sizes, guide = FALSE) +
  scale_discrete_manual(aesthetics = "stroke", values = Stroke, guide = FALSE) +
  scale_fill_manual(values = Color, guide = FALSE) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = Fills,
                    name = "Objective Difference",
                    breaks = c(1, 2, 3),
                    labels = c("< 2", "NA", "> 2")) +
  geom_point(aes(color = factor(Algorithm, levels = Names),
                 shape = factor(Convergence, levels = c(0, 1)),
                 size = factor(Obj_diff, levels = c(0, 1, 3)),
                 stroke = factor(Obj_diff, levels = c(0, 1, 3)),
                 fill = factor(Obj_diff, levels = c(0, 1, 3)))) +
  coord_cartesian(ylim = c(mibm, mabm)) +
  labs(title = "Median microbenchmark times for data fitting",
       subtitle = paste(
         "Dashed lines represent mean benchmark times",
         "Visible points have an objective difference greater than 1e-4",
         sep = ";\n"),
       x = "Method", y = "Time (milliseconds)") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = Sizes[c(2, 3)])),
         fill = guide_legend(order = 2,
                             override.aes = list(size = Sizes[c(2, 3)],
                                                 shape = c(21, 21),
                                                 fill = Fills[c(2, 3)]))) +
  theme_bw() +
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 10, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20,
                                    margin = margin(5, 10, 5, 5, "pt")),
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))

## ----fit-fev, eval=TRUE, fig.height=6-----------------------------------------
fit_fev <- melt(fit, id.vars = c("Algorithm", "Convergence", "Obj_diff"),
                measure.vars = "FuncEvals", value.name = "FuncEvals")[,-4]

mife <- min(fit[fit[["Algorithm"]] != "rtdists", "FuncEvals"])
mafe <- max(fit[fit[["Algorithm"]] != "rtdists", "FuncEvals"])

ggplot(fit_fev, aes(x = factor(Algorithm, levels = Names),
                    y = FuncEvals)) +
  geom_violin(trim = TRUE, alpha = 0.5,
              aes(color = factor(Algorithm, levels = Names),
                  fill = factor(Algorithm, levels = Names))) +
  geom_boxplot(width = 0.2, outlier.shape = NA,
               fill = "white", alpha = 0.4,
               aes(color = factor(Algorithm, levels = Names))) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .5, linetype = "dashed",
               color = Color) +
  scale_color_manual(values = Color, guide = FALSE) +
  scale_shape_manual(values = Shape,
                     name = "Convergence Code",
                     breaks = c(0, 1),
                     labels = c("Success", "Failure")) +
  scale_size_manual(values = Sizes, guide = FALSE) +
  scale_discrete_manual(aesthetics = "stroke", values = Stroke, guide = FALSE) +
  scale_fill_manual(values = Color, guide = FALSE) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = Fills,
                    name = "Objective Difference",
                    breaks = c(1, 2, 3),
                    labels = c("< 2", "NA", "> 2")) +
  geom_point(aes(color = factor(Algorithm, levels = Names),
                 shape = factor(Convergence, levels = c(0, 1)),
                 size = factor(Obj_diff, levels = c(0, 1, 3)),
                 stroke = factor(Obj_diff, levels = c(0, 1, 3)),
                 fill = factor(Obj_diff, levels = c(0, 1, 3)))) +
  coord_cartesian(ylim = c(mife, mafe)) +
  labs(title = "Function evaluations for data fitting",
       subtitle = paste(
         "Dashed lines represent mean benchmark times",
         "Visible points have an objective difference greater than 1e-4",
         sep = ";\n"),
       x = "Method", y = "Number of function evaluations") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = Sizes[c(2, 3)])),
         fill = guide_legend(order = 2,
                             override.aes = list(size = Sizes[c(2, 3)],
                                                 shape = c(21, 21),
                                                 fill = Fills[c(2, 3)]))) +
  theme_bw() +
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 10, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20,
                                    margin = margin(5, 10, 5, 5, "pt")),
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))

## ----session-info, collapse=TRUE----------------------------------------------
sessionInfo()

