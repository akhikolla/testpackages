## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)

## ----setup--------------------------------------------------------------------
library(fastpos)
set.seed(19950521)

## ----message=TRUE, warning=TRUE, paged.print=TRUE-----------------------------
find_critical_pos(rhos = seq(.1, .7, .1), sample_size_max = 1e3,
                  n_studies = 10e3)

## ----fig.height=4.8, fig.width=6.4--------------------------------------------
pop <- create_pop(0.5, 1000000)
pos <- simulate_pos(x_pop = pop[,1],
                    y_pop = pop[,2],
                    n_studies = 10000,
                    sample_size_min = 20,
                    sample_size_max = 1000,
                    replace = T,
                    lower_limit = 0.4,
                    upper_limit = 0.6)
hist(pos, xlim = c(0, 1000), xlab = c("Point of stability"),
     main = "Histogram of points of stability for rho = .5+-.1")
quantile(pos, c(.8, .9, .95), na.rm = T)

## ----parallel1, message=FALSE, warning=FALSE----------------------------------
onecore <- function() {find_critical_pos(0.5)}
multicore <- function() {find_critical_pos(0.5, n_cores = future::availableCores())}
microbenchmark::microbenchmark(onecore(), multicore(), times = 10)

## ----parallel2, message=FALSE, warning=FALSE----------------------------------
onecore <- function() {find_critical_pos(0.5, n_studies = 1e5)}
multicore <- function() {find_critical_pos(0.5, n_studies = 1e5,
                                           n_cores = future::availableCores())}
microbenchmark::microbenchmark(onecore(), multicore(),
                               times = 10)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  speedup <- function(n_max, n_min){
#    (n_max*(n_max+1)/2-n_min*(n_min-1)/2)/(2*n_max-n_min)
#  }
#  speedup(1000, 20)
#  speedup2 <- function(n_max, n_min, cpos){
#    (n_max*(n_max+1)/2-n_min*(n_min-1)/2)/(2*n_max-cpos)
#  }
#  speedup2(1000, 20, 119)
#  # get the 50% quantile to estimate how long on average it takes for fastpos to stop.

