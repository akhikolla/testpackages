## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>")
library(baggr)
library(ggplot2)
baggr_schools <- baggr(schools, model = "rubin", pooling = "partial")
# baggr_plot(baggr_schools)
# baggr_compare(schools)
# my_baggr_plot <- baggr_compare(schools)

## -----------------------------------------------------------------------------
prepare_ma(microcredit_simplified, outcome = "consumerdurables")

## -----------------------------------------------------------------------------
schools

## ----eval=FALSE---------------------------------------------------------------
#  baggr_schools <- baggr(schools, model = "rubin", pooling = "partial")

## -----------------------------------------------------------------------------
print(baggr_schools)

## ----eval=FALSE---------------------------------------------------------------
#  baggr(schools, "rubin",
#        prior_hypermean = normal(-5, 10),
#        prior_hypersd   = uniform(0, 5))

## ----eval=FALSE---------------------------------------------------------------
#  custom_priors <- list( hypermean = cauchy(0,25), hypersd = normal(0,30))
#  baggr(schools, "rubin", pooling = "partial", prior = custom_priors)

## ----eval=FALSE---------------------------------------------------------------
#  baggr_schools <- baggr(schools, model = "rubin", pooling = "partial",
#                         iter = 10000, chains = 8)

## -----------------------------------------------------------------------------
pooling(baggr_schools)

## ----fig.width=4--------------------------------------------------------------
plot(baggr_schools, order = FALSE)

## ----fig.width = 4------------------------------------------------------------
effect_plot(baggr_schools)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  my_baggr_comparison <- baggr_compare(schools)

## ---- echo=FALSE, eval = TRUE, include = FALSE--------------------------------
my_baggr_comparison <- baggr_compare(schools)

## ----fig.width=5, fig.height=4, echo = TRUE-----------------------------------
plot(my_baggr_comparison)[[1]] + 
  ggtitle("8 schools: model comparison")

## ---- eval=F, echo=T----------------------------------------------------------
#  baggr_schools_v2 <- baggr(schools, prior_hypermean = normal(10, 2.5))

## ---- eval=T, include=F-------------------------------------------------------
baggr_schools_v2 <- baggr(schools, prior_hypermean = normal(10, 2.5))

## ---- fig.width=5, fig.height=5-----------------------------------------------
effect_plot("Default model" = baggr_schools, "normal(10, 2.5)" = baggr_schools_v2) +
  coord_cartesian(xlim = c(-10, 30)) + theme(legend.position = "top")
baggr_compare("Default model" = baggr_schools, "normal(10, 2.5)" = baggr_schools_v2)

## ---- fig.width=5-------------------------------------------------------------
forest_plot(baggr_schools)

## ---- fig.width=5-------------------------------------------------------------
forest_plot(baggr_schools, show = "both")

## ----loocv, echo = T, results = 'hide', warning = F, message = F--------------
loocv_res <- loocv(schools, return_models = FALSE, "rubin", pooling = "partial")

## -----------------------------------------------------------------------------
loocv_res

## -----------------------------------------------------------------------------
names(attributes(loocv_res))
attr(loocv_res, "df")

## ---- echo = FALSE, include = FALSE, results = 'hide'-------------------------
fit1 <- baggr(data = schools[1:7,], test_data = schools[8,], 
              model = "rubin", pooling = "partial")
fit2 <- baggr(data = schools[1:7,], test_data = schools[8,], 
              model = "rubin", pooling = "full")

## -----------------------------------------------------------------------------
fit1$mean_lpd
fit2$mean_lpd

## ----results = 'hide'---------------------------------------------------------
cv_1 <- loocv(data = schools, 
              model = "rubin", 
              pooling = "partial")
cv_2 <- loocv(data = schools, 
              model = "rubin", 
              pooling = "full")

## -----------------------------------------------------------------------------
loo_compare(cv_1, cv_2)

