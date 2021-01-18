## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install-cran, eval = FALSE-----------------------------------------------
#  install.packages("bootUR")

## ----install-github, eval = FALSE---------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("smeekes/bootUR")

## ----install-vign, eval = FALSE-----------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("smeekes/bootUR", build_vignettes = TRUE, dependencies = TRUE)

## ----load---------------------------------------------------------------------
library(bootUR)

## ----missing------------------------------------------------------------------
data("MacroTS")
check_missing_insample_values(MacroTS)

## ----start_end----------------------------------------------------------------
sample_check <- find_nonmissing_subsample(MacroTS)
# Provides the number of the first and last non-missing observation for each series:
sample_check$range 
# Gives TRUE if the time series all start and end at the same observation:
sample_check$all_equal

## ----plot_na, fig.height = 4, fig.width = 7-----------------------------------
plot_missing_values(MacroTS, show_names = TRUE, axis_text_size = 5, legend_size = 6)

## ----adf----------------------------------------------------------------------
set.seed(155776)
GDP_NL <- MacroTS[, 4]
adf_out <- boot_df(GDP_NL, B = 399, boot = "SB", dc = 2, detr = c("OLS", "QD"), verbose = TRUE)

## ----union--------------------------------------------------------------------
union_out <- boot_union(GDP_NL, B = 399, boot = "SWB", verbose = TRUE)

## ----panel--------------------------------------------------------------------
panel_out <- paneltest(MacroTS, boot = "DWB", B = 399, verbose = TRUE)

## ----iADF---------------------------------------------------------------------
iADF_out <- iADFtest(MacroTS[, 1:5], boot = "MBB", B = 399, verbose = TRUE, union = FALSE, 
                     dc = 2, detr = "OLS")

## ----BSQT---------------------------------------------------------------------
N <- ncol(MacroTS)
# Test each unit sequentially
BSQT_out1 <- BSQTtest(MacroTS, q = 0:N, boot = "AWB", B = 399, verbose = TRUE)
# Split in four equally sized groups (motivated by the 4 series per country)
BSQT_out2 <- BSQTtest(MacroTS, q = 0:4 / 4, boot = "AWB", B = 399, verbose = TRUE)

## ----bFDR---------------------------------------------------------------------
N <- ncol(MacroTS)
bFDR_out <- bFDRtest(MacroTS[, 1:10], level = 0.1, boot = "BWB", B = 399, verbose = TRUE)

## ----orders-------------------------------------------------------------------
out_orders <- order_integration(MacroTS[, 11:15], test = "bFDRtest", B = 399)
# Orders
out_orders$order_int
# Differenced data
stationary_data <- out_orders$diff_data

## ----plot_orders, fig.width = 7, fig.height = 4-------------------------------
plot_order_integration(out_orders$order_int)

