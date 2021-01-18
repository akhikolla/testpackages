## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(knitr)
options(scipen=1, digits=3)

## -----------------------------------------------------------------------------
library(leontief)

X <- transaction_matrix
w <- wage_demand_matrix[, "wage"]
c <- wage_demand_matrix[, "household_consumption"]
d <- wage_demand_matrix[, "final_total_demand"]
e <- employment_matrix[, "employees"]

## -----------------------------------------------------------------------------
A <- input_requirement(X, d)
A_aug <- augmented_input_requirement(X,w,c,d)
rownames(A_aug) <- c(rownames(X), "wage_over_demand")
colnames(A_aug) <- c(rownames(X), "consumption_over_demand")
kable(A_aug)

## -----------------------------------------------------------------------------
B <- output_allocation(X, d)
rownames(B) <- rownames(X)
colnames(B) <- rownames(X)
kable(B)

## -----------------------------------------------------------------------------
L <- leontief_inverse(A); 
rownames(L) <- rownames(X)
colnames(L) <- rownames(X)
kable(L)

## -----------------------------------------------------------------------------
eq <- equilibrium_output(L, d)
rownames(eq) <- rownames(X)
colnames(eq) <- "output"
kable(eq)

## -----------------------------------------------------------------------------
out <- output_multiplier(L)

## -----------------------------------------------------------------------------
inc <- income_multiplier(L, w/d)

## -----------------------------------------------------------------------------
emp <- employment_multiplier(L, e/d)

## -----------------------------------------------------------------------------
sm <- round(cbind(out,inc,emp),4)
rownames(sm) <- rownames(X)
colnames(sm) <- c("output_multiplier", "income_multiplier", "employment_multiplier")
kable(sm)

## -----------------------------------------------------------------------------
bl <- backward_linkage(A)
fl <- forward_linkage(A)
bfl <- cbind(bl,fl)
rownames(bfl) <- rownames(X)
colnames(bfl) <- c("backward_linkage", "forward_linkage")
kable(bfl)

## -----------------------------------------------------------------------------
bl <- power_dispersion(L)
bl_cv <- power_dispersion_cv(L)
bl_t <- cbind(bl,bl_cv)
rownames(bl_t) <- rownames(X)
colnames(bl_t) <- c("power_dispersion", "power_dispersion_cv")
kable(bl_t)

## -----------------------------------------------------------------------------
sl <- sensitivity_dispersion(L)
sl_cv <- sensitivity_dispersion_cv(L)
sl_t <- cbind(sl,sl_cv)
rownames(sl_t) <- rownames(X)
colnames(sl_t) <- c("power_dispersion", "power_dispersion_cv")
kable(sl_t)

## -----------------------------------------------------------------------------
mp <- multiplier_product_matrix(L)
rownames(mp) <- rownames(X)
colnames(mp) <- rownames(X)
kable(mp)

## -----------------------------------------------------------------------------
bli <- backward_linkage(A_aug)
fli <- forward_linkage(A_aug)
bfli <- cbind(bli,fli)
rownames(bfli) <- c(rownames(X), "wage")
# wie = with induced effect
colnames(bfli) <- c("backward_linkage_wie", "forward_linkage_wie")
kable(bfli)

