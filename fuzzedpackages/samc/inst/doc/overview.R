## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridis")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

library("raster")
library("samc")
library("viridis")

## ---- fig.show='hold'---------------------------------------------------------
str(samc::ex_res_data)
str(samc::ex_abs_data)
str(samc::ex_occ_data)


plot(raster(samc::ex_res_data, xmn = 1, xmx = ncol(samc::ex_res_data), ymn = 1, ymx = nrow(samc::ex_res_data)),
     main = "Example Resistance Data", xlab = "x", ylab = "y", col = viridis(256))

plot(raster(samc::ex_abs_data, xmn = 1, xmx = ncol(samc::ex_abs_data), ymn = 1, ymx = nrow(samc::ex_abs_data)),
     main = "Example Absorption Data", xlab = "x", ylab = "y", col = viridis(256))

plot(raster(samc::ex_occ_data, xmn = 1, xmx = ncol(samc::ex_occ_data), ymn = 1, ymx = nrow(samc::ex_occ_data)),
     main = "Example Occupancy Data", xlab = "x", ylab = "y", col = viridis(256))

## ----table1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'--------
tabl <- "
| Function | Equation | Description |
|:---------|:---------|:------------|
| `cond_passage()` | $\\tilde{t} = \\tilde{B}_j^{-1}\\tilde{F}\\tilde{B}_j{\\cdot}1$ | Mean first conditional passage time |
| `dispersal()` | $\\tilde{D}_{jt}=({\\sum}_{n=0}^{t-1}\\tilde{Q}^n)\\tilde{q}_j$ | Probability of an individual visiting a location, if starting at any other location, before or at time *t* |
| | $\\psi^T\\tilde{D}_{jt}$ | Probability of an individual visiting a location, before or at time *t*, regardless of initial location |
| | $D=(F-I)diag(F)^{-1}$ | Probability of an individual visiting a location |
| | $\\psi^TD$ | Probability of an individual visiting a location, regardless of initial location |
| `distribution()` | $Q^t$   | Probability of an individual being at a location at time *t* |
| | $\\psi^TQ^t$ | Probability of an individual being at a location at time *t*, regardless of initial location |
| `mortality()` | $\\tilde{B}_t = (\\sum_{n=0}^{t-1} Q^n) \\tilde{R}$ | Probability of an individual experiencing mortality at a location before or at time *t* |
| | $\\psi^T \\tilde{B}_t$ | Probability of an individual experiencing mortality at a location, before or at time *t*, regardless of initial location |
| | $B = F \\tilde{R}$ | Probability of an individual experiencing mortality at a location |
| | $\\psi^T B$ | Probability of an individual experiencing mortality at a location, regardless of initial location |
| `survival()` | $z=(I-Q)^{-1}{\\cdot}1=F{\\cdot}1$ | Expected life expectancy of an individual |
| | ${\\psi}^Tz$ | Overall life expectancy, regardless of initial location |
| `visitation()` | $F = (I-Q)^{-1}$ | Expected number of times an individual visits a location |
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

