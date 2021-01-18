## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas0)

## -----------------------------------------------------------------------------
eq <- "x^2 + 4 + 2*x + 2*x"
yacas(eq)
yacas(paste0("Simplify(", eq, ")"))
yacas(paste0("Factor(", eq, ")"))

## -----------------------------------------------------------------------------
Simplify(yacas(eq))
Factor(yacas(eq))

## -----------------------------------------------------------------------------
yeq <- yacas(eq)
Simplify(yeq)
Factor(yeq)

## -----------------------------------------------------------------------------
res <- Factor(yeq)
as.character(res)
TeXForm(res)

## -----------------------------------------------------------------------------
xs <- Sym("x")
eqs <- xs^2 + 4 + 2*xs + 2*xs
Simplify(eqs)
Factor(eqs)
res <- Factor(eqs)
as.character(res)
TeXForm(res)

## -----------------------------------------------------------------------------
yacas("p := x^2 + 4 + 2*x + 2*x")
ps <- Sym("p")
ps
TeXForm(Simplify(ps))

