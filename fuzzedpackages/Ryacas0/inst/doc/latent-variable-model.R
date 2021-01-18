## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas0)
library(Matrix)

## -----------------------------------------------------------------------------
L1chr <- diag(4)
L1chr[2:4, 1] <- "-a"
L1 <- as.Sym(L1chr)
L1

## -----------------------------------------------------------------------------
L1inv <- Simplify(Inverse(L1))
K1 <- Simplify(Transpose(L1) * L1)
V1 <- Simplify(L1inv * Transpose(L1inv))

## ---- results="asis"----------------------------------------------------------
cat(
  "\\begin{align} 
    K_1 &= ", TeXForm(K1), " \\\\ 
   V_1 &= ", TeXForm(V1), " 
  \\end{align}", sep = "")

## -----------------------------------------------------------------------------
L2chr <- diag(4)
L2chr[2:4, 1] <- c("-a1", "-a2", "-a3")
L2 <- as.Sym(L2chr)
L2
Vechr <- diag(4)
Vechr[cbind(1:4, 1:4)] <- c("w1", "w2", "w2", "w2")
Ve <- as.Sym(Vechr)
Ve

## -----------------------------------------------------------------------------
L2inv <- Simplify(Inverse(L2))
K2 <- Simplify(Transpose(L2) * Inverse(Ve) * L2)
V2 <- Simplify(L2inv * Ve * Transpose(L2inv))

## ---- results="asis"----------------------------------------------------------
cat(
  "\\begin{align} 
    K_2 &= ", TeXForm(K2), " \\\\ 
   V_2 &= ", TeXForm(V2), " 
  \\end{align}", sep = "")

