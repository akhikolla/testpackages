## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("TORCH_TEST", unset = "0"), "1")
)

## ----setup--------------------------------------------------------------------
#  library(torch)

## -----------------------------------------------------------------------------
#  torch_tensor(c(1,2,3))
#  
#  # conform to row-major indexing used in torch
#  torch_tensor(matrix(1:10, ncol = 5, nrow = 2, byrow = TRUE))
#  torch_tensor(array(runif(12), dim = c(2, 2, 3)))

## -----------------------------------------------------------------------------
#  torch_tensor(1, dtype = torch_long())
#  torch_tensor(1, device = "cpu", dtype = torch_float64())

## -----------------------------------------------------------------------------
#  x <- torch_randn(5, 3)
#  x

## -----------------------------------------------------------------------------
#  x <- torch_ones(2, 4, dtype = torch_int64(), device = "cpu")
#  x

## -----------------------------------------------------------------------------
#  x <- torch_tensor(1)
#  y <- x$to(dtype = torch_int32())
#  x
#  y

## ---- eval = FALSE------------------------------------------------------------
#  x <- torch_tensor(1)
#  y <- x$cuda())

