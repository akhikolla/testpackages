## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("TORCH_TEST", unset = "0"), "1")
)

## ----setup--------------------------------------------------------------------
#  library(torch)

## -----------------------------------------------------------------------------
#  x <- torch_randn(10, 10)
#  torch_save(x, "tensor.pt")
#  x_ <- torch_load("tensor.pt")
#  
#  torch_allclose(x, x_)

## -----------------------------------------------------------------------------
#  module <- nn_module(
#    "my_module",
#    initialize = function() {
#      self$fc1 <- nn_linear(10, 10)
#      self$fc2 <- nn_linear(10, 1)
#    },
#    forward = function(x) {
#      x %>%
#        self$fc1() %>%
#        self$fc2()
#    }
#  )
#  
#  model <- module()
#  torch_save(model, "model.pt")
#  model_ <- torch_load("model.pt")
#  
#  # input tensor
#  x <- torch_randn(50, 10)
#  torch_allclose(model(x), model_(x))

## ----eval = FALSE-------------------------------------------------------------
#  state_dict <- load_state_dict(fpath)
#  model <- Model()
#  model$load_state_dict(state_dict)

