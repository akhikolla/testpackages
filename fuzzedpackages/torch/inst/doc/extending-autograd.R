## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("TORCH_TEST", unset = "0"), "1")
)

## ----setup--------------------------------------------------------------------
#  library(torch)

## -----------------------------------------------------------------------------
#  linear <- autograd_function(
#    forward = function(ctx, input, weight, bias = NULL) {
#      ctx$save_for_backward(input = input, weight = weight, bias = bias)
#      output <- input$mm(weight$t())
#      if (!is.null(bias))
#        output <- output + bias$unsqueeze(0)$expand_as(output)
#  
#      output
#    },
#    backward = function(ctx, grad_output) {
#  
#      s <- ctx$saved_variables
#  
#      grads <- list(
#        input = NULL,
#        weight = NULL,
#        bias = NULL
#      )
#  
#      if (ctx$needs_input_grad$input)
#        grads$input <- grad_output$mm(s$weight)
#  
#      if (ctx$needs_input_grad$weight)
#        grads$weight <- grad_output$t()$mm(s$input)
#  
#      if (!is.null(s$bias) && ctx$needs_input_grad$bias)
#        grads$bias <- grad_output$sum(dim = 0)
#  
#      grads
#    }
#  )

## -----------------------------------------------------------------------------
#  mul_constant <- autograd_function(
#    forward = function(ctx, tensor, constant) {
#      ctx$save_for_backward(constant = constant)
#      tensor * constant
#    },
#    backward = function(ctx, grad_output) {
#      v <- ctx$saved_variables
#      list(
#        tensor = grad_output * v$constant
#      )
#    }
#  )

## -----------------------------------------------------------------------------
#  x <- torch_tensor(1, requires_grad = TRUE)
#  o <- mul_constant(x, 2)
#  o$backward()
#  x$grad

