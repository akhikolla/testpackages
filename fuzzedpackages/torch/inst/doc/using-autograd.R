## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("TORCH_TEST", unset = "0"), "1")
)

## ----setup--------------------------------------------------------------------
#  library(torch)

## -----------------------------------------------------------------------------
#  x <- torch_ones(2,2, requires_grad = TRUE)

## -----------------------------------------------------------------------------
#  y <- x$mean()

## -----------------------------------------------------------------------------
#  y$grad_fn

## -----------------------------------------------------------------------------
#  y$backward()

## -----------------------------------------------------------------------------
#  x$grad

## -----------------------------------------------------------------------------
#  x1 <- torch_ones(2,2, requires_grad = TRUE)
#  x2 <- torch_tensor(1.1, requires_grad = TRUE)
#  y <- x1 * (x2 + 2)
#  y$retain_grad()
#  z <- y$pow(2) * 3
#  z$retain_grad()
#  out <- z$mean()

## -----------------------------------------------------------------------------
#  # how to compute the gradient for mean, the last operation executed
#  out$grad_fn
#  # how to compute the gradient for the multiplication by 3 in z = y$pow(2) * 3
#  out$grad_fn$next_functions
#  # how to compute the gradient for pow in z = y.pow(2) * 3
#  out$grad_fn$next_functions[[1]]$next_functions
#  # how to compute the gradient for the multiplication in y = x * (x + 2)
#  out$grad_fn$next_functions[[1]]$next_functions[[1]]$next_functions
#  # how to compute the gradient for the two branches of y = x * (x + 2),
#  # where the left branch is a leaf node (AccumulateGrad for x1)
#  out$grad_fn$next_functions[[1]]$next_functions[[1]]$next_functions[[1]]$next_functions
#  # here we arrive at the other leaf node (AccumulateGrad for x2)
#  out$grad_fn$next_functions[[1]]$next_functions[[1]]$next_functions[[1]]$next_functions[[2]]$next_functions

## -----------------------------------------------------------------------------
#  out$backward()
#  z$grad
#  y$grad
#  x2$grad
#  x1$grad

## -----------------------------------------------------------------------------
#  ### generate training data -----------------------------------------------------
#  # input dimensionality (number of input features)
#  d_in <- 3
#  # output dimensionality (number of predicted features)
#  d_out <- 1
#  # number of observations in training set
#  n <- 100
#  # create random data
#  x <- torch_randn(n, d_in)
#  y <- x[,1]*0.2 - x[..,2]*1.3 - x[..,3]*0.5 + torch_randn(n)
#  y <- y$unsqueeze(dim = 1)
#  ### initialize weights ---------------------------------------------------------
#  # dimensionality of hidden layer
#  d_hidden <- 32
#  # weights connecting input to hidden layer
#  w1 <- torch_randn(d_in, d_hidden, requires_grad = TRUE)
#  # weights connecting hidden to output layer
#  w2 <- torch_randn(d_hidden, d_out, requires_grad = TRUE)
#  # hidden layer bias
#  b1 <- torch_zeros(1, d_hidden, requires_grad = TRUE)
#  # output layer bias
#  b2 <- torch_zeros(1, d_out,requires_grad = TRUE)
#  ### network parameters ---------------------------------------------------------
#  learning_rate <- 1e-4
#  ### training loop --------------------------------------------------------------
#  for (t in 1:200) {
#  
#      ### -------- Forward pass --------
#      y_pred <- x$mm(w1)$add(b1)$clamp(min = 0)$mm(w2)$add(b2)
#      ### -------- compute loss --------
#      loss <- (y_pred - y)$pow(2)$mean()
#      if (t %% 10 == 0) cat(t, as_array(loss), "\n")
#      ### -------- Backpropagation --------
#      # compute the gradient of loss with respect to all tensors with requires_grad = True.
#      loss$backward()
#  
#      ### -------- Update weights --------
#  
#      # Wrap in torch.no_grad() because this is a part we DON'T want to record for automatic gradient computation
#      with_no_grad({
#  
#        w1$sub_(learning_rate * w1$grad)
#        w2$sub_(learning_rate * w2$grad)
#        b1$sub_(learning_rate * b1$grad)
#        b2$sub_(learning_rate * b2$grad)
#  
#        # Zero the gradients after every pass, because they'd accumulate otherwise
#        w1$grad$zero_()
#        w2$grad$zero_()
#        b1$grad$zero_()
#        b2$grad$zero_()
#  
#      })
#  
#  }
#  

