## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("TORCH_TEST", unset = "0"), "1")
)

## ----setup--------------------------------------------------------------------
#  library(torch)

## ---- eval = FALSE------------------------------------------------------------
#  # ...
#  ds <- mnist_dataset(
#    dir,
#    download = TRUE,
#    transform = function(x) {
#      x <- x$to(dtype = torch_float())/256
#      x[newaxis,..]
#    }
#  )
#  
#  dl <- dataloader(ds, batch_size = 32, shuffle = TRUE)
#  
#  for (b in enumerate(dl)) {
#    # ...
#  

## -----------------------------------------------------------------------------
#  library(palmerpenguins)
#  library(magrittr)
#  
#  penguins

## -----------------------------------------------------------------------------
#  penguins_dataset <- dataset(
#  
#    name = "penguins_dataset",
#  
#    initialize = function() {
#      self$data <- self$prepare_penguin_data()
#    },
#  
#    .getitem = function(index) {
#  
#      x <- self$data[index, 2:-1]
#      y <- self$data[index, 1]$to(torch_long())
#  
#      list(x, y)
#    },
#  
#    .length = function() {
#      self$data$size()[[1]]
#    },
#  
#    prepare_penguin_data = function() {
#  
#      input <- na.omit(penguins)
#      # conveniently, the categorical data are already factors
#      input$species <- as.numeric(input$species)
#      input$island <- as.numeric(input$island)
#      input$sex <- as.numeric(input$sex)
#  
#      input <- as.matrix(input)
#      torch_tensor(input)
#    }
#  )
#  

## -----------------------------------------------------------------------------
#  tuxes <- penguins_dataset()
#  tuxes$.length()
#  tuxes$.getitem(1)

## -----------------------------------------------------------------------------
#  dl <-tuxes %>% dataloader(batch_size = 8)

## -----------------------------------------------------------------------------
#  dl$.length()

## -----------------------------------------------------------------------------
#  iter <- dl$.iter()
#  b <- iter$.next()
#  b

## -----------------------------------------------------------------------------
#  net <- nn_module(
#    "PenguinNet",
#    initialize = function() {
#      self$fc1 <- nn_linear(7, 32)
#      self$fc2 <- nn_linear(32, 3)
#    },
#    forward = function(x) {
#      x %>%
#        self$fc1() %>%
#        nnf_relu() %>%
#        self$fc2() %>%
#        nnf_log_softmax(dim = 1)
#    }
#  )
#  
#  model <- net()
#  

## -----------------------------------------------------------------------------
#  optimizer <- optim_sgd(model$parameters, lr = 0.01)

## -----------------------------------------------------------------------------
#  for (epoch in 1:10) {
#  
#    l <- c()
#  
#    for (b in enumerate(dl)) {
#      optimizer$zero_grad()
#      output <- model(b[[1]])
#      loss <- nnf_nll_loss(output, b[[2]])
#      loss$backward()
#      optimizer$step()
#      l <- c(l, loss$item())
#    }
#  
#    cat(sprintf("Loss at epoch %d: %3f\n", epoch, mean(l)))
#  }
#  

