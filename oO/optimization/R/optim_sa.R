#--------------------------------------------------#
#### Simulated Annealing Optimization Algorithm ####
#--------------------------------------------------#

# Updated 20.11.2016 #

optim_sa <- function (fun, start, maximization = FALSE, trace = FALSE ,lower, upper, control = list()) {

  #------------------#
  ## Initialisation ##
  #------------------#

  ## User declarations ##

  # start: Vector with starting values
  # fun: Function to be optimized
  # trace: Trace Matrix
  # lower, upper: Boundaries of the x variables
  # control: List with optional assignments


  ## Declaration of variables, vectors and matrices ##


  # Default control assignments:
  con <- list(
    vf      = NULL,
    rf      = 1,
    dyn_rf  = TRUE,
    t0      = 1000,
    nlimit  = 100,
    r       = 0.6,
    k       = 1,
    t_min   = 0.1,
    maxgood = 100,
    stopac  = 30,
    ac_acc  = NA
  )
  corr_names <- names(con) # Saving the correct names for later consistency checks.

  con[names(control)] <- control # Overwriting the default assignments with user declarations.

  # Saving the list objects as stand alone variables.
  vf      <- con$vf
  rf      <- con$rf
  dyn_rf  <- con$dyn_rf
  t0      <- con$t0
  nlimit  <- con$nlimit
  r       <- con$r
  k       <- con$k
  t_min   <- con$t_min
  maxgood <- con$maxgood
  stopac  <- con$stopac
  ac_acc  <- con$ac_acc


  # Declaration of Intern variables:
  fun_length   <- length(start) # Number of parameters of loss function.
  temp         <- t0 # Starting temperature.
  para_0       <- start # Storing vector for the current parameter combination. Initialy with starting values.
  para_opt     <- start # Storing vector for current best parameter combination. Initialy with starting values.
  para_i       <- rep(NA, fun_length) # Vector with length = number of parameters
  loss_i       <- NA # Result of the loss function
  loss_0       <- fun(start) # Result of the loss function at initially the parameter comb.
  delta        <- NA
  goodcounter  <- 0
  ac           <- 0
  n_outer      <- 0 # Counter vector for the outer while loop
  n_inner      <- 0 # Counter vector for the inner loop (necessary since the inner loop has stop criteria and is not always of length 'nlimit').
  n_oob        <- rep(0, fun_length) # Counter vector for parameters out of bounds after changing function.
  ratio_noob   <- rep(0, fun_length) # Ratio between iterations of parameters NOT out of bounds and total number of iterations in the inner loop.
  vf_user      <- FALSE # Bool variable: Is there a user declared variance function?

  # The objective result at start.
  if (is.na(fun(start))) {
    stop ("The objective result at initial variables combination must be valid.")
  } else {
    loss_opt <- fun(start)
  }

  if (trace) {
    trace_array <- matrix(nrow = 0, ncol = 5 + fun_length, dimnames = list(character (0),
                                                                           c("iteration_outer",
                                                                             "function_value",
                                                                             paste("x", c(1 : fun_length), sep="_"),
                                                                             "n_iterations_inner",
                                                                             "temperature",
                                                                             "good_counter")
                                                                           )
                          ) # Saving matrix for the trace.
  } else {
    trace_array <- NULL # Empty variable. Necessary because the output is expected to have the same dimension even if trace is FALSE.
  }


  ## Check for consistency of user declarations ##

  # Check if there are start values and boundaries.
  if (!exists("start") &&  (length(lower) == 0 || length(upper) == 0)) {
    stop ("Starting values or boundaries are not defined.")
  }

  if (length(lower) != length(upper) || length(start) != length(upper) || length(start) != length(lower)) {
    stop ("Starting values or boundaries vector do not have the same length.")
  }

  if (any(start < lower | start > upper)) {
    start[start < lower | start > upper] <- apply(cbind(lower[start < lower | start > upper], upper[start < lower | start > upper]), 1, mean)
    warning ("At least one starting value was out of bounds. These were set to mean of their boundaries.", call. = FALSE)
    }

  if (mode(vf) == "function") {
    var_func <- vf
    vf_user <- TRUE
  } else {
    #var_func <- function() {}
     var_func <- function (para_0, fun_length, rf) {
       ret_var_func <- para_0 + runif(fun_length, 0.000001, rf) *  ((rbinom(fun_length, 1, 0.5) * -2) + 1)
       return (ret_var_func)
       }
  }

  if (r >=1) {
    r <- 0.9
    warning("r must be < 1. It is set to 0.9", call. = FALSE)
  }

  # The user can choose wheather rf is a scalar or a vector. A scalar will be extended (repeated) to the required length.
  # A vector of wrong length conditions a warning message and only the first entry is considered.
  if (length(rf) == 1) {
    rf <- rep(rf,fun_length)
  } else {
      if (!length (rf) == fun_length) {
        rf <- rep(rf[1], fun_length)
        warning ("rf was of wrong length. Only first value was conisidered.", call. = FALSE)
        }
  }

  # Are there wrong variable names in the control list?
  testnames <- names (control)
  if (length(testnames[!testnames %in% corr_names])!=0) {
    warning ("Control contains wrong names.", call. = FALSE)
    }

  # If ac_acc is not initialized by user, it will be stated relatively to the dimension of the initial y.
  if (is.na (ac_acc)) {ac_acc <- fun(start) / 10000}

  #----------------#
  ## Optimization ##
  #----------------#

  # Calling the Cpp source
  result <- .Call('_optimization_main_loop',
              temp = temp,
              t_min = t_min,
              r = r,
              fun_length = fun_length,
              nlimit = nlimit,
              para_0 = para_0,
              para_i = para_i,
              var_func = var_func,
              vf_user = vf_user,
              trace = trace,
              rf = rf,
              lower = lower,
              upper = upper,
              fun = fun,
              loss_0 = loss_0,
              k = k,
              loss_opt = loss_opt,
              para_opt = para_opt,
              dyn_rf = dyn_rf,
              maxgood = maxgood,
              ac_acc,
              stopac = stopac,
              maximization = maximization,
              PACKAGE = 'optimization')

  #----------------------------------------#
  ## Postprocessing and output generation ##
  #----------------------------------------#

  if(trace) {
    df_para <- t(as.data.frame(result[["para"]]))
    df_rf <- t(as.data.frame(result["rf"]))
    rownames(df_para) <- NULL; rownames(df_rf) <- NULL
    trace_array <- cbind(result[["n_outer"]], result[["loss"]], df_para, result[["n_inner"]], result[["temp"]], result[["goodcounter"]], df_rf)
    colnames(trace_array) <- c("n_outer", "loss", paste("x", c(1 : fun_length), sep = "_"), "n_inner", "temperature", "goodcounter", paste("rf", c(1 : fun_length), sep = "_"))
  } else {
    trace_array <- NULL
  }

  output <- list(par            = result[["para_opt"]],
                 function_value = result[["loss_opt"]],
                 trace          = trace_array,
                 fun            = fun,
                 start          = start,
                 lower          = lower,
                 upper          = upper,
                 control        = list(
                   rf      = rf,
                   t0      = t0,
                   nlimit  = nlimit,
                   r       = r,
                   k       = k,
                   t_min   = t_min,
                   maxgood = maxgood,
                   stopac  = stopac,
                   ac_acc  = ac_acc
                   )
                 )

  #class(output) <- append(class(output), "optim_nmsa")
  class(output) <- "optim_nmsa"
  return (output)
}
