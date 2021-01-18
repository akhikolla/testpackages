##################################
####    evaluate arguments    ####
##################################


#----    eval_arguments_retrospective    ----

# Given the arguments of the function retrospective, evaluate that they are
# correctly specified with no conflicts.

eval_arguments_retrospective <- function(effect_size, sample_n1, sample_n2,
                                         effect_type, test_method, sig_level,
                                         ratio_sd, B, tl, tu, B_effect,
                                         display_message, ...){
  # Check inputs arguments
  if(!is.function(effect_size) && !is_single_numeric(effect_size))
    stop("Argument 'effect_size' has to be a single numeric value or a function")

  if(!is_single_numeric(sample_n1) || sample_n1 <= 1 )
    stop("Argument 'sample_n1' has to be a single integer value grater than 1")

  if(!is.null(sample_n2) && (!is_single_numeric(sample_n2) || sample_n2 <= 1))
    stop("If specified, argument 'sample_n2' has to be a single integer value grater than 1")

  if(!is_single_numeric(sig_level) || sig_level >= 1 || sig_level <= 0)
    stop("Argument 'sig_level' has to be a single value between 0 and 1")

  if(!is_single_numeric(ratio_sd) || ratio_sd <= 0)
    stop("Argument 'ratio_sd' has to be a single finite number grater than 0")

  if(!is_single_numeric(B) || B <= 1)
    stop("Argument 'B' has to be a single integer value grater than 1")

  if(!is_single_numeric(tl, infinite = TRUE ))
    stop("Argument 'tl' has to be a single numeric value")

  if(!is_single_numeric(tu, infinite = TRUE))
    stop("Argument 'tu' has to be a single numeric value")

  if(!is_single_numeric(B_effect) || B_effect <= 1)
    stop("Argument 'B_effect' has to be a single integer value grater than 1")

  if(!is.logical(display_message))
    stop("Argument 'display_message' has to be logical")

  # Check coherence effect_type, test_method and others
  if(effect_type == "correlation" && test_method != "pearson")
    stop("If  'effect_type = correlation', argument 'test_method' has to be 'pearson'")

  if(effect_type == "cohen_d" && test_method == "pearson")
    stop("No appropriate 'test_method' for 'effect_type = cohen_d'")

  if(test_method == "one_sample" && !is.null(sample_n2))
    stop("If 'test_method = one_sample', argument 'sample_n2' must be set to NULL")

  if(test_method == 'paired' && ((is.null(sample_n2) || sample_n1!=sample_n2)))
    stop("If 'test_method = paired', arguments 'sample_n1' and 'sample_n2' must be equal")

  if(test_method %in% c("two_sample", "welch") && is.null(sample_n2))
    stop("Argument 'sample_n2' is required for the specified 'test_method'")

  if(!isTRUE(all.equal(ratio_sd, 1)) && test_method != "welch")
    stop("Argument 'ratio_sd' is required only for 'test_method = welch'")

  if(test_method == "welch" && isTRUE(all.equal(ratio_sd, 1)))
    stop("Argument 'ratio_sd' can not be 1 for 'test_method = welch'\n  Consider 'test_method = two_sample' instead")
}


#----    eval_arguments_prospective    ----

# Given the arguments of the function prospective, evaluate that they are
# correctly specified with no conflicts.

eval_arguments_prospective <- function(effect_size, power, ratio_n,
                                       effect_type, test_method,
                                       sig_level, ratio_sd, B, tl, tu,
                                       B_effect, sample_range, tol,
                                       display_message, ...){
  # Check inputs arguments
  if(!is.function(effect_size) && !is_single_numeric(effect_size))
    stop("Argument 'effect_size' has to be a single numeric value or a function")

  if(!is_single_numeric(power) || power >= 1 || power <= 0)
    stop("Argument 'power' has to be a single value between 0 and 1")

  if(!is.null(ratio_n) && (!is_single_numeric(ratio_n) || ratio_n < 0))
    stop("If specified, argument 'ratio_n' has to be a single positive value")

  if(!is_single_numeric(sig_level) || sig_level >= 1 || sig_level <= 0)
    stop("Argument 'sig_level' has to be a single value between 0 and 1")

  if(!is_single_numeric(ratio_sd) || ratio_sd <= 0)
    stop("Argument 'ratio_sd' has to be a single finite number grater than 0")

  if(!is_single_numeric(B) || B <= 1)
    stop("Argument 'B' has to be a single integer value grater than 1")

  if(!is_single_numeric(tl, infinite = TRUE ))
    stop("Argument 'tl' has to be a single numeric value")

  if(!is_single_numeric(tu, infinite = TRUE))
    stop("Argument 'tu' has to be a single numeric value")

  if(!is_single_numeric(B_effect) || B_effect <= 1)
    stop("Argument 'B_effect' has to be a single integer value grater than 1")

  if(length(sample_range)!=2L || sum(!is.finite(sample_range)))
    stop("Argument 'sample_range' has to be a length-2 numeric vector")
  if(sample_range[1] <= 1 || sample_range[1] > sample_range[2])
    stop("Argument 'sample_range' minimum has to be grater than 1 and less than sample range maximum")

  if(!is_single_numeric(tol) || tol <= 0 || tol >= 1)
    stop("Argument 'tol' has to be a single value between 0 and 1")

  if(!is.logical(display_message))
    stop("Argument 'display_message' has to be logical")

  # Check coherence effect_type, test_method and others
  if(effect_type == "correlation" && test_method != "pearson")
    stop("If  'effect_type = correlation', argument 'test_method' has to be 'pearson'")

  if(effect_type == "cohen_d" && test_method == "pearson")
    stop("No appropriate 'test_method' for 'effect_type = cohen_d'")

  if(test_method == "paired" && (is.null(ratio_n) || ratio_n != 1))
    stop("If 'test_method = paired', argument 'ratio_n' has to be 1")

  if(test_method == "one_sample" && !is.null(ratio_n))
    stop("If 'test_method = one_sample', argument 'ratio_n' must be set to NULL")

  if(test_method %in% c("two_sample", "welch") && is.null(ratio_n))
    stop("Argument 'ratio_n' is required for the specified 'test_method'")

  if(!isTRUE(all.equal(ratio_sd, 1)) && test_method != "welch")
    stop("Argument 'ratio_sd' is required only for 'test_method = welch'")

  if(test_method == "welch" && isTRUE(all.equal(ratio_sd, 1)))
    stop("Argument 'ratio_sd' can not be 1 for 'test_method = welch'\n  Consider 'test_method = two_sample' instead")

}


#----    eval_effect_size    ----

# Given the effect type, effect size values, lower and upper truncation limits,
# and number of sampled effects, evaluate that arguments are correctly specified
# with no conflicts and sample effects according to definition.

eval_effect_size <- function(effect_type, effect_size,
                             tl = -Inf, tu = Inf, B_effect = 250){
  correlation <- effect_type == "correlation"
  sample_fun <- is.function(effect_size)

  if(!sample_fun){
    res <- list(effect_function = "single_value",
               effect_summary = summary(effect_size),
               effect_samples = effect_size)

    if(correlation && (effect_size < -1 || effect_size > 1))
      stop("If 'effect_type = correlation', argument 'effect_size' must be between -1 and 1")
  } else {

    if(correlation && (tl < -1 || tu > 1)){
      tl <- max(-1, tl)
      tu <- min(tu, 1)
      message(paste("If 'effect_type = correlation', effect_size distribution is truncated between",
                    tl,"and", tu))
    }

    res <- sample_effect(FUN = effect_size, B_effect = B_effect,
                         tl = tl, tu = tu)
    res <- c(res,
            tl = tl,
            tu = tu)
  }

  return(res)
}


#----    eval_samples    ----

# Given the sample size ratio between the first and second group and the current
# sample size value, compute the sample size of the two groups.

eval_samples <- function(ratio_n, current_n){


  if(is.null(ratio_n)){
    sample_n1 <- current_n
    sample_n2 <- NULL
  } else {
    sample_n1 <- round(current_n * ratio_n,0)
    sample_n2 <- current_n
  }

  return(list(sample_n1 = sample_n1, sample_n2 = sample_n2))
}


#----    eval_test_method    ----

# Given the effect type, effect size value, test method, sample in the first and
# second group (when needed), alternative hypothesis, significance level,
# standard deviation ratio between the two groups, evaluate that arguments are
# correctly specified with no conflicts. One loop in which observation are
# sampled and appropriate test is run.

eval_test_method <- function(effect_type, effect_target, test_method,
                             sample_n1, sample_n2 = NULL,
                             alternative, sig_level, ratio_sd = 1, ...){

  # Define conf.level according to sig_level
  conf.level <- 1 - sig_level

  # Set correct alternative
  if(alternative == "two_sided")
    alternative <- "two.sided"

  # Cohen d
  if(effect_type == "cohen_d"){

    correct_diff <- compute_correct_diff(effect_target, test_method, ratio_sd)
    groups <- sample_groups(sample_n1, correct_diff, sample_n2, ratio_sd)

    paired <- FALSE
    var.equal <- FALSE
    if(test_method == 'paired'){
      paired <- TRUE
    } else if(test_method == 'two_sample'){
      var.equal <- TRUE
    }

    t.test(groups$x, groups$y, paired = paired, var.equal = var.equal,
           alternative = alternative, conf.level = conf.level)


  } else if (effect_type == "correlation"){

    groups <- sample_obs_cor(sample_n1, effect_target)
    cor.test(groups$x, groups$y, method = test_method,
             alternative = alternative, conf.level = conf.level)
    }
}


#----    eval_rgn_function    ----

# Given the function to sample effect size values, evaluate the the function
# returns the required number of numeric values.

eval_rgn_function <- function(FUN, n = 10){

  args <- list(x = n)

  if(names(formals(FUN))!="x")
    names(args) <- names(formals(FUN))

  out <- do.call(FUN, args)

  n_out <- length(out) == as.integer(n)
  res <- is.numeric(out) && n_out

  return(res)
}

#----    eval_effect_type    ----

# Given the test_method return the corresponding effect_type. In the case of
# "two_sample", "welch", "paired", or "one_sample" test return "choen_d. In the
# case of "pearson" test, return "correlation".

eval_effect_type <- function(test_method = c("pearson", "two_sample", "welch",
                                             "paired", "one_sample")){
  mean_diff <- c("two_sample", "welch", "paired", "one_sample")
  correlation <- c("pearson")

  if(test_method %in% mean_diff){
    effect_type <- "cohen_d"
  } else if (test_method %in% correlation){
    effect_type <- "correlation"
  }

  return(effect_type)
}


#----

