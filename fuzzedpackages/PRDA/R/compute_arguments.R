#################################
####    Compute arguments    ####
#################################


#----    compute_eigen_matrix    ----

# Compute the Eigen Matrix that is used to simulate observation from a bivariate
# normal distribution with a given correlation value (i.e., effect_target).

compute_eigen_matrix <- function(effect_target){

  Sigma <- matrix(c(1,effect_target,effect_target,1), ncol = 2)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  Eign_matrix <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), 2)

  return(Eign_matrix)
}


#----    compute_errors    ----

# Given the sampled p-values, estimated effects, true effect value, significance
# level, and number of iterations, compute the power level, Type M error, and
# Type S error.

compute_errors <- function(p_values, estimates, true_value, sig_level,B){

  sig_p.value <- p_values < sig_level
  sum_sig_p <- sum(sig_p.value)

  power <- sum_sig_p / B
  typeS <- 0

  if(true_value >= 0) {
    typeS <- sum(sig_p.value & (estimates < 0)) / sum_sig_p
  } else {
    typeS <- sum(sig_p.value & (estimates > 0)) / sum_sig_p}

  typeM <- mean(abs(estimates[sig_p.value])) / abs(true_value)

  res <- data.frame(power = power,
                    typeM = typeM,
                    typeS = typeS)

  return(res)
}


#----    compute_df    ----

# Given the effect type, sample size of the first group and second group (when
# needed), standard deviation ratio between the two groups, and test method,
# compute the degrees of freedom of the test.

compute_df <- function(effect_type, sample_n1, sample_n2 = NULL,
                       ratio_sd =1, test_method){
  df <- NULL
  if(effect_type == "cohen_d"){
    if (test_method == "two_sample") {
      df <- sample_n1 + sample_n2 - 2L
    } else if (test_method == "welch"){
      df1 <- sample_n1-1L
      df2 <- sample_n2-1L
      var1 <- ratio_sd^2
      var2 <- 1
      ratio1 <- var1/sample_n1
      ratio2 <- var2/sample_n2
      df <- (ratio1 + ratio2)^2/(ratio1^2/df1 + ratio2^2/df2)
    } else if (test_method %in% c("one_sample","paired")){
      df <- sample_n1-1L
    }
  } else if(effect_type == "correlation"){
    df <- sample_n1 - 2L
  }

  return(df)
}


#----    compute_critical_t    ----

# Given the degrees of freedom, significance level, and alternative hypothesis,
# compute the critical value of the t-statistic.

compute_critical_t <- function(df, sig_level, alternative = "two_sided"){

  critical_t <- NULL

  if(alternative == "two_sided"){
    critical_t <- qt(1-sig_level/2, df)
  } else if(alternative == "greater"){
    critical_t <- qt(1-sig_level, df)
  } else if(alternative == "less"){
    critical_t <- qt(sig_level, df)
  }

  return(critical_t)
}


#----    compute_critical_effect    ----

# Given the effect type, sample size of the first group and second group (when
# needed), test method, significance level, alternative hypothesis, standard
# deviation ratio between the two groups, and value of the null hypothesis,
# compute the degrees of freedom of the test and critical effect value (i.e.,
# the minimum absolute effect size value that would result significant).

compute_critical_effect <- function(effect_type, sample_n1, sample_n2 = NULL,
                                    test_method, sig_level, alternative,
                                    ratio_sd = 1, mu = 0, ...){

  df <- compute_df(effect_type = effect_type,
                  sample_n1 = sample_n1,
                  sample_n2 = sample_n2,
                  test_method = test_method,
                  ratio_sd = ratio_sd)

  critical_t <- compute_critical_t(df, sig_level, alternative)

  critical_effect <- NULL

  if(effect_type == "cohen_d"){
    if (test_method == "one_sample") {
      critical_effect <- critical_t /sqrt(sample_n1)
    } else if (test_method == "paired"){
      critical_effect <- critical_t /sqrt(sample_n1) + mu
    } else if (test_method == "two_sample") {
      critical_effect <- critical_t * sqrt(
        (sample_n1+sample_n2)/(sample_n1*sample_n2)) + mu
    } else if (test_method == "welch") {
      var1 <- ratio_sd^2
      var2 <- 1
      critical_effect <- critical_t * sqrt(2/(sample_n1 * sample_n2) * (
        var1*sample_n2 + var2*sample_n1)/(var1 + var2)) + mu
    }
  } else if(effect_type == "correlation"){
    critical_effect <- critical_t /sqrt(sample_n1-2+critical_t^2)
  }

  res <- list(df = df, critical_effect = critical_effect)

  return(res)
}


#----    compute_correct_diff    ----

# Given the effect size value (i.e, Cohen's d), test method, and standard
# deviation ratio between the two groups, compute the correct mean difference
# between groups.

compute_correct_diff <- function(effect_target, test_method, ratio_sd){
  if(test_method == "paired"){
    correct_diff <- effect_target * sqrt(2)
  } else if(test_method == "welch"){
    correct_diff <- effect_target * sqrt((ratio_sd^2 + 1)/2)
  } else {
    correct_diff <- effect_target
  }

  return(correct_diff)
}



#----



