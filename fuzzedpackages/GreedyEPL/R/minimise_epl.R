MinimiseEPL <- function(sample_of_partitions, pars = list())
{
  N <- ncol(sample_of_partitions)
  niter <- nrow(sample_of_partitions)
  if (missing(pars) || is.null(pars$weights)) weights = rep(1,niter) else weights = pars$weights
  # cat("\nThe following weigths are being used ", weights)
  if (missing(pars) || is.null(pars$Kup))  Kup = N else Kup = pars$Kup
  if (missing(pars) || is.null(pars$decision_init)) 
  {
    # cat("\nCreating a random starting partition with ", Kup, " groups")
    decision_init = sample(x = 1:Kup, size = N, replace = T)
  }
  else decision_init = as.numeric(pars$decision_init)
  if (missing(pars) || is.null(pars$loss_type)) loss_type = "VI" else loss_type = pars$loss_type
  
  if (sum(decision_init <= 0) > 0) stop("Negative entries in decision_init")

  for (iter in 1:niter) sample_of_partitions[iter,] = CollapseLabels(decision = sample_of_partitions[iter,])
  decision_init = CollapseLabels(decision = decision_init)
  
  if (loss_type == "VI") output <- p__MinimiseAverageVI(sample_of_partitions = sample_of_partitions - 1, weights = weights, decision_init = decision_init - 1)
  else if (loss_type == "B") output <- p__MinimiseAverageB(sample_of_partitions = sample_of_partitions - 1, weights = weights, decision_init = decision_init - 1)
  else if (loss_type == "NVI") output <- p__MinimiseAverageNVI(sample_of_partitions = sample_of_partitions - 1, weights = weights, decision_init = decision_init - 1)
  else if (loss_type == "NID") output <- p__MinimiseAverageNID(sample_of_partitions = sample_of_partitions - 1, weights = weights, decision_init = decision_init - 1)
  else stop("Loss function not recognised")
  
  list(EPL_stored_values = as.numeric(output$EPL_stored_values), EPL = output$EPL, decision = CollapseLabels(as.numeric(output$decision + 1)))
}




