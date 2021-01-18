put_labels_to_coefficients <-function(beta, labels){
  
  beta_ <- as.data.frame(beta)
  
  names(beta_)<- labels
  
  beta_
}

is_there_parameter <-function(parameters, par_name)
{
  sum(match(names(parameters),par_name,nomatch=0))
}

