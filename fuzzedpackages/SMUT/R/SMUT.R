SMUT=function(G,mediator,outcome,outcome_type="continuous",method="score",
              approxi=TRUE, debug=FALSE){

  # p_value_theta is the p value for testing theta in the outcome model
  # outcome model: outcome ~ intercept + G*gamma + mediator*theta + error
  p_value_theta=Testing_coefficient_of_mediator(G,mediator,outcome,outcome_type,method,approxi,debug)

  # p_value_beta is the p value for testing beta in the mediator model
  # mediator model: mediator ~ intercept + G*beta + error
  obj=SKAT_Null_Model(mediator~1,out_type="C")
  p_value_beta=SKAT( as.matrix(G), obj,impute.method="bestguess")$p.value

  # p_value_IUT is the p value for testing theta*beta based on intersection-union test
  p_value_IUT=max(c(p_value_theta,p_value_beta))
  return(list("p_value_IUT"=p_value_IUT,"p_value_theta"=p_value_theta,"p_value_beta"=p_value_beta))
}


