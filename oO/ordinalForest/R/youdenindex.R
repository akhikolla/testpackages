youdenindex <-
function(ytest, ytestpred, categ) {
  
  # Number of true positives:
  ntp <- sum((ytest==categ) & (ytestpred==categ))
  # Number of positives:
  np <- sum(ytest==categ)
  
  # Number of true negatives:
  ntn <- sum((ytest!=categ) & (ytestpred!=categ))
  # Number of negatives:
  nn <- sum(ytest!=categ)
  
  # Sensitivity and specificity:
  sens <- ntp/np; spec <- ntn/nn
  
  # Return value of Youden's J statistic:
  return(sens + spec - 1)
  
}
