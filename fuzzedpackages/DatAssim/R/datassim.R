#' Data Assimilation 
#' 
#' This function estimates a variable of interest through Data Assimilation 
#' technique by incorporating results from previous assessments.
#'
#' @param X Matrix of predictions, with \code{n}  number of rows as the number of observations, 
#' and \code{t} number of columns as the number of time points from which data were collected.
#' @param Var Matrix of corresponding prediction variances, same dimension as \code{X}.
#' @param Corr Matrix or value of correlations between observations from different time points, by default \code{Corr} = 0.
#'
#' @return 
#' \code{$weights} \tab Estimated Kalman weights according to Eq.[7] in Ehlers et al. (2017).
#' 
#' \code{$PreDA} \tab Predicted values through Data Assimilation according Eq.[5] in Ehlers et al. (2017).
#' 
#' \code{$VarDA} \tab Corresponding estiamted variances according Eq.[6] in Ehlers et al. (2017).
#' 
#' Correlation matrix \code{$Correlation}.
#' @examples
#' Prediction = dataset$Prediction; # Predicted proportions of broadleaf trees 
#' Variance = dataset$Variance; # Corresponding prediction variances
#' 
#' datassim(X = Prediction, Var = Variance); # Corr = 0 by default
#' datassim(Prediction, Variance, 0.5); # Corr = 0.5
#' 
#' Corr = cor(Prediction);
#' datassim(Prediction, Variance, Corr);
#'
#' @references Ehlers, S., Saarela, S., Lindgren, N., Lindberg, E., Nyström, M., Grafström, A.,
#' Persson, H., Olsson, H. & Ståhl, G. (2017). Assessing error correlations in remote sensing-based 
#' predictions of forest attributes for improved data assimilation. Remote Sensing. (In Press)
#' @export
datassim <- function(X, Var, Corr=0) {
  n = nrow(X);
  t = ncol(X);
  X = as.matrix(X);
  C = as.matrix(Corr);
  Var = as.matrix(Var);
  if(ncol(C)==1) {
    Corr = matrix(C[,1], t, t);
    diag(Corr) = 1;
  } else {
    Corr = C;
  }
  return(datassimcpp(X, Var, Corr));
}