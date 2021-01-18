#' Ordinal Forests: Prediction and Variable Ranking with Ordinal Target Variables
#'
#' The ordinal forest (OF) method allows ordinal regression with high-dimensional
#' and low-dimensional data. After having constructed an OF prediction rule using a training dataset, 
#' it can be used to predict the values of the ordinal target variable for new observations.
#' Moreover, by means of the (permutation-based) variable importance measure of OF, it is also possible to rank the covariates 
#' with respect to their importances in the prediction of the values of the ordinal target 
#' variable. \cr
#' OF is presented in Hornung (2020).
#'
#' Starting with package version 2.4, it is also possible to obtain class probability 
#' predictions in addition to the class point predictions and variable importance values
#' based on the class probabilities through using the (negative) ranked probability score (Epstein, 1969)
#' as performance function (\code{perffunction="probability"}, new default). Using the ranked probability score in the variable importance can be expected to deliver more stable variable rankings, because the ranked probability score accounts for the ordinal scale of the dependent variable.  In situations in which there is no need for predicting class probabilities, but simply
#' class predictions are sufficient, other performance functions may be more suitable. See the documentation of the \code{\link{ordfor}} function for further details.
#'
#' For a brief, practice-orientated introduction to OF see: \code{\link{ordfor}}
#'
#' The main functions are: \code{\link{ordfor}} (construction of OF prediction rules) and
#' \code{\link{predict.ordfor}} (prediction of the values of the target variable values of new observations).
#'
#' NOTE: \pkg{ordinalForest} uses R code and C++ code from the R package \pkg{ranger} for the involved regression forests.
#' \pkg{ordinalForest} does, however, not depend on \pkg{ranger} or import \pkg{ranger}, because it was necessary to
#' copy the C++ code and parts of the R code from \pkg{ranger} to \pkg{ordinalForest} instead. The reason for this
#' is that \pkg{ranger}'s C++ code had to be altered in part in order to implement ordinal forest.
#'
#' @references
#' \itemize{
#'   \item Hornung R. (2020) Ordinal Forests. Journal of Classification 37, 4–17. <\doi{10.1007/s00357-018-9302-x}>.
#'   \item Epstein E.S. (1969) A scoring system for probability forecasts of ranked categories, Journal of Applied Meteorology. 8(6), 985-987.
#'   }
#' 
#' @examples
#' \dontrun{
#' # Illustration of the key functionalities of the package:
#' ##########################################################
#' 
#' # Load example dataset:
#' 
#' data(hearth)
#' 
#' # Inspect the data:
#' table(hearth$Class)
#' dim(hearth)
#'
#' head(hearth) 
#'
#' 
#' # Split into training dataset and test dataset:
#' 
#' set.seed(123)
#' trainind <- sort(sample(1:nrow(hearth), size=floor(nrow(hearth)*(2/3))))
#' testind <- setdiff(1:nrow(hearth), trainind)
#' 
#' datatrain <- hearth[trainind,]
#' datatest <- hearth[testind,]
#' 
#' 
#' # Construct OF prediction rule using the training dataset (default 
#' # perffunction = "probability" corresponding to the 
#' # (negative) ranked probability score as performance function):
#' 
#' ordforres <- ordfor(depvar="Class", data=datatrain, nsets=1000, ntreeperdiv=100, 
#'   ntreefinal=5000, perffunction = "equal")
#' ordforres
#' 
#' # Study variable importance values:
#' sort(ordforres$varimp, decreasing=TRUE)
#'
#' # Take a closer look at the top variables:
#' boxplot(datatrain$oldpeak ~ datatrain$Class, horizontal=TRUE)
#' fisher.test(table(datatrain$exang, datatrain$Class))
#' 
#' # Predict values of the ordinal target variable in the test dataset:
#' 
#' preds <- predict(ordforres, newdata=datatest)
#' preds
#' 
#' # Compare predicted values with true values:
#' table(data.frame(true_values=datatest$Class, predictions=preds$ypred))
#' } 
#'
#' @name ordinalForest-package
#' @aliases ordinalForest
NULL