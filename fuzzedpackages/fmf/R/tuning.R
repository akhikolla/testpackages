#'Tuning For Fast Class Noise Detector with Multi-Factor-Based Learning
#'
#' This function tunes the hyper-parameters the threshold and the k of k-NN
#' @param formula a formula describing the classification variable and the attributes to be used.
#' @param data,x data frame containing the tranining dataset to be filtered.
#' @param knn_k range of the total number of nearest neighbors to be used.The default is 3:5.
#' @param classColumn positive integer indicating the column which contains the
#' (factor of) classes. By default, a dataframe built from 'data' using the variables indicated in 'formula' and The first column is the response variable, thus no need to define the classColumn.
#' @param boxplot_range range of box and whisker diagram. The default is seq(0.8,1.2,0.1).
#' @param repeats the number of cross-validation. The default is 10.
#' @param method the classifier to be used to compute the accuracy. The valid methods are svm (default) and c50.
#' @param iForest compute iForest score or not. The dafault is TRUE.
#' @param threads the number of cores to be used in parallel
#' @param ... Optional parameters to be passed to other methods.
#' @return An object of class \code{filter}, which is a list with two components:
#' \itemize{
#'    \item \code{summary} is the a vector of values when different hyper-parameter is set.
#'    \item \code{call} contains the original call to the filter.
#' }
#' @author Wanwan Zheng
#' @import caret
#' @import C50
#' @importFrom("stats", "as.formula", "model.frame")
#' @importFrom("kernlab", "ksvm", "predict")
#' @examples
#' \donttest{
#' data(iris)
#' out = tuning(Species~.,iris)
#' }
#' @name tuning
#'
#' @export
tuning = function(x, ...){
  UseMethod("tuning")
}

#' @rdname tuning
#' @export
tuning.formula = function(formula,
                         data,
                           ...){
  if(!is.data.frame(data)){
    stop("data argument must be a data.frame")
  }
  modFrame = model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
  attr(modFrame,"terms") = NULL

  ret = tuning.default(x = modFrame,...,classColumn = 1)
  ret$call =  match.call(expand.dots = TRUE)
  ret$call[[1]] = as.name("tuning")
  return(ret)
}

#' @rdname tuning
#' @export
tuning.default = function(x,
                          knn_k=seq(3,7,2),
                          classColumn=1,
                          boxplot_range=seq(0.1,1.1,0.2),
                          repeats = 10,
                          method = "svm",
                          iForest = TRUE,
                          threads = 1,
                           ...)
{
  if(!is.data.frame(x)){
    stop("data argument must be a data.frame")
  }
  if(!classColumn%in%(1:ncol(x))){
    stop("class column out of range")
  }
  formu = as.formula(paste(names(x)[classColumn],"~.",sep = ""))
  cv = caret::createFolds(1:nrow(x),repeats)
  method = match.arg(method, c('svm', 'c50'), FALSE)
  summary =  matrix(NA, nrow = 0, ncol = 6)
  colnames(summary) = c("k","range","noiseNum","noiseScore","accuracy","kappa")

  for (kv in knn_k) {
    for(rv in boxplot_range){
      num.res = score.res = acc.res = kappa.res = as.numeric()
      for (v in 1:repeats) {
        train.t = x[-cv[[v]],]
        test.t  = x[cv[[v]],]
        output = fmf(formu,train.t, knn = kv, boxplot_range = rv, iForest = iForest, threads = threads)
        score.res[v] = mean(output$noise_score[output$remIdx])
        num.res[v] = length(output$remIdx)
         switch(method,
           svm = {svm.res = kernlab::ksvm(formu, output$cleanData);acc = caret::confusionMatrix(table(test.t[,1],kernlab::predict(svm.res,test.t[,-1])));
           W = c(acc$overall[1],acc$overall[2])},
           c50 = {c50.res = C50::C5.0(formu, output$cleanData);acc = caret::confusionMatrix(table(test.t[,1],predict(c50.res,test.t[,-1])));
           W = c(acc$overall[1],acc$overall[2])}
        )
        acc.res[v] = W[1];kappa.res[v] = W[2];
      }
      summary= rbind(summary,c(kv,rv,mean(num.res),mean(score.res),mean(acc.res),mean(kappa.res)))
    }
  }

  ##### Building the 'filter' object ###########
  #message("result:", "\n")
  #message(summary)
  call = match.call()
  call[[1]] = as.name("tuning")
  ret = list(summary = summary,
              call = call)
  class(ret) = "filter"
  return(ret)
}
