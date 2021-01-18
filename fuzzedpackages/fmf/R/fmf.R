#'Fast Class Noise Detector with Multi-Factor-Based Learning
#'
#' This function computes the noise score for each observation
#' @param formula a formula describing the classification variable and the attributes to be used.
#' @param data,x data frame containing the tranining dataset to be filtered.
#' @param knn total number of nearest neighbors to be used.The default is 5.
#' @param classColumn positive integer indicating the column which contains the
#' (factor of) classes. By default, a dataframe built from 'data' using the variables indicated in 'formula' and The first column is the response variable, thus no need to define the classColumn.
#' @param boxplot_range range of box and whisker diagram. The dafault is 1.
#' @param iForest compute iForest score or not. The dafault is TRUE.
#' @param threads the number of cores to be used in parallel.
#' @param ... optional parameters to be passed to other methods.
#' @return an object of class \code{filter}, which is a list with four components:
#' \itemize{
#'    \item \code{cleanData} is a data frame containing the filtered dataset.
#'    \item \code{remIdx} is a vector of integers indicating the indexes for
#'    removed instances (i.e. their row number with respect to the original data frame).
#'    \item \code{noise_score} is a vector of values indicating the optential of being a noise.
#'    \item \code{call} contains the original call to the filter.
#' }
#' @author Wanwan Zheng
#' @import solitude
#' @import Rcpp
#' @import RcppArmadillo
#' @import e1071
#' @import caret
#' @importFrom("stats", "as.formula", "model.frame")
#'
#' @examples
#' 
#' data(iris)
#' out = fmf(Species~.,iris)
#' 
#'@name fmf
#'@export


fmf = function(x, ...){
  UseMethod("fmf")
}

#' @rdname fmf
#' @export
fmf.formula = function(formula,
                       data,
                       ...){
  if(!is.data.frame(data)){
    stop("data argument must be a data.frame")
  }
  modFrame = model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
  attr(modFrame,"terms") = NULL
  ret = fmf.default(x=modFrame,...,classColumn=1)
  ret$call = match.call(expand.dots = TRUE)
  ret$call[[1]] = as.name("fmf")
  cleanData = data
  ret$cleanData = cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
  return(ret)
}

#' @rdname fmf
#' @export
fmf.default = function(x,
                       knn=5,
                       classColumn=1,
                       boxplot_range=1,
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
  X = x[,-classColumn]
  Y = as.numeric(x[,classColumn])
  distance = c('euclidean', 'canberra', 'mahalanobis','schi',"pearson_correlation")
  list_idx = list_dis = as.list(NULL)
  error = c()

  for(dd in 1:length(distance)){
    out = kernelknn(X, TEST_data = NULL, Y, k = nrow(X)-1, method = distance[dd],
                    threads = threads, regression = FALSE, Levels=Y)
    pre = apply(out$idx[,1:knn],1,function(x){
      tab = table(x)
      names(tab)[match(max(tab),tab)]
    })
    list_idx = c(list_idx, list(out$idx))
    list_dis = c(list_dis, list(out$dis))
    error = c(error,caret::confusionMatrix(table(Y,as.numeric(pre)))$overall[1])
  }
  error = replace(error, which(is.na(error)), 0)
  for(i in 1:length(distance)){
    message(distance[i], ":", round(error[i],4))
  }
  selected = match(max(error),error)
  message("********************************************")
  message("Best index:", distance[selected])
  message("********************************************")
  Mat = list_idx[[selected]]
  Dis = list_dis[[selected]]

  Tab_Y = as.matrix(table(Y))
  len_min = min(Tab_Y)
  Mat_lab = cbind(Y, Mat)
  Rownames = as.integer(rownames(Tab_Y))

  list_score = sapply(1:nrow(Tab_Y), function(i, iforest = iForest){
    Log = Mat_lab[,1]==Rownames[i]
    subdata_idx = Mat_lab[Log,1:len_min]
    subdata = as.data.frame(cbind(subdata_idx, Dis[Log,1:(len_min-1)]))
    real = rep(Rownames[i],len_min-1)

    score = apply(subdata,1, function(x, true=real,length=len_min){
      pre = x[2:length]
      log = pre[pre!=true]
      Len= length(log)/(length-1)
      if((Len == 0) ||(Len == 1)){
        Len
      }else{
        Len * (1-entropy(table(as.numeric(log))))
      }
    })
    score = normalized(score)

    density = apply(subdata,1, function(x, true=real,length=len_min){
      pre = x[2:length]
      den = x[-c(1:length)]
      Len= length(pre[pre!=true])

      if(Len == 0){
        0
      }else if (Len == length-1){
        1
      }else{
        density_hit = mean(den[pre==true])
        density_diff = mean(den[pre!=true])
        1-e1071::sigmoid(density_diff/density_hit)
      }
    })
    density = normalized(density)

    if(iforest == TRUE){
      isf = solitude::isolationForest$new(sample_size=floor(nrow(X[Log,])/3))
      isf$fit(X[Log,])
      iForest.score = normalized(isf$predict(X[Log,])$anomaly_score)
      return((score+density+iForest.score)/3)
    }else{
      return((score+density)/2)
    }
  })

  if(is.matrix(list_score)){
    list_score = as.vector(list_score)
  }else{
    list_score = unlist(list_score)
  }

  #' @keywords internal
  #' @importFrom("graphics", "boxplot")

  outlier = unique(boxplot(list_score,plot=FALSE,range = boxplot_range)$out)
  isNoise = logical(length(list_score))
  if(length(outlier) > 0){
    for(i in 1:length(outlier)){
      isNoise[!isNoise] = (list_score[!isNoise] == outlier[i])
    }
  }

  ##### Building the 'filter' object ###########
  remIdx = which(isNoise)
  #if(length(remIdx)>0){
    #cat("noise:",remIdx, "\n")
  #}
  cleanData = x[!isNoise,]
  call = match.call()
  call[[1]] = as.name("fmf")
  ret = list(cleanData = cleanData,
             remIdx = remIdx,
             noise_score = list_score,
             call = call)
  class(ret) = "filter"
  return(ret)
}
