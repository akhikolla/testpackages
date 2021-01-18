# function to conduct k-fold cross validation.

.Mesh <- function(rownamesY, kFold)
{
  numSamples <- length(rownamesY)
  res <- NULL
  subSampleSize <- floor(numSamples/kFold)
  for (i in 1:kFold)
  {
    start <- (i-1)*subSampleSize + 1
    if(i < kFold)
      end <- i*subSampleSize
    else
      end <- numSamples
    if(i == 1)
      res <- list(c(start:end))
    else
      res[[i]] <- c(start:end) 
  }
  res
}
.whichpart <- function(x, n=30) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}

## example function for return the top 4 hit base on P value for each feature.
.selectCrit <- function(J, P)
{
  pcut <- 0.5
  hit <- NULL
  for(i in 1:ncol(P)){
    if( i == 1)
      hit <- list(rownames(P)[.whichpart(P[,i], 4)])
    else
      hit[[i]] = rownames(P)[.whichpart(P[,i], 4)]
  }
  res <- do.call(cbind, hit)
  colnames(res) <- colnames(P)
  res
}

fastJT.select <- function(Y, X, cvMesh = NULL, kFold = 10L, selCrit = NULL,  outTopN = 15L, numThreads = 1L)
{
  # get the number of the sample
  numSamples <- nrow(Y)
  
  # get the ids to remove to form each subsample.
  if(is.null(cvMesh))
    subRmIds <- .Mesh(rownames(Y),kFold)		
  else
    subRmIds <- cvMesh(rownames(Y),kFold)
  
  # run test on the subsamples and produce return.
  res <- NULL
  pvalue <- NULL
  hits <- NULL
  if(!is.null(selCrit)){
    for (i in 1:length(subRmIds)){
      if(i == 1){
        res <- list(fastJT(Y[-subRmIds[[1]],], X[-subRmIds[[1]],],
                           NA, numThreads, standardized = TRUE))
        pvalue <- list(pvalues(res[[i]]))
        hits <- list(selCrit(res[[i]]$J, pvalue[[i]]))
        res[[i]] <- res[[i]]$J
      }else{    
        res[[i]] <- fastJT(Y[-subRmIds[[i]],], X[-subRmIds[[i]],],
                           NA, numThreads, standardized = TRUE)
        pvalue[[i]] <- pvalues(res[[i]])
        hits[[i]] <- selCrit(res[[i]]$J, pvalue[[i]])
        res[[i]] <- res[[i]]$J
      }	
    }
  }
  
  if(is.null(selCrit)){
    for (i in 1:length(subRmIds)){
      if(i == 1){
        res <- list(fastJT(Y[-subRmIds[[1]],], X[-subRmIds[[1]],],
                           outTopN, numThreads, standardized = TRUE))
        pvalue <- list(pvalues(res[[i]]))
        hits <- list(res[[i]]$XIDs)
        res[[i]] <- res[[i]]$J
      }else{
        res[[i]] <- fastJT(Y[-subRmIds[[i]],], X[-subRmIds[[i]],],
                           outTopN, numThreads, standardized = TRUE)
        pvalue[[i]] <- pvalues(res[[i]])
        hits[[i]] <- res[[i]]$XIDs
        res[[i]] <- res[[i]]$J      
      } 
    }
  }
  
  result <- NULL
  
  # return sevaral lists containing all the results for the cross validation.
  
  result$J <- res
  result$Pval <- pvalue
  result$XIDs <- hits
  
  result
}





