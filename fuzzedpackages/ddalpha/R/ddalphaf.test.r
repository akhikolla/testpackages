ddalphaf.test <- function(learn, learnlabels, test, testlabels, disc.type = c("LS", "comp"), ...){
  ops <- options(warn = -1)
  on.exit(options(ops))
  
  disc.type <- match.arg(disc.type)
  
  ftrain = switch(disc.type,
              "LS" = ddalphaf.train,
              "comp" = compclassf.train
  )
  fclassify = switch(disc.type,
                  "LS" = ddalphaf.classify,
                  "comp" = compclassf.classify
  )
  
  
  tryCatch({
    time <- system.time(
      ddalpha <- ftrain(learn, learnlabels, ...)
    )
    cc = fclassify(objectsf = test,ddalphaf = ddalpha)
    if (is.numeric(testlabels[[1]])){
      if(is.factor(cc[[1]]) || is.character(cc[[1]])){
        cc <- unlist(lapply(cc, as.character))
        cc[cc == "Ignored"] <- NA
      }
      equal = (cc == testlabels)
    } else {
      cc <- unlist(lapply(cc, as.character)) 
      equal = (cc == as.character(testlabels))      
    }
    
    if(!(T %in% equal) && !(F %in% equal)) 
    { return(NA)}
    error = sum(!equal,na.rm = T)/(sum(!equal,na.rm = T)+sum(equal,na.rm = T))
    return(list(error = error, correct = sum(equal,na.rm = T), incorrect = sum(!equal,na.rm = T), 
                total = length(cc)-sum(is.na(equal)), ignored = sum(is.na(equal)), n = length(cc),
                time = time[1])) 
  }
  #  tryCatch({}
  , error = function(e) {
    print ("ERROR T")
    print (e)
  }, finally = {          
  })
  return (NA)
}


ddalphaf.getErrorRateCV <- function(dataf, labels, numchunks = 10, disc.type = c("LS", "comp"),  ...){
  n = length(dataf)
  numchunks = min(n, numchunks)
  chunksize = ceiling(n/numchunks)
  
  sample = seq(from = 1, by = numchunks, length.out = chunksize)   
  
  errors = 0
  total = 0
  times = c()
  
  for (i in 1:numchunks){
    sample = sample[sample<=n]
    learn = dataf[-sample]
    test  = dataf[sample]
    learnlabels = labels[-sample]
    testlabels  = labels[sample]
    
    el = ddalphaf.test(learn, learnlabels, test, testlabels, disc.type, ...)
    if(is.list(el)){
      errors = errors + el$incorrect
      total = total + el$total
      times = c(times,el$time)
    }
    
    sample = sample+1
  }
  
  return (list(errors = errors/total, time = mean(times), time_sd = sd(times)))  
}


ddalphaf.getErrorRatePart <- function(dataf, labels, size = 0.3, times = 10, disc.type = c("LS", "comp"),  ...){
  
  if (!is.numeric(size) || size <=0 || size >= length(dataf)) stop("Wrong size of excluded sequences")
  
  if(size < 1)
    size = max(1, size*length(dataf)) # at least 1 point
  
  size = as.integer(size)
  
  indexes = 1:length(dataf)
  
  errors = c()
  total = 0
  time = c()
  
  for (i in 1:times){
    samp = sample(indexes, size)
    learn = dataf[-samp]
    test  = dataf[samp]
    learnlabels = labels[-samp]
    testlabels  = labels[samp]
    
    el = ddalphaf.test(learn, learnlabels, test, testlabels, disc.type, ...)
    if(is.list(el)){
      errors = c(errors,el$incorrect/el$total)
      time = c(time,el$time)
    }
  }
  
  return (list(errors = mean(errors), errors_sd = sd(errors), errors_vec = errors, time = mean(time), time_sd = sd(time)))  
}