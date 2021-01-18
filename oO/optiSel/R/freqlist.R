
"freqlist"<-function(...){
  dots <- list(...)
  for(i in c("map", "Indiv", "thisBreed", "minSNP", "minL", "unitL")){
    attributes(dots)[[i]]<-attributes(dots[[1]])[[i]]
    for(k in 1:length(dots)){
      if(!identical(attributes(dots[[k]])[[i]], attributes(dots)[[i]])){
        stop(paste("Attribute", i , "is different for component 1 and component",k,".\n"))
      }
      attributes(dots[[k]])[[i]]<-NULL
    }
  }

  if(is.null(names(dots))){names(dots)<-rep("",length(dots))}
  Names <- names(dots)
  for(k in 1:length(dots)){
    if(is.null(attributes(dots[[k]])[["refBreeds"]])){
      stop(paste("Attribute 'refBreeds' is missing for component",k,"\n"))
    }
    if(Names[k]==""){Names[k] <- paste(attributes(dots[[k]])[["refBreeds"]],collapse=".")}
  }
  names(dots) <- Names
  for(k in 1:length(dots)){
    if(!("freq" %in% names(dots[[k]]))){stop("Component 'freq' is missing.")}
    dots[[k]]<-dots[[k]]$freq
  }
  class(dots) <- "HaploFreq"
  dots
}