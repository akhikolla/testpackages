lavaan.model <- function(y.name, x.name, med.name, medcov)
  {
      n.med <- length(med.name)
  
    mod <- paste(y.name, " ~ ", x.name, "\n", sep="")
      
    for(i in 1:n.med){
      mod <- paste(mod, y.name, " ~ ", med.name[i], " \n ", sep="")
    }

    for(i in 1:n.med){
      mod <- paste(mod, med.name[i], " ~ ", x.name, "\n ", sep="")
    }
    

    for(i in 1:n.med){
      mod <- paste(mod, med.name[i], "  ~~ ",round(medcov[i,i],3), " * ", med.name[i], " \n ", sep="")
    }
    for(i in 1:(n.med-1)){
      for(j in (i+1):n.med){
        mod <- paste(mod, med.name[i],"  ~~ ",round(medcov[i,j],3), " * ", med.name[j], " \n ", sep="")
      }
    }
    return(mod)
  }
