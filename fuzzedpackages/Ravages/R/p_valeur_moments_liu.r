p.valeur.moments.liu <- function(Q, mu, sigma, skewness, kurtosis, estimation.pvalue) {
  # l'astuce pour le cas kurtosis < 3 vient de SKAT
  # ca revient a  prendre une approximation normale
  if(!is.na(kurtosis)){ #tester si kurtosis a NA: peut etre le cas si seulement variants monomorphes
    if(kurtosis < 0){
      chi2val <- 1e5 + sqrt(2*1e5 / sigma**2)*(Q - mu)
      p.value <- pchisq(chi2val, 1e5, lower.tail=FALSE)
    }else{
      s1 <- skewness / sqrt(8)
      s2 <- kurtosis / 12
      if(s1^2 > s2){
        params <- Liu.sitA(s1, s2)
        if(params["l"] < 0 | params["d"] < 0){ #situation B
          params <- Liu.sitB(s1, s2, estimation.pvalue = estimation.pvalue)
         }      
      }else{
        params <- Liu.sitB(s1, s2, estimation.pvalue = estimation.pvalue)
      }
      muX <- params["l"]+params["d"]
      sigmaX <- sqrt(2)*params["a"]
      Q.norm <- (Q-mu)/sigma
      Q.norm1 <- Q.norm*sigmaX + muX
  
      #Calcul de la p-valeur 
      p.value <- pchisq(Q.norm1, df=params["l"], ncp=params["d"], lower.tail=FALSE)
    }
  }else{
    p.value <- NA
  }
  return(p.value)
}


Liu.sitA <- function(s1, s2){
  a <- 1/(s1-sqrt(s1**2-s2))
  d <- s1*a**3-a**2
  l=a**2-2*d
  return(c(a = as.numeric(a), d = as.numeric(d), l = as.numeric(l)))
}

Liu.sitB <- function(s1, s2, estimation.pvalue){
  if(estimation.pvalue == "skewness"){
    a=1/s1
    d=0
    l=1/s1**2
    if(d < 0 | l <0){ #solution de secours: use kurtosis
      l = 1/s2
      a = sqrt(l)
      d = 0
    }
  }else{
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  return(c(a = as.numeric(a), d = as.numeric(d), l = as.numeric(l)))
}
