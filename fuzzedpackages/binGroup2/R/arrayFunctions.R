###############################################################################
# This file contains functions to calculate the expected testing expenditure 
# of array testing, and the individual specific measures of accuracy Pooling 
# Sensitivity, Pooling Specificity, Pooling Positive Predicitive Value, and 
# Pooling Negative Predicitive Value.
# Last modified date: 8-15-2011
# Author: Chris McMahan
###############################################################################
# Brianna Hitt - 12-18-19
# Revised to correct the calculation of p.r and p.c, and 
# to allow for Se/Sp to vary across stages of testing

Array.Measures<-function(p,se,sp){
  d<-dim(p)
  J<-d[1]
  K<-d[2]
  
  #######################################################
  # Finds probability that each individual tests positive
  e.ind<-matrix(-100,nrow=J,ncol=K)
  r<-matrix(0:J,ncol=1,nrow=(J+1))
  RA<-apply(r,1,rc.pos,p=p,id=1)
  c<-matrix(0:K,ncol=1,nrow=(K+1))
  CB<-apply(c,1,rc.pos,p=p,id=2)
  
  for(j in 1:J){
    for(k in 1:K){
      
      p.temp1<-p
      p.temp1[,k]<-rep(0,J)
      RAJ<-apply(r,1,rc.pos,p=p.temp1,id=1)
      # below is Chris McMahan's original code
      # p.c<-prod((1-p[,k]))
      # Brianna Hitt - 12-18-19
      # below is the corrected code, where p.c is a product over the rows
      # see McMahan et al. (2012) for details
      p.c <- prod((1-p[j,]))
      # P(all Rj = 0, Ck = 1)
      # below is Chris McMahan's original code
      # a1<-sum((1-sp)*(1-se)^r*sp^(J-r)*p.c*RAJ + 
      #           se*(1-se)^r*sp^(J-r)*((RA-p.c*RAJ)))
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      a1 <- sum((1-sp[1])*(1-se[1])^r*sp[1]^(J-r)*p.c*RAJ + 
                  se[1]*(1-se[1])^r*sp[1]^(J-r)*((RA-p.c*RAJ)))
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # a1 <- sum((1-sp[2])*(1-se[1])^r*sp[1]^(J-r)*p.c*RAJ + 
      #             se[2]*(1-se[1])^r*sp[1]^(J-r)*((RA-p.c*RAJ)))
            
      p.temp1<-p
      p.temp1[j,]<-rep(0,K)
      CBK<-apply(c,1,rc.pos,p=p.temp1,id=2)
      # below is McMahan's original code
      # p.r<-prod((1-p[j,]))
      # Brianna Hitt - 12-18-19
      # below is the corrected code, where p.r is a product over the columns
      # see McMahan et al. (2012) for details
      p.r <- prod((1-p[,k]))
      # P(Rj = 1, all Ck = 0)
      # below is Chris McMahan's original code
      # a2<-sum((1-sp)*(1-se)^c*sp^(K-c)*p.r*CBK + 
      #           se*(1-se)^c*sp^(K-c)*((CB-p.r*CBK)))
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      a2 <- sum((1-sp[1])*(1-se[1])^c*sp[1]^(K-c)*p.r*CBK + 
                  se[1]*(1-se[1])^c*sp[1]^(K-c)*((CB-p.r*CBK)))
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # a2 <- sum((1-sp[1])*(1-se[2])^c*sp[2]^(K-c)*p.r*CBK + 
      #             se[1]*(1-se[2])^c*sp[2]^(K-c)*((CB-p.r*CBK)))
      
      # P(Rj = 1, Ck = 1)
      # below is Chris McMahan's original code
      # a3<-se^2+(1-se-sp)^2*prod(1-p[,k])*prod(1-p[j,])/(1-p[j,k])+
      #   (se*(1-sp)-se^2)*(prod(1-p[j,])+prod(1-p[,k]))
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      a3 <- se[1]^2 + (1-se[1]-sp[1])^2*prod(1-p[,k])*prod(1-p[j,])/(1-p[j,k]) + 
        (se[1]*(1-sp[1])-se[1]^2)*(prod(1-p[j,])+prod(1-p[,k]))
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # a3 <- se[1]*se[2] + (se[2]-se[2]*sp[1])*prod(1-p[,k]) +
      #   (se[1]-se[1]*sp[2])*prod(1-p[j,]) -
      #   se[1]*se[2]*(prod(1-p[j,])+prod(1-p[,k])) +
      #   ((1-sp[1])*(1-sp[2])-(se[1]-se[1]*sp[2])-(se[2]-se[2]*sp[1])+
      #      se[1]*se[2])*(prod(1-p[,k])*prod(1-p[j,])/(1-p[j,k]))
      
      e.ind[j,k]<-(a1+a2+a3)
    }
  }
  # Brianna Hitt - 04.02.2020
  # changed from "T" to "ET"
  ET <- sum(e.ind)+J+K
  
  ###############################################
  # Finds Pooling Sensitivity for each individual
  pse.ind<-matrix(-100,nrow=J,ncol=K)
  # below is Chris McMahan's original code
  # c0<-1-(se+(1-se-sp)*apply((1-p),2,prod))
  # r0<-1-(se+(1-se-sp)*apply((1-p),1,prod))
  # below assumes that se, sp for row and column testing are equal, 
  #  where se, sp is specified as c(row/column testing, individual)
  # P(Ck' = 0)
  c0<-1-(se[1]+(1-se[1]-sp[1])*apply((1-p),2,prod))
  # P(Rj' = 0)
  r0<-1-(se[1]+(1-se[1]-sp[1])*apply((1-p),1,prod))
  # below allows se, sp for row and column testing to differ,
  #   where se, sp is specified as c(row testing, column testing, 
  #   individual testing)
  # c0<-1-(se[2]+(1-se[2]-sp[2])*apply((1-p),2,prod))
  # r0<-1-(se[1]+(1-se[1]-sp[1])*apply((1-p),1,prod))
  
  for(j in 1:J){
    for(k in 1:K){
      # below is Chris McMahan's original code
      # pse.ind[j,k]<-se^3+se^2*(1-se)*(prod(c0[-k])+prod(r0[-j]))
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      pse.ind[j,k] <- se[2]*se[1]^2 + 
        se[2]*se[1]*(1-se[1])*(prod(c0[-k])+prod(r0[-j]))
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # pse.ind[j,k]<-se[3]*se[2]*se[1] + se[3]*se[1]*(1-se[2])*prod(c0[-k]) + 
      #   se[3]*se[2]*(1-se[1])*prod(r0[-j])
    }
  }
  
  ###############################################
  # Finds Pooling Specificity for each individual
  psp.ind<-matrix(-100,nrow=J,ncol=K)
  r<-matrix(0:J,ncol=1,nrow=(J+1))
  c<-matrix(0:K,ncol=1,nrow=(K+1))
  
  for(j in 1:J){
    for(k in 1:K){
      p.temp1<-p
      p.temp1[,k]<-rep(0,J)
      p.temp2<-p
      p.temp2[j,k]<-0
      
      RAJ<-apply(r,1,rc.pos,p=p.temp1,id=1)
      RAj<-apply(r,1,rc.pos,p=p.temp2,id=1)
      # this is McMahan's original code
      # p.c<-prod((1-p.temp2[,k]))
      # Brianna Hitt - 12-18-19
      # below is the corrected code, where p.c is a product over the rows
      p.c <- prod((1-p.temp2[j,]))
      # P(Yjk = 1, all Rjk = 0, Ck = 1 | ~Yjk = 0) / (1 - Sp:I) = 
      #   P(all Rj = 0, Ck = 1 | ~Yjk = 0)
      # below is Chris McMahan's original code
      # a1<-sum((1-sp)*(1-se)^r*sp^(J-r)*p.c*RAJ + 
      #           se*(1-se)^r*sp^(J-r)*((RAj-p.c*RAJ)))
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      a1 <- sum((1-sp[1])*(1-se[1])^r*sp[1]^(J-r)*p.c*RAJ + 
                  se[1]*(1-se[1])^r*sp[1]^(J-r)*((RAj-p.c*RAJ)))
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # a1 <- sum((1-sp[2])*(1-se[1])^r*sp[1]^(J-r)*p.c*RAJ + 
      #             se[2]*(1-se[1])^r*sp[1]^(J-r)*((RAj-p.c*RAJ)))
      
      p.temp1<-p
      p.temp1[j,]<-rep(0,K)
      p.temp2<-p
      p.temp2[j,k]<-0
      
      CBK<-apply(r,1,rc.pos,p=p.temp1,id=2)
      CBk<-apply(r,1,rc.pos,p=p.temp2,id=2)
      # this is McMahan's original code
      # p.r<-prod((1-p.temp2[j,]))
      # Brianna Hitt - 12-18-19
      # below is the corrected code, where p.c is a product over the rows
      p.r <- prod((1-p.temp2[,k]))
      # P(Yjk = 1, Rj = 1, all Ck = 0 | ~Yjk = 0) / (1 - Sp:I) = 
      #   P(Rj = 1, all Ck = 0 | ~Yjk = 0)
      # below is Chris McMahan's original code
      # a2<-sum((1-sp)*(1-se)^c*sp^(K-c)*p.r*CBK + 
      #           se*(1-se)^c*sp^(K-c)*((CBk-p.r*CBK)))
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      a2 <- sum((1-sp[1])*(1-se[1])^c*sp[1]^(K-c)*p.r*CBK + 
                  se[1]*(1-se[1])^c*sp[1]^(K-c)*((CBk-p.r*CBK)))
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # a2 <- sum((1-sp[1])*(1-se[2])^c*sp[2]^(K-c)*p.r*CBK + 
      #             se[1]*(1-se[2])^c*sp[2]^(K-c)*((CBk-p.r*CBK)))
      
      # P(Yjk = 1, Rj = 1, Ck = 1 | ~Yjk = 0) / (1 - Sp:I)
      # below is Chris McMahan's original code
      # a3<-(se+(1-se-sp)*prod((1-p[,k]))/(1-p[j,k]))*
      #   (se+(1-se-sp)*prod((1-p[j,]))/(1-p[j,k]))
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      a3 <- (se[1]+(1-se[1]-sp[1])*prod((1-p[,k]))/(1-p[j,k]))*
        (se[1]+(1-se[1]-sp[1])*prod((1-p[j,]))/(1-p[j,k])) 
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # a3 <- (se[1]+(1-se[1]-sp[1])*prod((1-p[,k]))/(1-p[j,k]))*
      #   (se[2]+(1-se[2]-sp[2])*prod((1-p[j,]))/(1-p[j,k]))

      # PSp = 1 - (1 - Sp:I)*(a1 + a2 + a3)
      # below is Chris McMahan's original code
      # psp.ind[j,k]<-1-(1-sp)*(a1+a2+a3)
      # below assumes that se, sp for row and column testing are equal, 
      #  where se, sp is specified as c(row/column testing, individual)
      psp.ind[j,k]<-1-(1-sp[2])*(a1+a2+a3)
      # below allows se, sp for row and column testing to differ,
      #   where se, sp is specified as c(row testing, column testing, 
      #   individual testing)
      # psp.ind[j,k]<-1-(1-sp[3])*(a1+a2+a3)
    }
  }
  
  #############################################################################
  # Finds Pooling Positive (Negative) Predicitive Value for each individual
  ppv.ind<-p*pse.ind/(p*pse.ind+(1-p)*(1-psp.ind))
  npv.ind<-(1-p)*psp.ind/((1-p)*psp.ind + p*(1-pse.ind))
  
  # Brianna Hitt - 04.02.2020
  # changed "T" to "ET" in returned list of values
  return(list("ET"=ET,"PSe"=pse.ind,"PSp"=psp.ind, "PPV"=ppv.ind, "NPV"=npv.ind))
}

##################################################################
##################################################################
# Support Functions needed for Array.Measures()               ####
##################################################################
##################################################################

rc.pos<-function(n,p,id){
  d<-dim(p)
  M<-d[id]
  p.pos<-1-apply((1-p),id,prod)
  p.frac<-p.pos/(1-p.pos)
  if(n==0){
    res<-prod(1-p.pos)
  }
  if(n>0){
    if(n==1){
      res<-sum(p.frac)*prod(1-p.pos)
    }
    if(n>1){
      temp<-cumsum(p.frac[M:n])
      temp<-temp[(M-n+1):1]
      if(n>2){
        for(i in 1:(n-2)){
          temp<-p.frac[(n-i):(M-i)]*temp
          temp<-cumsum(temp[(M-n+1):1])
          temp<-temp[(M-n+1):1]
        }
      }
      temp<-sum(p.frac[1:(M-n+1)]*temp)
      res<-temp*prod(1-p.pos)
    }
  }
  return(res)
}


###############################################################################
###############################################################################
###### Spiral or Gradient Array Construction Function                 #########
###############################################################################
###############################################################################

#' @title Arrange a matrix of probabilities for informative array testing
#'
#' @description Arrange a vector of individual risk probabilities in a matrix 
#' for informative array testing without master pooling.
#'
#' @param prob.vec vector of individual risk probabilities, of length nr*nc.
#' @param nr number of rows in the array.
#' @param nc number of columns in the array.
#' @param method character string defining the method to be used for matrix 
#' arrangement. Options include spiral ("\kbd{sd}") and gradient ("\kbd{gd}") 
#' arrangement. See McMahan et al. (2012) for additional details.
#'
#' @return A matrix of probabilities arranged according to the specified 
#' method.
#' 
#' @author This function was originally written by Christopher McMahan for 
#' McMahan et al. (2012). The function was obtained from 
#' \url{http://chrisbilder.com/grouptesting}.
#'
#' @references
#' \insertRef{McMahan2012b}{binGroup2}
#'
#' @seealso
#' \code{\link{expectOrderBeta}} for generating a vector of individual risk 
#' probabilities.
#'
#' @examples
#' # Use the gradient arrangement method to create a matrix
#' #   of individual risk probabilities for a 10x10 array.
#' # Depending on the specified probability, alpha level,
#' #   and overall group size, simulation may be necessary 
#' #   in order to generate the vector of individual 
#' #   probabilities. This is done using the expectOrderBeta() 
#' #   function and requires the user to set a seed in order 
#' #   to reproduce results.
#' set.seed(1107)
#' p.vec1 <- expectOrderBeta(p=0.05, alpha=2, grp.sz=100)
#' informativeArrayProb(prob.vec=p.vec1, nr=10, nc=10, 
#'                        method="gd")
#'
#' # Use the spiral arrangement method to create a matrix
#' #   of individual risk probabilities for a 5x5 array.
#' set.seed(8791)
#' p.vec2 <- expectOrderBeta(p=0.02, alpha=0.5, grp.sz=25)
#' informativeArrayProb(prob.vec=p.vec2, nr=5, nc=5, 
#'                        method="sd")

informativeArrayProb<-function(prob.vec, nr, nc, method = "sd"){
  
  prob.vec<-sort(prob.vec,decreasing=TRUE)
  if(method=="sd"){
    array.probs<-prob.vec[1:2]
    prob.vec<-prob.vec[-(1:2)]
    array.probs<-cbind(array.probs,sort(prob.vec[1:2],decreasing=FALSE))
    prob.vec<-prob.vec[-(1:2)]
    
    max.iter<-max(nr,nc)
    
    for(i in 1:max.iter){
      if(nrow(array.probs) < nr){
        array.probs<-rbind(array.probs,prob.vec[1:(ncol(array.probs))])
        prob.vec<-prob.vec[-(1:(ncol(array.probs)))]
      }
      if(ncol(array.probs) < nc){
        array.probs<-cbind(array.probs,sort(prob.vec[1:(nrow(array.probs))],
                                            decreasing=FALSE))
        prob.vec<-prob.vec[-(1:(nrow(array.probs)))]
      }
    }
  }
  
  if(method=="gd"){
    array.probs<-matrix(prob.vec,ncol=max(nr,nc),nrow=min(nr,nc),byrow=FALSE)
  }
  return(array.probs)
}
