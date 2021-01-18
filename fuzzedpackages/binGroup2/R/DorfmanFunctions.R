###############################################################################
# Used to find the characteristics of an Informative Dorfman decoding process 
# specified by method= OD (Optimal Dorfman), TOD (Thresholded Optimal Dorfman), 
# and PSOD (Pool Specific Optimal Dorfman) given 
# p (a vector of all subjects' infection probabilities), 
# se (sensitivity of diagnostic test),
# sp (specificity of diagnostic test), 
# max.pool (maximum allowable pool size),
# thresh.pool (initial pool size used for TOD if threshold not specified), and 
# threshold (threshold value for TOD, note if a threshold value is not 
#   specified one is found algorithmically). 
# Function returns 
# p.star (threshold value used, only applies when using TOD),
# res.e (the expected expenditure of the decoding process),
# res.v (the variance of expenditure of the decoding process), and 
# res.mat a matrix of summary measures which includes each subject's 
#   infection probability , 
# pool (pool to which they belong),
# PSe (pooling sensitivity), 
# PSp (pooling specificity), 
# PPV (pooling positive predictive value), and 
# NPV (pooling negative predictive value).
#
#EXAMPLE: opt.info.dorf(prob=rbeta(1000,1,10), se = 1, sp = 1, method ="OD", 
#                       max.pool=15, thresh.pool=8, threshold=NULL)

opt.info.dorf<-function(prob, se = 1, sp = 1, method ="OD", max.pool=15, 
                        thresh.pool=8, threshold=NULL){
  
  # Saves original ordering
  ind.order<-order(prob)
  
  # Orders subjects, required under all Informative measures
  prob<-sort(prob)
  
  # Determines number of subjects being screened, and sets up vectors for 
  # storing summary measures. Also initializes the threshold p.star
  N<-length(prob)
  pool.id<-rep(-1000,N)
  PSe<-rep(-100,N)
  PSp<-rep(-100,N)
  PPV<-rep(-100,N)
  NPV<-rep(-100,N)
  p.star<-threshold
  
  
  # If method is TOD this finds the threshold value and divides the 
  # subjects into the high and low risk classes
  if(method=="TOD"){
    if(is.null(p.star)){
      p.star<-thresh.val.dorf(p=prob, psz=thresh.pool, se=se, sp=sp)
    }
    if(p.star < 1){
      N.high<-length(prob[prob > p.star])
      N<-N-N.high
    }
  }
  
  
  # Finds optimal pool size to be used with OD or, finds optimal pool 
  # size to be used to decode low risk class in TOD
  if(method=="OD" | method=="TOD"){
    if(N==0){psz<-NULL}
    if(N > 0){
      res.psz<-opt.pool.size(p=prob[1:N] ,max.p=max.pool, se=se, sp=sp)
      J<-ceiling(N/res.psz)
      rem<-N-(J-1)*res.psz
      if(rem!=0){
        psz<-c(rep(res.psz,(J-1)),rem)
      }  
      if(rem==0){
        psz<-rep(res.psz,J)
      }
    }
    if(method=="TOD"){
      if(p.star < 1){
        psz<-c(psz,rep(1,N.high))
        J<-length(psz)
        N<-N+N.high
      }
    }
  }
  
  
  # Finds pool sizes to be used with PSOD
  if(method=="PSOD"){
    psz<-pool.specific.dorf(p=prob , max.p=max.pool , se=se , sp=sp)
    J<-length(psz)
  }
  
  # Finds measures pool by pool
  psz<-c(psz,0)
  lower<-1
  upper<-psz[1]
  vec.e<-rep(-1000,J)
  vec.v<-rep(-1000,J)
  for(i in 1:J){
    p.pool<-prob[lower:upper]
    pool.id[lower:upper]<-rep(i,length(p.pool))
    
    res<-characteristics.pool(p=p.pool,se=se,sp=sp)
    vec.e[i]<-res$e
    vec.v[i]<-res$v
    
    res.acc<-accuracy.dorf(p=p.pool,se=se,sp=sp)
    PSe[lower:upper]<-res.acc$PSe
    PSp[lower:upper]<-res.acc$PSp
    PPV[lower:upper]<-res.acc$PPV
    NPV[lower:upper]<-res.acc$NPV
    
    lower<-1+upper
    upper<-upper+psz[i+1]
  }
  
  # Finds total expectation and variation
  res.e<-sum(vec.e)
  res.v<-sum(vec.v)
  
  # Returns all subjects to original ordering, along with their 
  # corresponding measures
  prob<-prob[order(ind.order)]
  pool.id<-pool.id[order(ind.order)]
  PSe<-PSe[order(ind.order)]
  PSp<-PSp[order(ind.order)]
  PPV<-PPV[order(ind.order)]
  NPV<-NPV[order(ind.order)]
  
  res.mat<-matrix(c(pool.id, prob, PSe, PSp, PPV, NPV),
                  nrow=N,ncol=6,byrow=FALSE, 
                  dimnames=list(as.character(1:N) , 
                                c("pool","probability","PSe", 
                                  "PSp", "PPV", "NPV")))
  
  prob<-prob[ind.order]
  return(list("tv"=p.star, "e"=res.e, "v"=res.v, "summary"=res.mat))
}



######################################
### DORFMAN DECODING FUNCTIONS #######
######################################

###############################################################################
# Function for the expectation and variation in testing expenditure 
# of a pool of size greater than or equal to one
# here p is a vector of all subjects probability of infection,
# se is sensitivity, and
# sp is specificity.
# res.e is the expected expenditure of the pool and
# res.v is the variance of expenditure of the pool.

characteristics.pool<-function(p,se,sp){
  n<-length(p)
  if(n>1){
    # below is Chris McMahan's original code
    # prob<-se+(1-se-sp)*prod(1-p)
    # Brianna Hitt - 01-03-20
    # below allows Se, Sp to vary across stages of testing
    prob<-se[1]+(1-se[1]-sp[1])*prod(1-p)
    res.e<-1+n*prob
    res.v<-n^2*prob*(1-prob)
  }
  if(n==1){
    res.e<-1
    res.v<-0
  }
  return(list("e"=res.e, "v"=res.v))
}



###############################################################################
# Function that returns PSe, PSp, PPV, and NPV for all individuals
# belonging to a pool of size greater than or equal to one
# here p is a vector of all subject probabilities,
# se is sensitivity, and
# sp is specificity.

accuracy.dorf<-function(p,se,sp){
  cj<-length(p)
  se.vec<-rep(-100,cj)
  sp.vec<-rep(-100,cj)
  ppv.vec<-rep(-100,cj)
  npv.vec<-rep(-100,cj)
  
  if(cj==1){
    # below is Chris McMahan's original code
    # se.vec<-se
    # sp.vec<-sp
    # ppv.vec<-p*se/(p*se+(1-p)*(1-sp))
    # npv.vec<-(1-p)*sp/((1-p)*sp + p*(1-se))
    
    # Brianna Hitt - 01-03-20
    # below allows Se, Sp to vary across stages of testing
    se.vec<-se[2]
    sp.vec<-sp[2]
    # Brianna Hitt - 01-30-20
    # revised code below to use se.vec and sp.vec for efficiency
    ppv.vec<-p*se.vec/(p*se.vec+(1-p)*(1-sp.vec))
    npv.vec<-(1-p)*sp.vec/((1-p)*sp.vec + p*(1-se.vec))
  }
  
  if(cj>1){
    for(i in 1:cj){
      # below is Chris McMahan's original code
      # se.vec[i]<-se^2
      # sp.vec[i]<-1-(1-sp)*(se+(1-se-sp)*prod(1-p[-i]))
      # ppv.vec[i]<-p[i]*se^2/(p[i]*se^2+(1-p[i])*(1-sp.vec[i]))
      # npv.vec[i]<-(1-p[i])*sp.vec[i]/((1-p[i])*sp.vec[i] + p[i]*(1-se^2))
      
      # Brianna Hitt - 01-30-20
      # below allows Se, Sp to vary across stages of testing
      se.vec[i]<-se[1]*se[2]
      sp.vec[i]<-1-(1-sp[2])*(se[1]+(1-se[1]-sp[1])*prod(1-p[-i]))
      # Brianna Hitt - 01-30-20
      # revised code below to use se.vec and sp.vec for efficiency
      ppv.vec[i]<-p[i]*se.vec[i]/(p[i]*se.vec[i]+(1-p[i])*(1-sp.vec[i]))
      npv.vec[i]<-(1-p[i])*sp.vec[i]/((1-p[i])*sp.vec[i] + p[i]*(1-se.vec[i]))
    }
  }
  
  return(list("PSe"=se.vec, "PSp"=sp.vec, "PPV"=ppv.vec, "NPV"=npv.vec))
}



###############################################################################
# Pooling Algorithm Minimizes Individual 
# Risk on a per pool basis, so it allows for multiple 
# pooling sizes (psz) here p is a vector of all subject probabilities,
# se is sensitivity,
# sp is specificity, and
# max.p is the maximum allowable pool size.

pool.specific.dorf<-function(p,max.p,se,sp){
  p<-sort(p)
  N<-length(p)
  psz<-rep(-99,N)
  k<-0
  ind<-1
  
  while(k!=N){
    et<-1
    max=min(length(p),max.p)
    if(length(p)>1){
      et<-c(et,rep(99,max))
      for(i in 2:max){
        # below is Chris McMahan's original code
        # et[i]<-((i)*(se-(se+sp-1)*prod(1-p[1:(i)]))+1)/(i)
        
        # Brianna Hitt - 01-03-20
        # below allows Se, Sp to vary across stages of testing
        et[i]<-((i)*(se[1]-(se[1]+sp[1]-1)*prod(1-p[1:(i)]))+1)/(i)
      }}
    m<-1:max
    m<-m[order(et)[1]]
    psz[ind]<-m
    ind<-ind+1
    p<-p[-(1:m)]
    k<-k+m
  }
  psz<-psz[psz!=-99]
  return(psz)
}



###############################################################################
# Function for finding the optimal pool size (psz),                                               
# here p is a vector of all subject probabilities,
# se is sensitivity,
# sp is specificity, and
# max.p is the maximum allowable pool size.

opt.pool.size<-function(p, max.p, se=1, sp=1){
  
  N<-length(p)
  M<-min(N, max.p)
  p<-sort(p)
  e0<-N
  psz<-2
  L<-ceiling(N/psz)
  if(N<(L*psz)){psz.vec<-c(rep(psz,L-1),(N-psz*(L-1)))}
  if(N==(L*psz)){psz.vec<-rep(psz,L)}
  e1<-0
  sub.ind<-c(0,cumsum(psz.vec))
  for(m in 1:L){
    # below is Chris McMahan's original code
    # e1<-e1+(1+psz.vec[m]*(se+(1-se-sp)*prod(1-p[(sub.ind[m]+1):sub.ind[m+1]])))
    
    # Brianna Hitt - 01-03-20
    # below allows Se, Sp to vary across stages of testing
    e1<-e1+(1+psz.vec[m]*(se[1]+(1-se[1]-sp[1])*prod(1-p[(sub.ind[m]+1):sub.ind[m+1]])))
  }
  
  while(e1<e0 & psz<=max.p){
    e0<-e1
    psz<-psz+1
    L<-ceiling(N/psz)
    if(N<(L*psz)){psz.vec<-c(rep(psz,L-1),(N-psz*(L-1)))}
    if(N==(L*psz)){psz.vec<-rep(psz,L)}
    e1<-0
    sub.ind<-c(0,cumsum(psz.vec))
    for(m in 1:L){
      # below is Chris McMahan's original code
      # e1<-e1+(1+psz.vec[m]*(se+(1-se-sp)*prod(1-p[(sub.ind[m]+1):sub.ind[m+1]])))

      # Brianna Hitt - 01-03-20
      # below allows Se, Sp to vary across stages of testing
      e1<-e1+(1+psz.vec[m]*(se[1]+(1-se[1]-sp[1])*prod(1-p[(sub.ind[m]+1):sub.ind[m+1]])))
    }
  }
  
  return(psz-1)
}



###############################################################################
# Thresholding function for Dorfman Procedure, finds threshold (p.star).                                                                   
# Here p is a vector of all subject probabilities,
# se is sensitivity,
# sp is specificity, and
# psz is the initial pool size.

thresh.val.dorf<-function(p, psz, se=1, sp=1){
  N<-length(p)
  p<-sort(p)
  L<-ceiling(N/psz)
  p.star<-1
  upper<-N
  lower<-upper-psz+1
  lower<-max(1,lower)
  # below is Chris McMahan's original code
  # Ex.pool<-length(p[lower:upper])*(se+(1-se-sp)*prod(1-p[lower:upper]))+1

  # Brianna Hitt - 01-03-20
  # below allows Se, Sp to vary across stages of testing
  Ex.pool<-length(p[lower:upper])*(se[1]+(1-se[1]-sp[1])*prod(1-p[lower:upper]))+1
  
  if(Ex.pool > length(p[lower:upper])){
    while(Ex.pool > length(p[lower:upper]) & lower > 1){
      upper<-upper-psz
      lower<-lower-psz
      lower<-max(1,lower)
      # below is Chris McMahan's original code
      # Ex.pool<-length(p[lower:upper])*(se+(1-se-sp)*prod(1-p[lower:upper]))+1

      # Brianna Hitt - 01-03-20
      # below allows Se, Sp to vary across stages of testing
      Ex.pool<-length(p[lower:upper])*(se[1]+(1-se[1]-sp[1])*prod(1-p[lower:upper]))+1
    }
    p.star<-(p[upper]+p[upper+1])/2
  }
  if(lower==1){p.star<-0}
  return(p.star)
  
}



