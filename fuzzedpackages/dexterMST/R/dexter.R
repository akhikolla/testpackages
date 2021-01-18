

# internal dexter functions, slightly simplified if possible --------------
# copied from dexter to prevent notes
# TO~DO: check in tests that they always give equal results

# note: many of these can be removed when we remove the 0 categorie from dexter

dexter.IJ <- function(b, a, first, last, theta)
{
  nI=length(first)
  I = matrix(0,nI, length(theta))
  for (i in 1:nI)
  {
    aa = c(0, a[first[i]:last[i]])
    Fij = c(1, b[first[i]:last[i]])*exp(outer(aa,theta))
    Pi = apply(Fij,2,function(x)x/sum(x))
    M1 = Pi*aa
    M2 = M1*aa
    I[i,] =  colSums(M2) - colSums(M1)^2
  }
  colSums(I)
}

dexter.PairDIF_ <- function(par1,par2,cov1,cov2)
{
  labs=rownames(par1)
  D=kronecker(par2,t(par2),FUN="-")-kronecker(par1,t(par1),FUN="-") 
  var1=diag(cov1)
  var2=diag(cov2)
  S=(kronecker(var1,t(var1),FUN="+")-2*cov1)+(kronecker(var2,t(var2),FUN="+")-2*cov2)
  diag(S)=1
  D=D/sqrt(S)
  colnames(D)=labs; rownames(D)=labs
  return(D)
}

## produces a statistics for overall-DIF
dexter.OverallDIF_ <- function(par1,par2, cov1,cov2)
{
  r=1
  nI=length(par1)
  beta=par1-par2
  Sigma=cov1+cov2
  DIF_test=mahalanobis(beta[-r],rep(0,(nI-1)),Sigma[-r,-r])
  DIF_p=pchisq(DIF_test,(nI-1),lower.tail=FALSE)
  return(list(stat=DIF_test,df=nI-1, p=DIF_p))
}

############################################################
##    functions for reparameterizing to dexter/OPLM       ##
############################################################

# dexter first last 2 index
dexter.first_last2indx = function(first,last) unlist(apply(data.frame(first,last),1,function(x) x[1]:x[2]))

## Makes the reparameterization matrix from log(b) to beta
dexter.makeD <- function(a,first,last)
{
  k = length(a)
  D = matrix(0,k,k)
  tel=1
  for (i in 1:length(first))
  {
    for (j in 1:(last[i]-first[i]+1))
    {
      if (j==1){
        D[tel,tel]=-1/a[tel]
      }else
      {
        D[tel,tel-1]=-1/(a[tel-1]-a[tel])
        D[tel,tel]=1/(a[tel-1]-a[tel])        
      }
      tel=tel+1
    }
  }
  return(D)
}


# geminimaliseerde versie van dexter functie (zonder bayes)
dexter.toOPLM = function(a, b, first, last, H=NULL, fixed_b=NULL)
{
  ## for now remove zero category manually
  if (!is.null(H)) H=H[-first,-first]
  b=b[-first]
  if (!is.null(fixed_b)) fixed_b=fixed_b[-first]
  a=a[-first]
  new_first=first
  new_last=last-1
  for (i in 2:length(first))
  {
    ncat=last[i]-first[i]
    new_first[i]=new_last[i-1]+1
    new_last[i]=new_first[i]+ncat-1
  }
  first=new_first
  last=new_last
  logb=log(b)
  acov=NULL
  ########################
  
  ### CML; b is a single vector
  DD = dexter.makeD(a,first,last)
  if (is.null(fixed_b))
  {
    k=length(b)
    CC=matrix(-1/k,k,k); diag(CC)=(k-1)/k
    AA=CC%*%DD
    ## calculate Deltas and asymp. variance cov matrix
    #  Note: assumes that the first parameters is the reference
    delta=AA%*%logb
    if (!is.null(H))
    {
      acov=solve(H[-1,-1])
      acov=AA[,-1]%*%acov%*%t(AA[,-1])
    }
  }else # if there are fixed parameters we do not normalize
  {
    delta=DD%*%logb
    if (!is.null(H))
    {
      fixed_set=which(!is.na(fixed_b))
      acov=solve(H[-fixed_set,-fixed_set])
      acov=DD[,-fixed_set]%*%acov%*%t(DD[,-fixed_set])
    }
  }
  
  return(list(delta=delta, cov_delta=acov, a=a, first=new_first, last=new_last))
}

## This function expects category thresholds delta, a vector of item_category scores a,
#  and first and last. All without the zero category.
#  It returns dexter parameters b, as well as new a, first and last with the zero category.
dexter.toDexter <- function(delta, a, first, last, re_normalize=TRUE)
{
  if(!is.matrix(delta))
    delta = matrix(delta,nrow=1)
  
  if (re_normalize) 
    delta = delta-apply(delta,1,mean) 

  DDinv = solve(dexter.makeD(a,first,last))
  b = t(exp(apply(delta,1,function(d) DDinv%*%d)))
  
  colnames(b)=colnames(delta)
  
  new_first=first[1]
  new_last=last[1]+1
  if (length(first)>1)
  {
    for (i in 2:length(first))
    {
      nn=last[i]-first[i]
      new_first[i]=new_last[i-1]+1
      new_last=c(new_last,new_first[i]+nn+1)
    }
  }
  new_a=vector("numeric",length(a)+length(first))
  new_b = matrix(1,ncol=length(a)+length(first),nrow=nrow(b))
  new_a[-new_first]=a
  new_b[,-new_first]=b
  
  ## put everything in a (minimal) parms object
  est=list(b=new_b, a=new_a, beta=drop(delta))
  inputs=list(ssIS=list(item_score=new_a),ssI=list(first=new_first,last=new_last))
  parms = list(est=est, inputs=inputs)
  return(parms)
}

# communication with dexter -----------------------------------------------



toDexter = function(b, a, H=NULL, first, last, fixed_b=NULL)
{

  # new first,last
  n_last = last + 1:length(last)
  n_first=n_last
  n_first[1]=first[1]
  n_first[2:length(first)]=n_last[1:(length(first)-1)]+1
  
  # new b , a
  n_b=vector("numeric",n_last[length(n_last)])+1
  n_a=vector("numeric",n_last[length(n_last)])
  n_b[-n_first]=b
  n_a[-n_first]=a
  indx=dexter.first_last2indx(n_first+1,n_last)
  
  if (!is.null(H))
  {
    n_H=matrix(0,length(n_b), length(n_b))
    n_H[indx,indx]=H
    diag(n_H)[n_first]=1
  }else
  {
    n_H=NULL
  }
  
  if (!is.null(fixed_b))
  {
    n_fixed_b = rep(1,n_last[length(n_last)])
    n_fixed_b[-n_first] = fixed_b
  }else
  {
    n_fixed_b=NULL
  }
  OPCML_out = dexter.toOPLM(a=n_a, b=n_b, first=n_first, last=n_last, H=n_H, fixed_b = n_fixed_b)
  return(list(b=n_b, a=n_a, H=n_H, beta=drop(OPCML_out$delta), acov.beta=OPCML_out$cov_delta, first=n_first, last=n_last))
}