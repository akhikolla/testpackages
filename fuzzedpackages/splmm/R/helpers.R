

ZIdentity <- function(Z)
{
  ZId <- diag(dim(Z)[[1]])
  return(list(ZId))
}

covStartingValues <- function(xGroup,yGroup,zGroup,zIdGroup,b,ntot,N,lower=-10,upper=10)
{
  optimize1 <- function(x,y,b) y-x%*%b
  
  optimize2 <- function(gamma,zId,ZtZ)
  {
    H <- zId + exp(2*gamma)*ZtZ
    return(list(H=H))
  }
  
  optimize3 <- function(res,zId,ZtZ,gamma)
  {
    lambda <- optimize2(gamma,zId,ZtZ)
    logdetH <- determinant(lambda$H)$modulus
    quadH <- quad.form.inv(lambda$H,res)
    return(c(logdetH,quadH))
  }
  
  optimize4 <- function(gamma)
  {
    optH <- mapply(optimize3,resGroup,zIdGroup,ZtZ=ztzGroup,MoreArgs=list(gamma=gamma))
    H1 <- optH[1,]
    H2 <- optH[2,]
    
    fn <- ntot*log(sum(H2)) + sum(H1)
    fn
  }
  
  optimize5 <- function(z) tcrossprod(z)
  
  resGroup <- mapply(optimize1,x=xGroup,y=yGroup,MoreArgs=list(b=b),SIMPLIFY=FALSE)
  ztzGroup <- mapply(optimize5,z=zGroup,SIMPLIFY=FALSE)
  
  optRes <- optimize(f=optimize4,interval=c(lower,upper))
  
  gamma <- optRes$minimum
  
  quadH <- mapply(optimize3,resGroup,zIdGroup,ztzGroup,MoreArgs=list(gamma=gamma))[2,]
  
  sig <- sqrt(1/ntot*sum(quadH))
  tau <- exp(gamma)*sig
  objfct <- 1/2*(optRes$objective + ntot*(1-log(ntot)))
  
  return(list(tau=tau,sigma=sig,opt=objfct))
}

VInv <- function(Z,ZId,D,sigma)
{
  det.check = det(sigma^2*ZId + quad.tform(D,Z))
  
  if(det.check!=0){
    vInverse <- solve(sigma^2*ZId + quad.tform(D,Z),tol = 1e-40)
  }else{
    vInverse <- ginv(sigma^2*ZId + quad.tform(D,Z)) 
  }
  
  
  return(list(vInverse))
}

nlogdet <- function(LGroup)
{
  nlogdetfun <- function(L)
  {
    -1/2*determinant(L)$modulus[1]
  }
  
  sum(mapply(nlogdetfun,LGroup))
}


HessianMatrix <- function(xGroup,LGroup,activeSet,N,hessian,mat)
{
  for (i in 1:N)
  {  
    mat[i,] <- diag(t(xGroup[[i]][,activeSet,drop=FALSE])%*%LGroup[[i]]%*%xGroup[[i]][,activeSet,drop=FALSE])
  }
  hessian[activeSet] <- apply(mat,2,sum)
  return(hessian)
}

as1 <- function(xGroup,LGroup,activeSet,N)
{
  fs <- function(x,l,a) {l%*%x[,a]}
  
  SGroup <- mapply(fs,xGroup,LGroup,MoreArgs=list(a=activeSet),SIMPLIFY=FALSE)
  
  return(SGroup)
}

as2 <- function(x,y,b,j,activeSet,group,sGroup)
{
  r <- y-x[,-c(j),drop=FALSE]%*%b[-j]
  rGroup <- split(r,group)
  
  as3 <- function(s,r,j) crossprod(r,s[,j])
  
  ma <- mapply(as3,sGroup,rGroup,MoreArgs=list(j=match(j,activeSet)))
  
  sumMa <- sum(ma)
  
  return(sumMa)
}

SoftThreshold <- function(z,g)
{
  sign(z)*max(abs(z)-g,0)
}


ArmijoRule_b <- function(xGroup,yGroup,LGroup,b,j,cut,HkOldJ,HkJ,JinNonpen,lambda,nonpen,penalty,ll1,ll2,converged,control)
{
  b.new <- b
  bJ <- b[j]
  grad <- -cut + bJ*HkOldJ
  if(JinNonpen){
    dk <- -grad/HkJ
    } else {
    
    if(penalty=="lasso"){
      dk <- MedianValue(grad,HkJ,lambda/weights[j],bJ)
    }else if(penalty=="scad"){
      dk <- ScadValue(grad,HkJ,lambda/weights[j],bJ)
    }
    
    
    }
  
  
  if (dk!=0)
  { 
    # calculate delta_k
    if (JinNonpen) deltak <- dk*grad + control$gamma*dk^2*HkJ
    else deltak <- dk*grad + control$gamma*dk^2*HkJ + lambda/weights[j]*(abs(bJ+dk)-abs(bJ))
    
    fctOld <- ObjFunction_b(xGroup=xGroup,yGroup=yGroup,LGroup=LGroup,b=b,
                          lambda=lambda,nonpen=nonpen,penalty=penalty,ll1=ll1,ll2=ll2)
    for (l in 0:control$maxArmijo)
    { 
      b.new[j] <- bJ + control$a_init*control$delta^l*dk
      
      fctNew <- ObjFunction_b(xGroup=xGroup,yGroup=yGroup,LGroup=LGroup,b=b.new,
                            lambda=lambda,nonpen=nonpen,penalty=penalty,ll1=ll1,ll2=ll2)
      addDelta <- control$a_init*control$delta^l*control$rho*deltak
      if (fctNew <= fctOld + addDelta)
      {
        b[j] <- bJ + control$a_init*control$delta^l*dk
        fct <- fctNew
        break
      }
      if (l==control$maxArmijo)
      {
        converged <- converged + 2
        fct <- fctOld
      }
    }
  } 
  return(list(b=b,fct=fct,converged=converged))
}


ObjFunction_b <- function(xGroup,yGroup,LGroup,b,lambda,nonpen,penalty,resGroup=NULL,ll1=ll1,ll2=NULL)
{
  ResAs <- function(x,y,b,activeSet) y-x[,activeSet,drop=FALSE]%*%b[activeSet,drop=FALSE]
  
  tResLRes <- function(L,res) 1/2*quad.form(L,res) 
  
  if (missing(resGroup))
  {
    activeSet <- which(b!=0)
    resGroup <- mapply(ResAs,x=xGroup,y=yGroup,MoreArgs=list(b=b,activeSet=activeSet),SIMPLIFY=FALSE)
  }
  
  if (missing(ll2)) ll2 <- nlogdet(LGroup=LGroup)
  
  p.term=0
  if(penalty=="lasso"){
    p.term=lambda*sum(abs(b[-nonpen]))
  }else if(penalty=="scad"){
    p.term = sum(scad_group(b[-nonpen], lambda))
  }
  
  ll <- ll1 + ll2 + sum(mapply(tResLRes,LGroup,resGroup)) + p.term
  
  return(Fc=ll)      
}


scad_group = function(beta, lam1, scada=3.7){
  beta=abs(beta)
  tmp=beta<lam1
  
  Dim = length(beta)
  s=rep(0,Dim)
  s[tmp]=lam1*beta[tmp]
  tmp=(beta>lam1)&(beta<=scada*lam1)
  s[tmp]=-(beta[tmp]*beta[tmp]-2*scada*lam1*beta[tmp]+
             lam1*lam1)/2/(scada-1)
  tmp=beta>=scada*lam1
  s[tmp]=(scada+1)*lam1*lam1/2
  
  return(s)
}

scad = function(bj, lambda, a=3.7){
  temp = abs(bj)
  if(temp<=lambda){
    return(lambda)
  }else if(temp>lambda&temp<=a*lambda){
    return((a*lambda-bj)/(a-1))
  }else{
    return(0)
  }
}

MedianValue <- function(grad,hessian,lambda,bj)
{
  median(c((lambda-grad)/hessian,-bj,(-lambda-grad)/hessian))
}

ScadValue <- function(grad,hessian,lambda,bj, a=3.7)
{
  median(c((lambda-grad)/(hessian*(1-1/a)),-bj,(-lambda-grad)/hessian*(1-1/a)))
}


ResAsSplit <- function(x,y,b,f,activeset)
{
  r <- y-x[,activeset,drop=FALSE]%*%b[activeset,drop=FALSE]
  resGroup <- split(r,f)
  return(resGroup)
}

D.Gradient <- function(xGroup,zGroup,LGroup,yGroup,b,N,activeSet){
  
  mat = matrix(0,nrow = q,ncol = q)
  for (i in 1:N) {
    ztvz = t(zGroup[[i]])%*%LGroup[[i]]%*%zGroup[[i]]
    v = t(zGroup[[i]])%*%LGroup[[i]]%*%(yGroup[[i]]-xGroup[[i]][,activeSet,drop=FALSE]%*%b[activeSet])
    
    mat = mat-ztvz+v%*%t(v)
    
  }
  return(mat)
}

D.HessianMatrix <- function(xGroup,zGroup,LGroup,yGroup,b,N,activeSet)
{
  mat = matrix(0,nrow = q^2,ncol = q^2)
  for (i in 1:N)
  {  
    ztvz = t(zGroup[[i]])%*%LGroup[[i]]%*%zGroup[[i]]
    v = t(zGroup[[i]])%*%LGroup[[i]]%*%(yGroup[[i]]-xGroup[[i]][,activeSet,drop=FALSE]%*%b[activeSet])
    mat = mat+-ztvz%x%ztvz+ztvz%x%(v%*%t(v))+t(ztvz%x%(v%*%t(v)))
  }
  
  return(t(mat))
}



MLsigma <- function(zGroup,zIdGroup,resGroup,q,ll1,ll4,true.sigma,D,trace,CovOpt,VarInt)
{
  
  ZPZtGroup <- lapply(zGroup,MLsigmaZPsiZt,Psi=D)
  
  if (CovOpt=="optimize")
  {
    optRes <- optimize(MLsigmaFct,interval=VarInt,ZPZtGroup=ZPZtGroup,zIdGroup=zIdGroup,resGroup=resGroup)
    sigma <- optRes$minimum
  } else if (CovOpt=="nlminb")
  {
    optRes <- nlminb(true.sigma,MLsigmaFct,ZPZtGroup=ZPZtGroup,zIdGroup=zIdGroup,resGroup=resGroup)
    sigma <- optRes$par
  }  
  
  ll <- ll1 + optRes$objective + ll4
  if (trace>3) print(ll)
  
  return(list(sigma=sigma,fct=ll))  
}


MLsigmaZPsiZt <- function(Z,Psi) ZPsiZt <- quad.tform(Psi,Z)

MLsigmaFct <- function(sigma,ZPZtGroup,zIdGroup,resGroup)
{
  LambdaGroup <- mapply(MLsigmaLambda,zIdGroup,ZPZtGroup,MoreArgs=list(sigma=sigma))
  
  ll2 <- mapply(MLpdSymObj,LambdaGroup,resGroup)
  1/2*sum(ll2)
}

MLsigmaLambda <- function(ZId,ZPZt,sigma)
{
  Lambda <- sigma^2*ZId + ZPZt
  return(list(Lambda))
}

MLpdSymFct <- function(thetak,zGroup,zIdGroup,resGroup,sigma,a,b,LPsi)
{
  LPsi[a,b] <- thetak
  Psi <- tcrossprod(LPsi)
  LambdaGroup <- mapply(MLpdSymLambda,Z=zGroup,ZId=zIdGroup,MoreArgs=list(sigma=sigma,Psi=Psi))
  
  ll2 <- mapply(MLpdSymObj,LambdaGroup,resGroup)
  1/2*sum(ll2)
  
}

MLpdSymLambda <- function(Z,ZId,Psi,sigma)
{
  Lambda <- sigma^2*ZId + quad.tform(Psi,Z)
  return(list(Lambda))
}

MLpdSymObj <- function(Lambda,res)  determinant(Lambda)$modulus + MyQuadFormInv(Lambda,res)

matsplitter <- function(M, q) {
  
  k <- kronecker(matrix(1:q^2, q, byrow = TRUE), matrix(1, q, q))
  lapply(split(M, k), matrix, nr = q)
  
  #splitMatrix <- function(mat, nrow) {
  #  split.data.frame(t(mat), ceiling(1:ncol(mat)/ncol(mat)*nrow))
  #}
  #sapply(splitMatrix(M, c), splitMatrix, r)
}

MLloglik <- function(xGroup,yGroup,LGroup,b,ntot,N,activeSet)
{
  
  l1 <- l2 <- numeric(N) ; ll2b <- 0
  for (i in 1:N)
  {
    l1[i] <- -determinant(LGroup[[i]])$modulus
    l2[i] <- quad.form(LGroup[[i]],yGroup[[i]]-xGroup[[i]]%*%b)
  }
  
  ll <- - 1/2*(sum(l1) + sum(l2) + ntot*log(2*pi) + ll2b)
  return(loglik=ll)    
}

MyQuadFormInv <- function(M, x) crossprod(x, solve(M, x, tol = 1e-60))






armijoRule_L <- function(xGroup,yGroup,zGroup,L,l,k,grad,hessian,b,sigma,zIdGrp,linNonpen,lambda,nonpen,penalty,ll1,converged,control)
{
  L.new <- L
  Llk <- L[l,k]
  
  if (linNonpen) {dk <- -grad/hessian} else {
    if(penalty=="lasso"){
      L2norm <- sqrt(sum(L[l,]^2))
      dk <- (-grad-lambda/L2norm*Llk)/(hessian+lambda/L2norm)
    }else if(penalty=="scad"){
      L2norm <- sqrt(sum(L[l,]^2))
      group_scad = scad_group(L2norm, lambda)
      dk <- (-grad-lambda/L2norm*Llk)/(hessian+lambda/group_scad)
    }
  }
  
  if (dk!=0)
  { 
    # calculate delta_k
    if (linNonpen) deltak <- dk*grad + control$gamma*dk^2*hessian
    else {
      
      L.tmp = L
      L.tmp[l,k] = L[l,k]+dk
      deltak <- dk*grad + control$gamma*dk^2*hessian + lambda*(sqrt(sum(L.tmp[l,]^2))-L2norm)
      }
    
    fctOld <- objFunction_L(xGroup=xGroup,yGroup=yGroup,zGroup=zGroup,zIdGrp=zIdGrp,b=b,L=L,sigma=sigma,lambda=lambda,nonpen=nonpen,penalty=penalty,ll1=ll1)
    
    for (j in 0:control$maxArmijo)
    { 
      L.new[l,k] <- Llk + control$a_init*control$delta^j*dk
      
      fctNew <- objFunction_L(xGroup=xGroup,yGroup=yGroup,zGroup=zGroup,zIdGrp=zIdGrp,b=b,L=L.new,sigma=sigma,lambda=lambda,nonpen=nonpen,penalty=penalty,ll1=ll1)
      
      addDelta <- control$a_init*control$delta^j*control$rho*deltak
      if (fctNew <= fctOld + addDelta)
      {
        L[l,k] <- Llk + control$a_init*control$delta^j*dk
        fct <- fctNew
        break
      }
      if (j==control$maxArmijo)
      {
        converged <- converged + 2
        fct <- fctOld
      }
      
    }
  } else{
    fct = objFunction_L(xGroup=xGroup, yGroup=yGroup, zGroup=zGroup, zIdGrp=zIdGrp,
                               b=b, L=L, nonpen=nonpen, sigma=sigma, lambda=lambda, penalty=penalty, ll1=ll1)
  }
  return(list(L=L,fct=fct,converged=converged))
}



objFunction_L <- function(xGroup,yGroup,zGroup,zIdGrp,b,L,sigma,lambda,nonpen,penalty,ll1=ll1,ll2=NULL)
{
  ResAs <- function(x,y,b,activeSet) y-x[,activeSet,drop=FALSE]%*%b[activeSet,drop=FALSE]
  
  tResLRes <- function(L,res) 1/2*quad.form(L,res) 
  
  LGroup = mapply(VInv,Z=zGroup,ZId=zIdGrp,MoreArgs=list(D=L%*%t(L),sigma=sigma))
  
  activeSet <- which(b!=0)
  resGroup <- mapply(ResAs,x=xGroup,y=yGroup,MoreArgs=list(b=b,activeSet=activeSet),SIMPLIFY=FALSE)
  
  
  ll2 <- nlogdet(LGroup=LGroup)
  
  L2norm = sum(apply(L[-nonpen, , drop = FALSE], 1, function(x) sqrt(sum(x^2))))
  
  pen.term.L = 0
  if(penalty=="lasso"){
    pen.term.L = lambda*L2norm
  }else if(penalty=="scad"){
    pen.term.L = sum(scad_group(L2norm, lambda))
  }
  
  
  ll <- ll1 + ll2 + sum(mapply(tResLRes,LGroup,resGroup)) + pen.term.L

  return(Fc=ll)      
}