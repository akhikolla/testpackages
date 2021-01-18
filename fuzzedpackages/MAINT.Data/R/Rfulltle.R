RunfooforallComb <- function(n,k,foo,ResType=c("none","list","vector","lstandvct","onlybest","matrix"),
  matncol=NULL,Retbstcrt=TRUE,RetbstSet=TRUE,...)
{
  ResType <- match.arg(ResType)
  nCk <- choose(n,k)
  if (ResType=="list")  {
    Res <- vector("list",nCk)
  }  else if (ResType=="vector")  {
    Res <- numeric(nCk)
  }  else if (ResType=="lstandvct")  { 
    Res <- list(PRes=vector("list",nCk),CmpCrt=numeric(nCk))
  }  else if (ResType=="matrix")  {
    if (is.null(matncol))  { stop("Argument matncol is missing\n") } 
    Res <- matrix(nrow=nCk,ncol=matncol)
  }  
  bstcrt <- -Inf

  for (i in 1:nCk)  {
    if (i==1)  {
      v <- 1:k
    }  else  {
      for (j in k:1)  {
        if (v[j] < n-(k-j))
        {
          v[j] <- v[j]+1
          if (j<k)  { for (j1 in (j+1):k) v[j1] <- v[j1-1]+1 }
          break
        }
      }
    }
    if (ResType=="none")
    {
      if (Retbstcrt || RetbstSet)  {
        fooRes <- foo(v,...)
      }  else   {
        foo(v,...)
      }
    }  else if (ResType=="list")  {
      Res[[i]] <- foo(v,...)
    }  else if (ResType=="vector")  {
      Res[i] <- foo(v,...)
    }  else if (ResType=="lstandvct") {
      fooRes <- foo(v,...)
      Res$PRes[[i]] <- fooRes$PRes
      Res$CmpCrt[i] <- fooRes$CmpCrt
    }  else if (ResType=="onlybest")  {
      fooRes <- foo(v,...)
      if (fooRes$CmpCrt > bstcrt)
      {
        Res <- fooRes$PRes
        bstcrt <- fooRes$CmpCrt
      }
    } else if (ResType=="matrix")  {
      Res[i,] <- foo(v,...)
      if (fooRes$CmpCrt > bstcrt)
      {
        if (Retbstcrt)  { bstcrt <- fooRes$CmpCrt }
        if (RetbstSet)  { bestv <- v }
      }
    }
  }

  if (ResType=="none")
  {
    if (Retbstcrt && RetbstSet)  {
      return(list(bstcrt=bstcrt,bestSet=bestv))
    }  else if (Retbstcrt && !RetbstSet)   {
      return(bstcrt)
    }  else if (!Retbstcrt && RetbstSet)  {
      return(bestv)
    }
  }  else {
    if (Retbstcrt && RetbstSet)  {
      return(list(bstcrt=bstcrt,bestSet=bestv,Res=Res))
    }  else if (Retbstcrt && !RetbstSet)  {
      return(list(bstcrt=bstcrt,Res=Res))
    }  else if (!Retbstcrt && RetbstSet)  {
      return(list(bestSet=bestv,Res=Res))
    }  else if (!Retbstcrt && !RetbstSet)  {
      return(Res)
    }
  }
}

Subsmplmle <- function(SmplInd,FullData,Config,SelCrit,...)
{
  Res <- IdtNmle(FullData[SmplInd,],CovCaseArg=FALSE,Config=Config,SelCrit=SelCrit,...) 
  list(CmpCrt=Res@logLiks[Res@BestModel])
}

Rfulltle <- function(Idt,k=ceiling((Idt@NObs+2*Idt@NIVar+1)/2),Config=2,SelCrit=c("BIC","AIC"),force=FALSE,...)
{
  if (!force)
  {
    maxnCk <- 10000
    nCk <- choose(Idt@NObs,k) 
    if (nCk> maxnCk)  {
      stop( paste("fulltle might take too long since",nCk,"different subsets need to be evaluated.\n",
        "To proceed anyway set the 'force' argument to TRUE, otherwise try the fasttle method instead.\n") )
    }
  }					
  bestsol <- RunfooforallComb(Idt@NObs,k,Subsmplmle,ResType="none",FullData=Idt,Config=Config,SelCrit=match.arg(SelCrit),...)
  list(LogLik=bestsol$bstcrt,Set=bestsol$bestSet)
}

