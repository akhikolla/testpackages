#' @title  Principal Component Analysis Based on Least Trimmed Squaers (LTS-PCA)
#' @description the function that computes LTS-PCA
#' @param x the input data matrix
#' @param q the dimension of the PC subspace
#' @param alpha the robust parameter which takes value between 0 to 0.5, default is 0.5
#' @param b.choice intial loading matrix; by default is NULL and the deterministic starting values will be computed by the algorithm
#' @param tol convergence criterion
#' @param N1  the number controls the updates for a without updating b in the concentration step
#' @param N2 the number controls outer loop in the concentration step
#' @param N2bis the number controls the outer loop for the selected b
#' @param Npc the number controls the inner loop
#' @return the object of class "ltspca" is returned \cr
#' \item{b}{the unnormalized loading matrix}
#' \item{mu}{the center estimate}
#' \item{ws}{if the observation in included in the h-subset \code{ws}=1; otherwise \code{ws}=0}
#' \item{best.cand}{the method which computes the best deterministic starting value in the concentration step}
#' @examples
#' \dontrun{
#' ltspcaM <- ltspca(x = x, q = 2, alpha = 0.5)
#' }
#' @author Cevallos Valdiviezo
#' @references Cevallos Valdiviezo, H., Van Aelst, S. (2019), `` Fast computation of robust subspace estimators'', \emph{Computational Statistics & Data Analysis}, 134, 171--185.
ltspca <- function(x, q, alpha=0.5, b.choice=NULL,tol=1e-6, N1=3, N2=2, N2bis=10, Npc=10){


  # Give an initial (clean) h-subset
  determin.start <- function(x, strategy, alpha){

    n <- nrow(x)
    h <- n - floor(n*alpha)

    Qn.x <- apply(x,2,robustbase::Qn)
    rm.x <- which(Qn.x != 0)
    R.x <- apply(x[,rm.x], 2, rank)
    Med.x <- apply(x[,rm.x],2,median)


    if(strategy=="spearman"){
      Y <- R.x
      Qn.y <- apply(Y,2,robustbase::Qn)
      rm.y <- which(Qn.y !=0)
      U <- sweep(Y[,rm.y], 2, (Qn.x[rm.x][rm.y]/Qn.y[rm.y]) , "*")
      rm(Y)
      median.coorwise.U <- apply(U,2,median)
      Qn.coorwise.U <- apply(U,2,robustbase::Qn)
      Z <- scale( U ,center= median.coorwise.U, scale= Qn.coorwise.U)
      obswise.sqnorm <- apply(Z^2, 1, sum)
      rm(Z)
      U <- x[which(obswise.sqnorm <= sort(obswise.sqnorm)[h]),]
    }else if(strategy=="tanh"){
      Z <- scale(x[,rm.x] , center= Med.x, scale= Qn.x[rm.x])
      Y <- apply(Z, 2, tanh)
      Qn.y <- apply(Y,2,robustbase::Qn)
      rm.y <- which(Qn.y !=0)
      U <- sweep(Y[,rm.y], 2, (Qn.x[rm.x][rm.y]/Qn.y[rm.y]), "*")
      rm(Y)
      median.coorwise.U <- apply(U,2,median)
      Qn.coorwise.U <- apply(U,2,robustbase::Qn)
      Z <- scale( U , center= median.coorwise.U, scale= Qn.coorwise.U)
      obswise.sqnorm <- apply(Z^2, 1, sum)
      rm(Z)
      U <- x[which(obswise.sqnorm <= sort(obswise.sqnorm)[h]),]
    }else if(strategy=="normalscore"){
      Y <- qnorm( (R.x - (1/3) ) / (n + (1/3) ) )
      Qn.y <- apply(Y,2,robustbase::Qn)
      rm.y <- which(Qn.y !=0)
      U <- sweep(Y[,rm.y], 2, (Qn.x[rm.x][rm.y]/Qn.y[rm.y]) , "*")
      rm(Y)
      median.coorwise.U <- apply(U,2,median)
      Qn.coorwise.U <- apply(U,2,robustbase::Qn)
      Z <- scale( U , center= median.coorwise.U, scale= Qn.coorwise.U)
      obswise.sqnorm <- apply(Z^2, 1, sum)
      rm(Z)
      U <- x[which(obswise.sqnorm <= sort(obswise.sqnorm)[h]),]
    }else if(strategy=="closesubset"){
      Z <- scale(x[,rm.x] , center= Med.x, scale= Qn.x[rm.x])
      obswise.sqnorm <- apply(Z^2, 1, sum)
      rm(Z)
      U <- x[which(obswise.sqnorm <= sort(obswise.sqnorm)[h]),]
    }else if(strategy=="spatialsign"){
      Z <- scale(x[,rm.x] , center= Med.x, scale=FALSE)
      obswise.norm <- sqrt( apply(Z^2, 1, sum) )
      U <- sweep(Z,1,obswise.norm,"/")
      median.coorwise.U <- apply(U,2,median)
      Qn.coorwise.U <- apply(U,2,robustbase::Qn)
      Z <- scale( U , center= median.coorwise.U, scale= Qn.coorwise.U)
      obswise.sqnorm <- apply(Z^2, 1, sum)
      rm(Z)
      U <- x[which(obswise.sqnorm <= sort(obswise.sqnorm)[h]),]
    }
    return(U)
  }

x<-as.matrix(x)
n<-dim(x)[1]
p<-dim(x)[2]
h<-n-floor(n*alpha)
keep_s <- Inf

cand <- c("spearman", "tanh", "normalscore", "closesubset", "spatialsign")

for(i in cand){

  # Deterministic starting value for B
  U <- determin.start(x=x, strategy=i, alpha=alpha)
  mu <- apply(U,2,mean)
  if(is.null(b.choice))
  {
    set.seed(321)
    bM <- mvtnorm::rmvnorm(p,mean = rep(0,q),sigma = diag(q))
    b <- qr.Q(qr(bM))
    update=findpcs0(x=U, B=b, mu=mu, Npc=Npc, tol=tol)
    b=update$Bmat
    mu=update$mu
    # b <- qr.Q(qr(b))
  }else{
    b <- b.choice
  }

  rm(U)
# Loop starting from this random B
## if (i==2){ print(mu);print(b)}
  update=findpcs2(x=x, B=b, mu=mu, h=h, s=0, N1=N1, N2=N2, Npc=Npc, tol=tol,
                 fixed=!is.null(b.choice))
  b=update$Bmat
  mu=update$mu
  s=update$scale

##  print(paste('scale',s))
##  print(b)

  # Storage of the Nkeep a and B with the lowest values of s
if(s < keep_s){
  opt_cand <- i
  keep_s <- s
  keep_mu <- mu
  keep_b <- b
}

} # End of loop over i in 1:Ncand

#print(paste('scales',keep_s))
# print(keep_b)

mu <- keep_mu
b <- keep_b
s <- keep_s

# Improvement of the best starting value

update=findpcs(x=x, B=b, mu=mu, h=h, s=s, N1=0, N2=N2bis, Npc=Npc, tol=tol,
              fixed=!is.null(b.choice))
      b=update$Bmat
      mu=update$mu
      s=update$scale
      ws <- rep(0,n)
      ws[update$index] <- 1

# Results
res<-list(b=b, mu=mu, best.obj=s, ws=ws, best.cand=opt_cand)
return(res)
}
