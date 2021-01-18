

# This function is for estimating a block diagonal Delta (structured Delta)
# There are b blocks, and each is unstructured made with p_i predictors
# The first block is p1 x p1; the second is p2 x p2, ...
# Sfit and Sres are for Sigma_fit and Sigma_res, and numdir is the number of directions, r for rank of fy.

alg_struct <- function(Sfit, Sres, numdir, r, b = c(p1, p2, p3))
{
	"%^%" =function(M,pow) 
	{	
		if ((dim(M)[1]==1) & (dim(M)[2]==1)) return( as.matrix(M^pow) );
		svdM<-svd(M); return(svdM$u %*% diag(c(svdM$d)^pow) %*% t(svdM$v));
	}

	myexp = function(vDelta0)
	{
 		term <- 0;
		if (numdir < r)
		{
			S0 <- solve(vDelta0);
			newMat <- (S0%^%(0.5))%*%Sfit%*%(S0%^%(0.5));
			eigenMat <- eigen(newMat);
			ubar <- eigenMat$vectors;
			lam <- eigenMat$values;

			for (i in (numdir+1):r)
			{
 				longMat <- (vDelta0%^%(0.5))%*%ubar[,i]%*%t(ubar[,i])%*%(vDelta0%^%(0.5));
				vlongMat <- matrix(longMat, ncol=1);
				term <- term + lam[i]*vlongMat;
			}
		}
		# print("k")
		# outdelta <- solve(t(Gtilde)%*%Gtilde + diag(k))%*%t(Gtilde)%*%(vecSres + term);
		outdelta <- solve(t(Gtilde)%*%Gtilde)%*%t(Gtilde)%*%(vecSres + term);
		return(outdelta)	
	}

	nblocks = length(b)

	p = sum(b); G = list(); k = 0; I = diag(p)

	for (i in 1:b[1])
	{
		for(j in i:b[1])
		{
			k <- k+1;
			G[[k]] <- I[,j]%*%t(I[,i]);
			G[[k]] <- G[[k]] + t(G[[k]])-diag(diag(G[[k]]));
		}
	}

	for (m in 2:nblocks)
	{
		for (i in (sum(b[1:(m-1)])+1):sum(b[1:m]))
		{
			for(j in i:sum(b[1:m]))
			{
				k <- k+1;
				G[[k]] <- I[,j]%*%t(I[,i]);
				G[[k]] <- G[[k]] + t(G[[k]])-diag(diag(G[[k]]));
			}
		}
	}

	Gtilde <- NULL;
	for (m in 1:k) Gtilde <- cbind(Gtilde, c(G[[m]]))
	vecSres <- matrix(Sres, ncol=1)

	# step 1
	#delta0 <- solve(t(Gtilde)%*%Gtilde + diag(k))%*%t(Gtilde)%*%vecSres;
	delta0 <- solve(t(Gtilde)%*%Gtilde)%*%t(Gtilde)%*%vecSres;


	#step 2
	vecDelta0 <- Gtilde%*%delta0;
	Delta0 <- matrix(vecDelta0, ncol=p)
	
	#step 3
	repeat
	{
		delta0n <- myexp(Delta0)#  vDelta <- Delta0
		if (sum(abs(Re(delta0n-delta0))) < 0.01) break
		vecDelta0n <- Gtilde%*%delta0n;
		Delta0 <- matrix(vecDelta0n, ncol=p)
		delta0 <- delta0n;
	}
	return(list(ndeltas=k, m=k, Delta=Re(Delta0))) # To replace m by ndeltas
}

pfc_struct = function(X, y, fy, numdir=2, structure="unstr", b=c(p1, p2, p3))
{
	"%^%"<-function(M, p) 
	{ 
		# This operator compute the M^p where M is a matrix;
		if (!is.matrix(M)) stop("Argument of power is not a matrix in pfc")
		if (prod(dim(M)==list(1,1))) return( as.matrix(M^p) );
		svdM <- svd(M); return(svdM$u%*%diag(c(svdM$d)^p)%*%t(svdM$v));
	}
	Trace<-function(X)
	{	
		if (!is.matrix(X)) stop("Argument to Trace is not a matrix in pfc");
		return(sum(diag(X))) 
	}
	vnames <- dimnames(X)[[2]];
	if (is.null(vnames)) vnames <- paste("X", 1:p, sep="");

	r <- dim(fy)[2]; X <- as.matrix(X); 
	op <- dim(X); n <- op[1]; p <- op[2]; 	
	
	Muhat <- apply(X, 2, mean);
	Xc <-  scale(X, TRUE, FALSE);

	P_F <- fy%*%solve(t(fy)%*%fy)%*%t(fy);
	Sigmahat <- cov(X)		#t(Xc)%*%Xc/n;
	Sfit <- t(Xc)%*%P_F%*%Xc/n; 
	Sres <- Sigmahat - Sfit;

	struct_list <- alg_struct(Sfit, Sres, numdir=numdir, r=r, b = b)

	Deltahat <- struct_list$Delta
	m <- struct_list$m # Check later if this is needed

	# NOTE: This assumes Deltahat is given by the Appendix B algorithm 
	# in "PFC for Dimension Reduction in Regression" and the max likelihood
	# is given by equation 4 (This is the assumption in the LRT on page 496 second paragraph)
	sqrt_Dhat <- Deltahat%^%0.5; 
	Inv_Sqrt_Dhat <- solve(sqrt_Dhat);
	lf_matrix <- Inv_Sqrt_Dhat%*%Sfit%*%Inv_Sqrt_Dhat;
	all_evalues <- eigen(lf_matrix, symmetric=T)$values;
	evalues <- all_evalues[1:numdir];

	Vhat <- eigen(lf_matrix, symmetric=T)$vectors;
	Vhati <- matrix(Vhat[,1:numdir], ncol=numdir);
	Gammahat <- (Deltahat%^%0.5)%*%Vhati%*%solve((t(Vhati)%*%Deltahat%*%Vhati)%^%0.5);  
	dimnames(Gammahat)<- list(vnames, paste("Dir", 1:numdir, sep=""));

	Khat<-diag(0, p); 
	if (numdir < min(ncol(fy),p)) {diag(Khat)[(numdir+1):min(ncol(fy), p )]<- all_evalues[(numdir+1):min(ncol(fy), p)]};
	dimnames(Deltahat) <- list(vnames, vnames);
	InvDeltahat <- solve(Deltahat)
	Betahat <- ((t(Vhati)%*%Deltahat%*%Vhati)%^%0.5)%*%t(Vhati)%*%solve(Deltahat%^%0.5)%*%t(Xc)%*%fy%*% solve(t(fy)%*%fy);

	temp0 <- -(n*p/2)*log(2*pi);
	temp1 <- -(n/2)*log(det(Deltahat)); 
	temp2 <- -(n/2)*Trace(InvDeltahat%*%Sres);
	temp3 <- 0 

#	if (numdir < min(ncol(fy),p)) 
#	temp3 <- -(n/2)*sum(all_evalues[(numdir+1):p]);

#	if (numdir < p) temp3 <- -(n/2)*sum(all_evalues[(numdir+1):p]);
	if (numdir < p) temp3 <- -(n/2)*sum(eigen(InvDeltahat%*%Sres)$values[(numdir+1):p]);

	loglik <- Re(temp0 + temp1 + temp2 + temp3);

	numpar <- p + (p-numdir)*numdir + numdir*ncol(fy) + sum(unlist(lapply(b, function(x) x*(x+1)/2)));
	aic <- -2*loglik + 2*numpar;
	bic <- -2*loglik + log(n)*numpar;

	return(list(m=m, Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, evalues=evalues, 
				loglik=loglik, aic=aic, bic=bic, numpar=numpar, numdir=numdir, sizes=b));
}



FindGroups <- function(fit0, X, y, fy, p, show_plt=FALSE)
{
	# obtain correlation
	sig_cor <- cov2cor(fit0$Deltahat)

	# build dissimilarity matrix
	sig_diss <- apply(sig_cor, 1, function(x) 1-abs(x))

	# perform agglomerative nesting with complete linkage (hierarchical clustering)
	library(cluster)
	hclust <- agnes(sig_diss, diss=TRUE, method="complete")
  
	xx <- as.dendrogram(hclust)

	# re-construct the data frame to group dimensions
	Xclust <- X[,hclust$order]

	# go through each cluster level and compare structured and unstructured
	clusters <- list()
	clusters[[1]] <- xx

	# go through all merge levels
	levels <- c()
	labels <- c()
	best_height <- 0  # height in clustering tree of best fit
  fits = list()
  
	for (c in 1:(p-1)){

	    # keep track of the sub-cluster with max disimilarity
	    max_idx <- 0
	    max_dis <- 0

	    for (i in 1:length(clusters)){
	    	tmp_dis <- attr(clusters[[i]],"height")

		    if(tmp_dis > max_dis){
		      max_dis = tmp_dis
		      max_idx = i
		    }
	    }

    # reform the clusters 
    tclusters <- list()
    if(length(clusters) == 1){
    	tclusters[[1]] <- clusters[[max_idx]][[1]]
	    tclusters[[2]] <- clusters[[max_idx]][[2]]
    }else if (max_idx == 1){
    	tclusters[[1]] <- clusters[[max_idx]][[1]]
      tclusters[[2]] <- clusters[[max_idx]][[2]]
      
	    for(i in 2:length(clusters)){
 	      tclusters[[i+1]] <- clusters[[i]]
	    }
    } else if(max_idx == length(clusters)){
	      for(i in 1:(length(clusters)-1)){
 	        tclusters[[i]] <- clusters[[i]]
	      }
    	tclusters[[length(clusters)]] <- clusters[[max_idx]][[1]]
      tclusters[[length(clusters)+1]] <- clusters[[max_idx]][[2]]
    }else{
      for(i in 1:(max_idx-1)){
 	        tclusters[[i]] <- clusters[[i]]
	      }
      tclusters[[max_idx]] <- clusters[[max_idx]][[1]]
      tclusters[[max_idx+1]] <- clusters[[max_idx]][[2]]
	    for(i in (max_idx+1):length(clusters)){
        tclusters[[i+1]] <- clusters[[i]]
      }
    }   
    clusters = tclusters

    # finally determine the cluster sizes to pass to the structured 
    # PFC routine
    sizes <- c()
    for(i in 1:length(clusters)){
	    sizes <- c(sizes, attr(clusters[i][[1]], "members"))
    }

    # save off the sizes for display
    #print(sizes)
    #labels <- c(labels,paste(sizes,sep='',collapse=''))

    # calculate structured fit 
    fits[[c]] = pfc_struct(Xclust, y, fy, numdir=fit0$numdir, structure="unstr", b=sizes)
  }
  return(list(fits=fits, clusters=hclust))
}

struct_test <- function(fit0, fitall)
{
  holder = data.frame(matrix(vector(), ncol=6, nrow=length(fitall)+1), stringsAsFactors=F)
  colnames(holder) = c("ngroups", "aics", "bics", "loglik", "numpar", "pvalues")

  holder[1,1] = 1
  holder[1,2] = fit0$aic
  holder[1,3] = fit0$bic
  holder[1,4] = fit0$loglik
  holder[1,5] = fit0$numpar
  holder[1,6] = NA
  
  for(i in 1:length(fitall)){
    fit = fitall[[i]]
    statistic <- 2*(fit0$loglik - fit$loglik);
    thedf <- fit0$numpar - fit$numpar; 
    pvalue <- 1 - pchisq(statistic, df=thedf);

    holder[i+1,1] = length(fit$sizes)
    holder[i+1,2] = fit$aic
    holder[i+1,3] = fit$bic
    holder[i+1,4] = fit$loglik
    holder[i+1,5] = fit$numpar
    holder[i+1,6] = pvalue
  }
  return(holder)
}

find_best_fit <- function(fit0, fits, clusters, p, show_plt=FALSE)
{

	use_unstruct <- FALSE
	aics <- c(fit0$aic)
	lhoods <- c(fit0$loglik)
	labels <- c('unstructured')
	bidx <- 1

	# start by comparing the unstructured fit to the highest cluster level
	fit <- fits[[1]]
    statistic <- 2*(fit0$loglik - fit$loglik);
    thedf <- fit0$numpar - fit$numpar; 
    pvalue <- 1 - pchisq(statistic, df=thedf);
	bfit <- list()

	print(statistic)
	print(pvalue)

	if(pvalue < .05 && statistic > 0)
	{
		use_unstruct <- TRUE
		bfit <- fit0
	}

	# now compare the highest cluster level to subsequent levels
	if(use_unstruct == FALSE){

	for(c in 1:length(fits))
	{
		
		# note this is valid because adjacent cluster levels are nested		
	    statistic <- 2*(fits[[1]]$loglik - fits[[c]]$loglik);
	    thedf <- fits[[1]]$numpar - fits[[c]]$numpar; 
	    pvalue <- 1 - pchisq(statistic, df=thedf);

		if(pvalue < .05){
			break
		}
		

		bfit <- fits[[c]]
		bidx <- bidx + 1

	}
	}

	for(c in 1:length(fits))
	{
		aics <- c(aics,fits[[c]]$aic)
		lhoods <- c(lhoods,fits[[c]]$loglik)
		labels <- c(labels,paste(fits[[c]]$sizes,sep=' ', collapse=' '))

	}

    cols <- rep('magenta', p)
    cols[bidx] <- 'blue'

	# plot clusters, aics, and likelihoods
    if(show_plt){
	par(mfrow=c(2,2)); 
	plot(clusters, nmax= 5)
    barplot(aics,main="Grouping vs. AIC", names=labels,col=cols,las=2)
    barplot(lhoods,main="Grouping vs. Loglikelihood", names=labels,col=cols,las=2)
    }

	return(bfit)

}

orthonorm <- function(u) 
{
  if (is.null(u)) return(NULL)
  if (!(is.matrix(u))) u <- as.matrix(u);
  dd <- dim(u); n <- dd[1]; p <-dd[2];
  
  if (prod(abs(La.svd(u)$d) > 1e-08) == 0) stop("collinears vectors in orthonorm")
  if (n < p)
  {
    warning("There are too much vectors to orthonormalize in orthonorm.")
    u <- as.matrix(u[, 1:p])
    n <- p
  }
  v <- u;
  if (p > 1)
  {
    for (i in 2:p)
    {
      coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[, 1:(i - 1)]));
      v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = n) %*% matrix(coef.proj, nrow = i - 1)
    }
  }
  coef.proj <- 1/sqrt(diag(crossprod(v)))
  return(t(t(v) * coef.proj))
}

angle <- function(U,V)
{
  # Computes the angle between two one-dimensional subspaces
  U = orthonorm(as.matrix(U)); 
  V = orthonorm(as.matrix(V)); 
  
  out <- t(U)%*%V;  
  if (any(out > 1)) {index1 = which(out > 1); out[index1]=1}
  if (any(out < -1)) {index2 = which(out < -1); out[index2]=-1}
  return(180*acos(out)/pi)
}


Angle <- function(U,V)
{
  # Computes the angle between two subspaces
  if(!is.matrix(U)) U = orthonorm(as.matrix(U)); 
  if(!is.matrix(V)) V = orthonorm(as.matrix(V)); 
  if(ncol(U) > ncol(V)) {temp = U; U = V; V = temp}
  
  M = matrix(ncol=ncol(V), nrow=ncol(U))
  for (i in 1:ncol(U)) { for (j in 1:ncol(V)) M[i,j] = t(U[,i])%*%V[,j] }
  if ((ncol(V) == 1) & (ncol(U) == 1)) out = drop(M) else out <- det(M%*%t(M)) 
  return(180*acos(out)/pi)
}


PAngle <- function(U,V)
{
  # Principal Angles between two subspaces ---
  angles <- function(u,W)
  { 
    out <- t(u)%*%W; 
    if (any(out > 1)) {index1 = which(out > 1); out[index1]=1}
    if (any(out < -1)) {index2 = which(out < -1); out[index2]=-1}
    return(180*acos(out)/pi)
  }
  if(!is.matrix(U)) U = orthonorm(as.matrix(U));
  if(!is.matrix(V)) V = orthonorm(as.matrix(V));
  if(ncol(U) > ncol(V)) {temp = U; U = V; V = temp}
  
  tempV = V
  thetas = vector()
  for (i in 1:ncol(U))
  {
    out_thetas = angles(U[,i], tempV);
    thetas[i] = min(out_thetas)[1]
    if (i != ncol(U)) tempV = tempV[,-which(out_thetas == thetas[i])]
  }
  return(sort(thetas, decreasing = TRUE))
}


Trace<-function(X)
{
  # This function returns the trace of a matrix
  if (!is.matrix(X)) stop("Not correct argument")
  return(sum(diag(X)))
}


"%^%"<-function(M,pow) 
{
  # This operator compute the M^pow where M isa matrix;
  if (!is.matrix(M)) stop("The argument is not a matrix")
  if ((dim(M)[1]==1) & (dim(M)[2]==1)) return( as.matrix(M^pow) );
  
  svdM = svd(M); 
  return(svdM$u%*%diag(c(svdM$d)^pow)%*%t(svdM$v));
}

  
