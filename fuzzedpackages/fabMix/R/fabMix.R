
update_all_y <- function(x_data, mu, SigmaINV, Lambda, z){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- dim(x_data)[1]
	# 2. gibbs for y_i|...
	n_k <- table(z)
	alive <- as.numeric(names(n_k))
	y_new <- array(data = 0, dim = c(n, q))
	for(k in alive){
		ind <- which(z == k)
		center_x <- x_data[ind, ] - matrix(mu[k,], nrow = as.numeric(n_k[as.character(k)]), ncol = p, byrow=TRUE)
		Alpha <- t(Lambda[k, , ]) %*% SigmaINV %*% Lambda[k, , ]
		diag(Alpha) <- diag(Alpha) + 1
		Alpha <- solve(Alpha)
		tmp <- Alpha %*% t(Lambda[k, , ]) %*% SigmaINV
		y_mean <- t(apply(center_x, 1, function(tk){return(tmp %*% tk)} ))
		if(q == 1){ y_mean <- t(y_mean) }
		y_new[ind, ] <- t( apply( y_mean, 1, function( tk ){ return( mvrnorm(n = 1, mu = tk, Sigma = Alpha) ) } ) )
	}

	results <- vector("list", length=1)
	results[[1]] <- y_new
	names(results) <- c("y")
	return(results)
}

update_all_y_Sj <- function(x_data, mu, SigmaINV, Lambda, z){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- dim(x_data)[1]
	# 2. gibbs for y_i|...
	n_k <- table(z)
	alive <- as.numeric(names(n_k))
	y_new <- array(data = 0, dim = c(n, q))
	for(k in alive){
		ind <- which(z == k)
		center_x <- x_data[ind, ] - matrix(mu[k,], nrow = as.numeric(n_k[as.character(k)]), ncol = p, byrow=TRUE)
		Alpha <- t(Lambda[k, , ]) %*% SigmaINV[k, ,] %*% Lambda[k, , ]
		diag(Alpha) <- diag(Alpha) + 1
		Alpha <- solve(Alpha)
		tmp <- Alpha %*% t(Lambda[k, , ]) %*% SigmaINV[k,,]
		y_mean <- t(apply(center_x, 1, function(tk){return(tmp %*% tk)} ))
		if(q == 1){ y_mean <- t(y_mean) }
		y_new[ind, ] <- t( apply( y_mean, 1, function( tk ){ return( mvrnorm(n = 1, mu = tk, Sigma = Alpha) ) } ) )
	}

	results <- vector("list", length=1)
	results[[1]] <- y_new
	names(results) <- c("y")
	return(results)
}


# compute sufficient statistics given y, z, K (and x_data)
compute_sufficient_statistics <- function(y, z, K, x_data){
	cluster_size <- numeric(K)
	p <- dim(x_data)[2]
	q <- dim(y)[2]
	sx  <- array(data = 0, dim = c(K,p))
	sy  <- array(data = 0, dim = c(K,q))
#	sxx <- array(data = 0, dim = c(K,p,p)) #this is not needed at all.
	sxx <- 0
	syy <- array(data = 0, dim = c(K,q,q))
	sxy <- array(data = 0, dim = c(K,p,q))
	for(k in 1:K){
		index <- which(z == k)
		cluster_size[k] <- length(index)
		if( cluster_size[k] > 0){
			sx[k,]  <- colSums(array(x_data[index,],dim = c(cluster_size[k],p)))
			sy[k,]  <- colSums(array(y[index,],dim = c(cluster_size[k],q)))
			if(cluster_size[k] > 1){
				syy[k,,] <- crossprod(y[index, ])
				sxy[k,,] <- crossprod( x_data[index, ],  y[index, ])
			}else{
				syy[k,,] <- y[index,] %*% t(y[index,])
				sxy[k,,] <- x_data[index,] %*% t(y[index,])
			}
		}

	}
	results <- vector("list", length=6)
	names(results) <- c("cluster_size","sx","sy","sxx","syy","sxy")
	results[[1]] <- cluster_size
	results[[2]] <- sx
	results[[3]] <- sy
	results[[4]] <- sxx
	results[[5]] <- syy
	results[[6]] <- sxy
	return(results)
}



compute_sufficient_statistics_given_mu <- function(y, z, K, x_data, mu){
	cluster_size <- numeric(K)
	p <- dim(x_data)[2]
	q <- dim(y)[2]
	sx  <- array(data = 0, dim = c(K,p))
	sy  <- array(data = 0, dim = c(K,q))
#	sxx <- array(data = 0, dim = c(K,p,p))
	sxx <- 0
	syy <- array(data = 0, dim = c(K,q,q))
	sxy <- array(data = 0, dim = c(K,p,q))
	for(k in 1:K){
		index <- which(z == k)
		cluster_size[k] <- length(index)
		if( cluster_size[k] > 0){
			sx[k,]  <- colSums(array(x_data[index,],dim = c(cluster_size[k],p)))
			sy[k,]  <- colSums(array(y[index,],dim = c(cluster_size[k],q)))

			if(cluster_size[k] > 1){
				xNEW <- matrix( apply(x_data[index, ], 1, function(tk){return(tk - mu[k,])}), nrow = cluster_size[k], ncol = p, byrow = TRUE)
				syy[k,,] <- crossprod(y[index, ])
				sxy[k,,] <- crossprod( xNEW,  y[index, ])
			}else{
				xNEW <- x_data[index, ] - mu[k, ]
				syy[k,,] <- y[index,] %*% t(y[index,])
				sxy[k,,] <- xNEW %*% t(y[index,])
			}
#			for( i in index){
#				syy[k,,] <- syy[k,,] + y[i,] %*% t(y[i,])
#				sxy[k,,] <- sxy[k,,] + (x_data[i,] - mu[k, ]) %*% t(y[i,])	#v3: edw afairw to mu
#			}
		}

	}
	results <- vector("list", length=6)
	names(results) <- c("cluster_size","sx","sy","sxx","syy","sxy")
	results[[1]] <- cluster_size
	results[[2]] <- sx
	results[[3]] <- sy
	results[[4]] <- sxx
	results[[5]] <- syy
	results[[6]] <- sxy
	return(results)
}


# dirichlet function
myDirichlet <- function(alpha){
        k <- length(alpha)
        theta <- rgamma(k, shape = alpha, rate = 1)
        return(theta/sum(theta))
}


#simulate z and mixture weights from standard Gibbs
update_z_b <- function(w, mu, Lambda, y, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	q <- dim(y)[2]
	probs <- array(data = 0, dim =c(n,K))
	x_var <- array(data = 0, dim = c(p,p))
	diag(x_var) <- 1/diag(SigmaINV)
	for(k in 1:K){
		center_x <- x_data - t(apply(y,1,function(tmp){return(mu[k,] + array(Lambda[k,,], dim = c(p,q)) %*% tmp)}))
		probs[,k] <- log(w[k])  + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){if(anyNA(tmp)){tmp <- rep(1,K)};return(sample(K,1,prob = tmp))})
#	apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}

update_z_b_Sj <- function(w, mu, Lambda, y, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	q <- dim(y)[2]
	probs <- array(data = 0, dim =c(n,K))
	for(k in 1:K){
		x_var <- array(data = 0, dim = c(p,p))
		diag(x_var) <- 1/diag(SigmaINV[k,,])
		center_x <- x_data - t(apply(y,1,function(tmp){return(mu[k,] + array(Lambda[k,,], dim = c(p,q)) %*% tmp)}))
		probs[,k] <- log(w[k])  + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){if(anyNA(tmp)){tmp <- rep(1,K)};return(sample(K,1,prob = tmp))})
#	apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}


# using dmvnorm from package mvtnorm
update_z4 <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	probs <- array(data = 0, dim =c(n,K))
	for(k in 1:K){
		center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
		lpdf <- log(w[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
		probs[,k] <- lpdf
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){if(anyNA(tmp)){tmp <- rep(1,K)};return(sample(K,1,prob = tmp))})
#	apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}

update_z4_Sj <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	probs <- array(data = 0, dim =c(n,K))
	for(k in 1:K){
		center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV[k,,])
		lpdf <- log(w[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
		probs[,k] <- lpdf
	}
	probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
	z <- apply(probs,1,function(tmp){if(anyNA(tmp)){tmp <- rep(1,K)};return(sample(K,1,prob = tmp))})
#	apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
	results <- vector("list", length=2)
	names(results) <- c("w","z")
	results[[1]] <- w
	results[[2]] <- z
	return(results)
}



# using the matrix inversion lemma
update_z2 <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
        probs <- array(data = 0, dim =c(n,K))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- t(Lambda[k,,]) %*% SigmaINV %*% Lambda[k,,]
		diag(x_var) <- diag(x_var) + 1
		x_var <- try(solve(x_var), TRUE)
                if(is.numeric(x_var) == TRUE){
			x_var <- Lambda[k,,] %*% x_var %*% t(Lambda[k,,])
			x_var <- SigmaINV %*% x_var %*% SigmaINV
			x_var <- SigmaINV - x_var
                        probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
                }
        }
        probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
        z <- apply(probs,1,function(tmp){if(anyNA(tmp)){tmp <- rep(1,K)};return(sample(K,1,prob = tmp))})
	#apply(probs,1,function(tmp){if(anyNA(tmp)){tmp <- rep(1,K)};return(sample(K,1,prob = tmp))})
        results <- vector("list", length=2)
        names(results) <- c("w","z")
        results[[1]] <- w
        results[[2]] <- z
        return(results)
}


update_z2_Sj <- function(w, mu, Lambda, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
        probs <- array(data = 0, dim =c(n,K))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
		x_var <- t(Lambda[k,,]) %*% SigmaINV[k,,] %*% Lambda[k,,]
		diag(x_var) <- diag(x_var) + 1
		x_var <- try(solve(x_var), TRUE)
                if(is.numeric(x_var) == TRUE){
			x_var <- Lambda[k,,] %*% x_var %*% t(Lambda[k,,])
			x_var <- SigmaINV[k,,] %*% x_var %*% SigmaINV[k,,]
			x_var <- SigmaINV[k,,] - x_var
                        probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
                }
        }
        probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
        z <- apply(probs,1,function(tmp){if(anyNA(tmp)){tmp <- rep(1,K)};return(sample(K,1,prob = tmp))})
#	apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
        results <- vector("list", length=2)
        names(results) <- c("w","z")
        results[[1]] <- w
        results[[2]] <- z
        return(results)
}





complete.log.likelihood <- function(x_data, w, mu, Lambda, SigmaINV, z){
	n <- length(z)
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]

	probs <- numeric(n)
	alive <- as.numeric(names(table(z)))
	for(k in alive){
		index <- which(z == k)
		center_x <- x_data[index,] - matrix(mu[k,], nrow = length(index), ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV)
		#x_var <- solve(x_var)
		x_var <- try(solve(x_var), TRUE)
		if(is.numeric(x_var) == TRUE){
			probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
		}
	}
	return(sum(probs))
}


complete.log.likelihood_Sj <- function(x_data, w, mu, Lambda, SigmaINV, z){
	n <- length(z)
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]

	probs <- numeric(n)
	alive <- as.numeric(names(table(z)))
	for(k in alive){
		index <- which(z == k)
		center_x <- x_data[index,] - matrix(mu[k,], nrow = length(index), ncol = p, byrow=TRUE)
		x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
		diag(x_var) <- diag(x_var) + 1/diag(SigmaINV[k,,])
		#x_var <- solve(x_var)
		x_var <- try(solve(x_var), TRUE)
		if(is.numeric(x_var) == TRUE){
			probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
		}
	}
	return(sum(probs))
}



update_SigmaINV_faster <- function(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- length(z)

        SigmaINV <- array(data = 0, dim = c(p,p))
        s <- numeric(p)
	alive <- as.numeric(names(table(z)))
	for (k in alive){
		ind <- which(z == k)
		n_k <- length(ind)
		tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE) + t(
			apply
			(
				array(y[ind, ], dim = c(n_k,q)), 1, 
					function(tk)
					{ 
						return(array(Lambda[ k, , ], dim = c(p,q)) %*% tk) 
					}
			)
		)
		s <- s + colSums((x_data[ind, ] - tmp)^2)
	}
        diag(SigmaINV) <- rgamma(p,shape = alpha_sigma + n/2, rate = beta_sigma + s/2)
        return(SigmaINV)
}



update_SigmaINV_faster_Sj <- function(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- length(z)

        SigmaINV <- array(data = 0, dim = c(K,p,p))
	for (k in 1:K){
	        s <- numeric(p)
		ind <- which(z == k)
		n_k <- length(ind)
                if(n_k > 0){
                        tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE)  + t(
			apply
			(
				array(y[ind, ], dim = c(n_k,q)), 1, 
					function(tk)
					{ 
						return(array(Lambda[ k, , ], dim = c(p,q)) %*% tk) 
					}
			)
		)
                        s <- colSums((x_data[ind, ] - tmp)^2)
                }
		diag(SigmaINV[k, , ]) <- rgamma(p,shape = alpha_sigma + n_k/2, rate = beta_sigma + s/2)
	}
        return(SigmaINV)
}



#new in version 3
update_SigmaINV_xCC <- function(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- length(z)

        SigmaINV <- array(data = 0, dim = c(p,p))	# this is redundant: SigmaINV = sI_p
        s <- 0
	alive <- as.numeric(names(table(z)))
	for (k in alive){
		ind <- which(z == k)
		n_k <- length(ind)
		tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE) + t(
			apply
			(
				array(y[ind, ], dim = c(n_k,q)), 1, 
					function(tk)
					{ 
						return(array(Lambda[ k, , ], dim = c(p,q)) %*% tk) 
					}
			)
		)
		s <- s + sum((x_data[ind, ] - tmp)^2)
	}
        diag(SigmaINV) <- rep(rgamma(1,shape = alpha_sigma + n*p/2, rate = beta_sigma + s/2), p)
        return(SigmaINV)
}


#new in version 3
update_SigmaINV_xUC <- function(x_data, z, y, Lambda, mu, K, alpha_sigma, beta_sigma){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	n <- length(z)

        SigmaINV <- array(data = 0, dim = c(K,p,p))
	for (k in 1:K){
	        s <- 0
		ind <- which(z == k)
		n_k <- length(ind)
                if(n_k > 0){
                        tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE)  + t(
			apply
			(
				array(y[ind, ], dim = c(n_k,q)), 1, 
					function(tk)
					{ 
						return(array(Lambda[ k, , ], dim = c(p,q)) %*% tk) 
					}
			)
		)
                        s <- sum((x_data[ind, ] - tmp)^2)
                }
		diag(SigmaINV[k, , ]) <- rep(rgamma(1,shape = alpha_sigma + n_k*p/2, rate = beta_sigma + s/2),p)
	}
        return(SigmaINV)
}


#update OmegaINV
update_OmegaINV <- function(Lambda, K, g, h){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	OmegaINV <- array(data = 0, dim = c(q,q))
	betaVector <- numeric(q)
	for(k in 1:K){
		betaVector <- betaVector + colSums(array(Lambda[k,,]^2,dim = c(p,q)))
	}
	diag(OmegaINV) <- rgamma(q,shape = g + K*p/2, rate = h + betaVector/2)
	return(OmegaINV)
}

#update OmegaINV
update_OmegaINV_Cxx <- function(Lambda, K, g, h){
	p <- dim(Lambda)[2]
	q <- dim(Lambda)[3]
	OmegaINV <- array(data = 0, dim = c(q,q))
	betaVector <- colSums(array(Lambda[1,,]^2,dim = c(p,q)))
	diag(OmegaINV) <- rgamma(q,shape = g + p/2, rate = h + betaVector/2)
	return(OmegaINV)
}

# new in version 4.1
readLambdaValues <- function(myFile,K,p,q){
	l <- array(data = NA, dim = c(K,p,q))
	kpq <- K*p*q
	myline <- as.numeric(read.table(myFile, colClasses = rep('numeric', kpq) ))
	if(length(myline) != kpq){stop("dimensions do not match")}
	iter <- 0
	for(k in 1:K){
		for(i in 1:p){
			for(j in 1:q){
				iter <- iter + 1
				l[k,i,j] <- myline[iter]				
			}
		}
	}
	return(l)
}


#-------------------------------------------------------------------------------------------
### functions for q_0
#-------------------------------------------------------------------------------------------

complete.log.likelihood_q0 <- function(x_data, w, mu, SigmaINV, z){
	n <- length(z)
	p <- dim(x_data)[2]
        probs <- numeric(n)
        alive <- as.numeric(names(table(z)))
        for(k in alive){
                index <- which(z == k)
                center_x <- x_data[index,] - matrix(mu[k,], nrow = length(index), ncol = p, byrow=TRUE)
                x_var <- SigmaINV[k,,]
                probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) 
        }
        return(sum(probs))
}


complete.log.likelihood_q0_sameSigma <- function(x_data, w, mu, SigmaINV, z){
	n <- length(z)
	p <- dim(x_data)[2]
        probs <- numeric(n)
        alive <- as.numeric(names(table(z)))
        x_var <- SigmaINV
	myConstant <- log(det(x_var)) 
        for(k in alive){
                index <- which(z == k)
                center_x <- x_data[index,] - matrix(mu[k,], nrow = length(index), ncol = p, byrow=TRUE)
                probs[index] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*myConstant
        }
        return(sum(probs))
}


compute_sufficient_statistics_q0 <- function(z, K, x_data){
	p <- dim(x_data)[2]
        cluster_size <- numeric(K)
        sx  <- array(data = 0, dim = c(K,p))
        sy  <- 0
#        sxx <- array(data = 0, dim = c(K,p,p))
	sxx <- 0
        syy <- 0
        sxy <- 0
        for(k in 1:K){
                index <- which(z == k)
                cluster_size[k] <- length(index)
                if( cluster_size[k] > 0){
                        sx[k,]  <- colSums(array(x_data[index,],dim = c(cluster_size[k],p)))
                #        for( i in index){
                 #               sxx[k,,] <- sxx[k,,] + x_data[i,] %*% t(x_data[i,])
                  #      }
                }

        }
        results <- vector("list", length=6)
        names(results) <- c("cluster_size","sx","sy","sxx","syy","sxy")
        results[[1]] <- cluster_size
        results[[2]] <- sx
        results[[3]] <- sy
        results[[4]] <- sxx
        results[[5]] <- syy
        results[[6]] <- sxy
        return(results)
}




compute_A_B_G_D_and_simulate_mu_Lambda_q0 <- function(SigmaINV, suff_statistics, K, priorConst1, T_INV, v_r){
	p <- dim(SigmaINV[1,,])[2]
        A <- array(data = 0, dim = c(K,p,p))
        B <- mu <- array(data = 0, dim = c(K,p))
        G <- 0
        D <- 0
#                               (2) mu|Lambda, sufficient statistics, ...
        mu <- array(data = 0, dim = c(K,p))
        Lambdas <- 0
        for(k in 1:K){
                diag(A[k,,]) <- 1/( suff_statistics$cluster_size[k]*diag(SigmaINV[k,,]) + diag(T_INV))
                B[k,] <- SigmaINV[k,,] %*% (suff_statistics$sx[k,] ) + priorConst1
                # this is for simulating mu_k 
                mu_mean <- A[k,,] %*% B[k,]
                mu[k,] <- mvrnorm(n = 1, mu = mu_mean, Sigma = A[k,,])  
        }
        results <- vector("list", length=6)
        results[[1]] <- A
        results[[2]] <- B
        results[[3]] <- G
        results[[4]] <- D
        results[[5]] <- Lambdas
        results[[6]] <- mu
        names(results) <- c("A","B","G","D","Lambdas","mu")
        return(results)
}



compute_A_B_G_D_and_simulate_mu_Lambda_q0_sameSigma <- function(SigmaINV, suff_statistics, K, priorConst1, T_INV, v_r){
	p <- dim(SigmaINV)[2]
        A <- array(data = 0, dim = c(K,p,p))
        B <- mu <- array(data = 0, dim = c(K,p))
        G <- 0
        D <- 0
#                               (2) mu|Lambda, sufficient statistics, ...
        mu <- array(data = 0, dim = c(K,p))
        Lambdas <- 0
        for(k in 1:K){
                diag(A[k,,]) <- 1/( suff_statistics$cluster_size[k]*diag(SigmaINV) + diag(T_INV))
                B[k,] <- SigmaINV %*% (suff_statistics$sx[k,] ) + priorConst1
                # this is for simulating mu_k 
                mu_mean <- A[k,,] %*% B[k,]
                mu[k,] <- mvrnorm(n = 1, mu = mu_mean, Sigma = A[k,,])  
        }
        results <- vector("list", length=6)
        results[[1]] <- A
        results[[2]] <- B
        results[[3]] <- G
        results[[4]] <- D
        results[[5]] <- Lambdas
        results[[6]] <- mu
        names(results) <- c("A","B","G","D","Lambdas","mu")
        return(results)
}




update_z_q0 <- function(w, mu, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
        probs <- array(data = 0, dim =c(n,K))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p , byrow = TRUE)
                probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% SigmaINV[k,,] %*% tmp) )}) + 0.5*sum(log(diag(SigmaINV[k,,]))) 
        }
        probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
        z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
        results <- vector("list", length=2)
        names(results) <- c("w","z")
        results[[1]] <- w
        results[[2]] <- z
        return(results)
}


update_z_q0_sameSigma <- function(w, mu, SigmaINV, K, x_data){
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
        probs <- array(data = 0, dim =c(n,K))
	myConstant <- sum(log(diag(SigmaINV)))
        for(k in 1:K){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p , byrow = TRUE)
                probs[,k] <- log(w[k]) -0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% SigmaINV %*% tmp) )}) + 0.5*myConstant 
        }
        probs <- array(t(apply(probs, 1, function(tmp){return(exp(tmp - max(tmp)))} )),dim = c(n,K))
        z <- apply(probs,1,function(tmp){return(sample(K,1,prob = tmp))})
        results <- vector("list", length=2)
        names(results) <- c("w","z")
        results[[1]] <- w
        results[[2]] <- z
        return(results)
}



update_SigmaINV_faster_q0 <- function(z, mu, K, alpha_sigma, beta_sigma, x_data){
	p <- dim(x_data)[2]
        SigmaINV <- array(data = 0, dim = c(K,p,p))
        alive <- as.numeric(names(table(z)))
        for (k in 1:K){
                s <- numeric(p)
                ind <- which(z == k)
                n_k <- length(ind)
                if(n_k > 0){
                        tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE) 
                        s <- colSums((x_data[ind, ] - tmp)^2)
                }
                diag(SigmaINV[k, , ]) <- rgamma(p,shape = alpha_sigma + n_k/2, rate = beta_sigma + s/2)
        }
        return(SigmaINV)
}


update_SigmaINV_faster_q0_sameSigma <- function(z, mu, K, alpha_sigma, beta_sigma, x_data){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
        SigmaINV <- array(data = 0, dim = c(p,p))
        alive <- as.numeric(names(table(z)))
        s <- numeric(p)

        for (k in alive){
                ind <- which(z == k)
                n_k <- length(ind)
		tmp <- matrix(mu[k, ], nrow = n_k, ncol = p, byrow = TRUE)
		s <- s + colSums((x_data[ind, ] - tmp)^2) 
        }
	diag(SigmaINV) <- rgamma(p,shape = alpha_sigma + n/2, rate = beta_sigma + s/2)
        return(SigmaINV)
}



overfitting_q0 <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q = 0, zStart, gibbs_z, lowerTriangular=TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- 0
	if (q > 0){ stop(paste0('q should be equal to ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}
	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	w.values <- numeric(K)
	# initial values
	iter <- 1
	if(start_values == FALSE){
		for(k in 1:K){
			diag(SigmaINV.values[k,,]) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		tmp1 <- read.table("muValues.txt")
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		w.values <- as.numeric(read.table("wValues.txt"))
		z <- as.numeric(read.table("zValues.txt"))
	}
	###############################################
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- vector('list',length = K)
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
	for (iter in 2:m){
#		2
		suf_stat <- compute_sufficient_statistics_q0(z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_q0(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat,
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
#		cat(paste0("HERE"), "\n")

#		3
		f2 <- update_z_q0(w = w.values, mu = array(mu.values,dim = c(K,p)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		5
		SigmaINV.values <- update_SigmaINV_faster_q0(x_data = x_data, z = z,  
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_q0(x_data = x_data, w = w.values, mu = mu.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				for(k in 1:K){
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
				}
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
#	for(k in 1:K){
#		close(regulonExpressionConnection[[k]])
#		close(LambdaConnection[[k]])
#	}
	setwd("../")

}



overfitting_q0_sameSigma <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q = 0, zStart, gibbs_z, lowerTriangular=TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- 0
	if (q > 0){ stop(paste0('q should be equal to ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
        SigmaINV.values <- array(data = 0, dim = c(p,p))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	w.values <- numeric(K)
	# initial values
	iter <- 1
	if(start_values == FALSE){
		diag(SigmaINV.values) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma)
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		tmp1 <- read.table("muValues.txt")
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
		}
		w.values <- as.numeric(read.table("wValues.txt"))
		z <- as.numeric(read.table("zValues.txt"))
	}
	###############################################
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- vector('list',length = K)
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
	for (iter in 2:m){
#		2
		suf_stat <- compute_sufficient_statistics_q0(z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_q0_sameSigma(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat,
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
#		cat(paste0("here2"),"\n")
#		3
		f2 <- update_z_q0_sameSigma(w = w.values, mu = array(mu.values,dim = c(K,p)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		cat(paste0("here3"),"\n")
#		5
		SigmaINV.values <- update_SigmaINV_faster_q0_sameSigma(x_data = x_data, z = z,  
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)
#		cat(paste0("here4"),"\n")
		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_q0_sameSigma(x_data = x_data, w = w.values, mu = mu.values, 
				SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
			}
		}
	}
	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	setwd("../")
}





#-------------------------------------------------------------------------------------------
# end of q0 functions
#-------------------------------------------------------------------------------------------


################################################################################################################
################################################################################################################

overfittingMFA <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular = TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		overfitting_q0_sameSigma(x_data = x_data, originalX = originalX, 
			outputDirectory = outputDirectory, Kmax = Kmax, 
			m = m, thinning = thinning, burn = burn, g = g, 
			h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
			beta_sigma = beta_sigma, start_values = start_values, 
			q = 0, zStart = zStart, gibbs_z = gibbs_z, lowerTriangular = lowerTriangular)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		diag(SigmaINV.values) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
#			Lambda.values[k, , ] <- matrix(as.matrix(read.table(paste0("LambdaValues",k,".txt"))),nrow = p, ncol = q, byrow=TRUE) 
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
				}
				cat('\n', file = LambdaConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}




overfittingMFA_Sj <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular=TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		overfitting_q0(x_data = x_data, originalX = originalX, 
			outputDirectory = outputDirectory, Kmax = Kmax, 
			m = m, thinning = thinning, burn = burn, g = g, 
			h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
			beta_sigma = beta_sigma, start_values = start_values, 
			q = 0, zStart = zStart, gibbs_z = gibbs_z, lowerTriangular = lowerTriangular)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
			diag(SigmaINV.values[k,,]) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_Sj(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y_Sj(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster_Sj(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_Sj(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
				}
				cat('\n', file = LambdaConnection, append = TRUE)
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}



#new in version 3 (kapote na ftiakseis ta lambda, edw den xreiazontai ola ta connections.)
overfittingMFA_CCU <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular = TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		cat(paste('q = 0 is not currently supported for CCU model.'),'\n')
		#overfitting_q0_sameSigma(x_data = x_data, originalX = originalX, 
		#	outputDirectory = outputDirectory, Kmax = Kmax, 
		#	m = m, thinning = thinning, burn = burn, g = g, 
		#	h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
		#	beta_sigma = beta_sigma, start_values = start_values, 
		#	q = 0, zStart = zStart, gibbs_z = gibbs_z)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		diag(SigmaINV.values) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1)
		for(r in 1:p){
			Lambda.values[1,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
		}
 		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- Lambda.values[1,r,1:v_r[r]]
			}
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV_Cxx(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics_given_mu(y = y, z = z, K = K, x_data = x_data, mu = mu.values)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_CCU(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
				}
				cat('\n', file = LambdaConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}



# new in version 3
overfittingMFA_CUU <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular = TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		cat(paste('q = 0 is not currently supported for CUU model.'),'\n')
		#overfitting_q0(x_data = x_data, originalX = originalX, 
		#	outputDirectory = outputDirectory, Kmax = Kmax, 
		#	m = m, thinning = thinning, burn = burn, g = g, 
		#	h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
		#	beta_sigma = beta_sigma, start_values = start_values, 
		#	q = 0, zStart = zStart, gibbs_z = gibbs_z)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(r in 1:p){
			Lambda.values[1,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
		}

		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- Lambda.values[1,r,1:v_r[r]]
			}
			diag(SigmaINV.values[k,,]) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV_Cxx(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics_given_mu(y = y, z = z, K = K, x_data = x_data, mu = mu.values)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_CUU(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y_Sj(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster_Sj(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_Sj(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
				}
				cat('\n', file = LambdaConnection, append = TRUE)
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}



#new in version 3 
overfittingMFA_CCC <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular=TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		cat(paste('q = 0 is not currently supported for CCC model.'),'\n')
		#overfitting_q0_sameSigma(x_data = x_data, originalX = originalX, 
		#	outputDirectory = outputDirectory, Kmax = Kmax, 
		#	m = m, thinning = thinning, burn = burn, g = g, 
		#	h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
		#	beta_sigma = beta_sigma, start_values = start_values, 
		#	q = 0, zStart = zStart, gibbs_z = gibbs_z)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		diag(SigmaINV.values) <- rep(rgamma(n = 1, shape = alpha_sigma, rate = beta_sigma), p) ## parameterization with mean = g/h and var = g/(h^2)
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1)
		for(r in 1:p){
			Lambda.values[1,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
		}
 		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- Lambda.values[1,r,1:v_r[r]]
			}
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV_Cxx(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics_given_mu(y = y, z = z, K = K, x_data = x_data, mu = mu.values)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_CCU(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_xCC(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
				}
				cat('\n', file = LambdaConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}


# new in version 3
overfittingMFA_CUC <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular = TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		cat(paste('q = 0 is not currently supported for CUC model.'),'\n')
		#overfitting_q0(x_data = x_data, originalX = originalX, 
		#	outputDirectory = outputDirectory, Kmax = Kmax, 
		#	m = m, thinning = thinning, burn = burn, g = g, 
		#	h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
		#	beta_sigma = beta_sigma, start_values = start_values, 
		#	q = 0, zStart = zStart, gibbs_z = gibbs_z)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(r in 1:p){
			Lambda.values[1,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
		}

		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- Lambda.values[1,r,1:v_r[r]]
			}
			diag(SigmaINV.values[k,,]) <- rep(rgamma(n = 1, shape = alpha_sigma, rate = beta_sigma),p) ## parameterization with mean = g/h and var = g/(h^2)
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
#		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
#		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV_Cxx(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics_given_mu(y = y, z = z, K = K, x_data = x_data, mu = mu.values)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_CUU(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y_Sj(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_xUC(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_Sj(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
				}
				cat('\n', file = LambdaConnection, append = TRUE)
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}



#new in version 3
overfittingMFA_UCC <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular=TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		#overfitting_q0_sameSigma(x_data = x_data, originalX = originalX, 
		#	outputDirectory = outputDirectory, Kmax = Kmax, 
		#	m = m, thinning = thinning, burn = burn, g = g, 
		#	h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
		#	beta_sigma = beta_sigma, start_values = start_values, 
		#	q = 0, zStart = zStart, gibbs_z = gibbs_z)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		diag(SigmaINV.values) <- rep(rgamma(n = 1, shape = alpha_sigma, rate = beta_sigma),p) ## parameterization with mean = g/h and var = g/(h^2)
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_xCC(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
				}
				cat('\n', file = LambdaConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}



#new in version 3
overfittingMFA_UUC <- function(x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular=TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if(q == 0){
	# redirecting everything to the corresponding function
		#overfitting_q0(x_data = x_data, originalX = originalX, 
		#	outputDirectory = outputDirectory, Kmax = Kmax, 
		#	m = m, thinning = thinning, burn = burn, g = g, 
		#	h = h, alpha_prior = alpha_prior, alpha_sigma = alpha_sigma, 
		#	beta_sigma = beta_sigma, start_values = start_values, 
		#	q = 0, zStart = zStart, gibbs_z = gibbs_z)
		return(doNothing <- 0) # exit
	}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
#	cat(paste0("a = ", alpha_prior[1], ", p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
#	diag(T_INV) <- (apply(x_data,2,max) - apply(x_data,2,min))^2
	diag(T_INV) <- diag(var(x_data))
#	diag(T_INV) <- rep(1,p)
#	diag(T_INV) <- rep(100,p)
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colSums(x_data)/n
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
#	print("here")
	#############################################
	#
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
			diag(SigmaINV.values[k,,]) <- rep(rgamma(n = 1, shape = alpha_sigma, rate = beta_sigma),p) ## parameterization with mean = g/h and var = g/(h^2)
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
#	u_v <- runif(1)
	for (iter in 2:m){
		
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_data)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_Sj(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_data)
		}else{
			f2 <- update_z_b_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), K = K, x_data = x_data)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y_Sj(x_data = x_data, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_xUC(x_data = x_data, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_Sj(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				for(k in 1:K){
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
				}
				cat('\n', file = LambdaConnection, append = TRUE)
				cat('\n', file = sigmainvConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(LambdaConnection)
	setwd("../")

}











overfittingMFA_missing_values <- function(missing_entries, x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular=TRUE){
	# missing_entries: list which contains the row number (1st entry) and column indexes (subsequent entries) for every row containing missing values

	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data, na.rm = TRUE))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data, na.rm = TRUE)
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		diag(SigmaINV.values) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
		# start from colmeans
		myColMeans <- colMeans(x_data, na.rm = TRUE)
		x_complete <- x_data
		xReplacedValues <- lapply(missing_entries, function(y){ 
					x_complete[y[1], y[-1]] <<- myColMeans[y[-1]] 
				}
			)
		x_mean <- x_complete
		
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		diag(SigmaINV.values) <- as.numeric(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
		x_complete <- as.matrix(read.table("x_complete.txt", header=TRUE))
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	xConnection <- file("x_complete.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
	for (iter in 2:m){
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_complete)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda(SigmaINV = SigmaINV.values, 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_complete)
		}else{
			f2 <- update_z_b(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = SigmaINV.values, K = K, x_data = x_complete)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y(x_data = x_complete, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster(x_data = x_complete, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)
#		6 {missing data}
		mySD <- 1/sqrt(diag(SigmaINV.values))
		xReplacedValues <- lapply(missing_entries, function(f){ 
					nMissing <- length(f) - 1
					j <- z[f[1]]
					myMean <- mu.values[j, f[-1]] + Lambda.values[j,f[-1], ] %*% t(t(y[f[1], ]))
					x_complete[f[1], f[-1]] <<- rnorm(n = nMissing, mean = myMean, sd = mySD[f[-1]]) 
				}
			)
		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				cat(diag(SigmaINV.values), file = sigmainvConnection, '\n', append = TRUE)
				write.table(x_complete, file = xConnection, row.names = FALSE, quote=F, append=TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
				}
				cat('\n', file = LambdaConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(xConnection)
	close(LambdaConnection)
	setwd("../")

}



overfittingMFA_Sj_missing_values <- function(missing_entries, x_data, originalX, outputDirectory, Kmax, m, thinning, burn, g, h, alpha_prior, alpha_sigma, beta_sigma, start_values, q, zStart, gibbs_z, lowerTriangular=TRUE){
	if(missing(originalX)){originalX <- x_data}
	if(missing(gibbs_z)){gibbs_z = 0.05}
	if(missing(zStart)){zStart = FALSE}
	if(missing(x_data)){stop('x_data not provided.')}
	if(missing(q)){stop('q not provided.')}
	p <- dim(x_data)[2]
	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if (q > ledermannBound){ stop(paste0('q should not exceed the Ledermann bound: ', ledermannBound)) }
	n <- dim(x_data)[1]
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	if(missing(Kmax)){Kmax <- 20}
	if(missing(m)){m <- 21000}
	if(missing(burn)){burn <- 1000}
	if(missing(thinning)){thinning <- 10}
	if(missing(g)){g <- 2}
	if(missing(h)){h <- 1}
	if(missing(alpha_prior)){alpha_prior <- 1*rep(1/Kmax,Kmax)}
	if(missing(alpha_sigma)){alpha_sigma <- 2}
	if(missing(beta_sigma)){beta_sigma <- 1}
	if(missing(start_values)){start_values <- FALSE}
	if( start_values == F ){
		dir.create(outputDirectory)
	}
	setwd(outputDirectory)
	K <- Kmax
	# prior parameters
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data, na.rm = TRUE))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data, na.rm = TRUE)
	priorConst1 <- T_INV %*% ksi
	sigma_y2 <- 1/1
	OmegaINV.constant <- array(data = 0, dim = c(q,q)); 
	diag(OmegaINV.constant) <- rep(g/h,q)
#	diag(OmegaINV.constant) <- rep(1000,q)
	OmegaINV.constantINITIAL <- OmegaINV.constant
	SigmaINV.values <- array(data = 0, dim = c(K,p,p))
	Lambda.values <- array(data = 0, dim = c(K,p,q))
	mu.values <- array(data = 0, dim = c(K,p))
	z <- numeric(n)
	y <- array(data = 0, dim = c(n,q))
	w.values <- numeric(K)
	#############################################
	# initial values
	iter <- 1
	if(start_values == FALSE){
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		#diag(SigmaINV.values) <- rgamma(n = p, shape = 1000, rate = 1) 
		for(k in 1:K){
			mu.values[k,] <- rnorm(p,mean = ksi, sd = sqrt( 1/diag(T_INV) ))
			for(r in 1:p){
				Lambda.values[k,r,1:v_r[r]] <- mvrnorm(n = 1, mu = rep(0,v_r[r]), Sigma = omega[1:v_r[r],1:v_r[r]]) 
			}
			diag(SigmaINV.values[k,,]) <- rgamma(n = p, shape = alpha_sigma, rate = beta_sigma) ## parameterization with mean = g/h and var = g/(h^2)
		}
		for(i in 1:n){
			y[i,] <- rnorm(q,mean = 0,sd = 1)
		}
		w.values <- myDirichlet(alpha_prior[1:K])
		z <- sample(K,n,replace = TRUE, prob = w.values)
		if( outputDirectory == 'alpha_1'){
			if(is.numeric(zStart)){
				z <- zStart
				cluster_size <- numeric(K)
				for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
				w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
			}
		}
		# start from colmeans
		myColMeans <- colMeans(x_data, na.rm = TRUE)
		x_complete <- x_data
		xReplacedValues <- lapply(missing_entries, function(y){ 
					x_complete[y[1], y[-1]] <<- myColMeans[y[-1]] 
				}
			)
		x_mean <- x_complete
	}else{
#		cat(paste0('reading starting values... '))	
		diag(OmegaINV.constant) <- as.numeric(read.table('omegainvValues.txt', colClasses = rep('numeric',q))[1,])
		omega <- OmegaINV.constant
		diag(omega) <- 1/(diag(omega))
		tmp1 <- read.table("muValues.txt", colClasses = rep('numeric',Kmax*p))
		tmp2 <- as.matrix(read.table("sigmainvValues.txt"))
		for(k in 1:K){
			mu.values[k, ] <- as.matrix(tmp1[ , k + Kmax*((1:p)-1)])
			diag(SigmaINV.values[k, , ]) <- as.matrix(tmp2[,((k-1)*p + 1):(k*p)]) 
		}
		Lambda.values <- readLambdaValues(myFile = "LambdaValues.txt", K = K, p = p, q = q)
		y <- matrix(as.matrix(read.table('yValues.txt', colClasses = rep('numeric',n*q))), nrow = n , ncol = q)
		w.values <- as.numeric(read.table("wValues.txt", colClasses = rep('numeric',Kmax)))
		z <- as.numeric(read.table("zValues.txt", colClasses = rep('numeric', n)))
#		cat(paste0('done.'),'\n')
		x_complete <- as.matrix(read.table("x_complete.txt", header=TRUE))
	}
	###############################################
	yd <- array(data = 0, dim = c(n,q))
	trueVar <- array(data = 0, dim = c(K,p,p))
	trueVar.values <- array(data = 0, dim = c(K,p,p))
	mhAR <- mhAR1 <- 0
	mhDeltaAR <- 0
	mIter <- 0
	# MCMC sampler
	cluster_size <- numeric(K)
	zOld <- z
	kValues <- numeric(m)
	kValues[iter] <- length(table(z))
	zConnection <- file("zValues.txt",open = "w")
	yConnection <- file("yValues.txt",open = "w")
	sigmainvConnection <- file("sigmainvValues.txt",open = "w")
	omegainvConnection <- file("omegainvValues.txt",open = "w")
	muConnection <- file("muValues.txt",open = "w")
	wConnection <- file("wValues.txt",open = "w")
	logLConnection <- file("k.and.logl.Values.txt",open = "w")
	xConnection <- file("x_complete.txt",open = "w")
	LambdaConnection <- file("LambdaValues.txt", open = "w")
	current_matrix <- vector("list", length = 4)
	names(current_matrix) <- c("A","B","G","D")
	kavatza <- 0
	mySD <- array(data = 0, dim = c(K,p))
	for (iter in 2:m){
	
#		1
		OmegaINV.constant <- update_OmegaINV(Lambda = array(Lambda.values,dim = c(K,p,q)), K = K, g = g, h = h)
#		2
		suf_stat <- compute_sufficient_statistics(y = y, z = z, K = K, x_data = x_complete)
		f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_Sj(SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), 
				suff_statistics = suf_stat, OmegaINV = OmegaINV.constant, 
				K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
		mu.values <- f2$mu
		Lambda.values <- f2$Lambdas
#		3
		u_v <- runif(1)
		if(u_v < gibbs_z){
			f2 <- update_z2_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), SigmaINV = SigmaINV.values, K = K, x_data = x_complete)
		}else{
			f2 <- update_z_b_Sj(w = w.values, mu = array(mu.values,dim = c(K,p)), Lambda = array(Lambda.values,dim = c(K,p,q)), y = y, 
						SigmaINV = array(SigmaINV.values,dim = c(K,p,p)), K = K, x_data = x_complete)
		}
		z <- f2$z
		kValues[iter] <- length(table(z))
		cluster_size <- numeric(K)
		for(k in 1:K){ index <- which(z == k);	cluster_size[k] <- length(index)}	
		w.values <- myDirichlet(alpha_prior[1:K] + cluster_size)
#		4
		y <- array(update_all_y_Sj(x_data = x_complete, mu = mu.values, SigmaINV = SigmaINV.values, Lambda = array(Lambda.values,dim = c(K,p,q)), z = z)$y, dim = c(n, q))
#		5
		SigmaINV.values <- update_SigmaINV_faster_Sj(x_data = x_complete, z = z, y = y, Lambda = array(Lambda.values,dim = c(K,p,q)), 
				mu = array(mu.values,dim = c(K,p)), K = K, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)
#		6 {missing data}
		for(k in 1:K){
			mySD[k,] <-  1/sqrt(diag(SigmaINV.values[k, , ]))
		}
		xReplacedValues <- lapply(missing_entries, function(f){ 
					nMissing <- length(f) - 1
					j <- z[f[1]]
					myMean <- mu.values[j, f[-1]] + Lambda.values[j,f[-1], ] %*% t(t(y[f[1], ]))
					x_complete[f[1], f[-1]] <<- rnorm(n = nMissing, mean = myMean, sd = mySD[j, f[-1]]) 
				}
			)

		if(iter %% thinning == 0){
			if(iter > burn){
				logLValues <- c(kValues[iter], complete.log.likelihood_Sj(x_data = x_data, w = w.values, mu = mu.values, Lambda = Lambda.values, SigmaINV = SigmaINV.values, z = z))
				cat(logLValues, file = logLConnection, '\n', append = TRUE)
				cat(z, file = zConnection, '\n', append = TRUE)
				cat(y, file = yConnection, '\n', append = TRUE)
				cat(w.values, file = wConnection, '\n', append = TRUE)
				cat(mu.values, file = muConnection, '\n', append = TRUE)
				cat(diag(OmegaINV.constant), file = omegainvConnection, '\n', append = TRUE)
				write.table(x_complete, file = xConnection, row.names = FALSE, quote=F, append=TRUE)
				for(k in 1:K){
					for(r in 1:p){
						cat(Lambda.values[k, r, ], " ", file = LambdaConnection, append = TRUE)
					}
					cat(diag(SigmaINV.values[k,,]), " ", file = sigmainvConnection, append = TRUE)
				}
				cat('\n', file = sigmainvConnection, append = TRUE)
				cat('\n', file = LambdaConnection, append = TRUE)
			}
		}
	}

	close(zConnection)
	close(yConnection)
	close(wConnection)
	close(muConnection)
	close(sigmainvConnection)
	close(omegainvConnection)
	close(logLConnection)
	close(xConnection)
	close(LambdaConnection)
	setwd("../")

}




log_dirichlet_pdf <- function(alpha, weights){
	normConstant <- sum( lgamma(alpha) ) - lgamma( sum(alpha) )
	pdf <- sum( (alpha - 1)*log(weights) ) - normConstant
	return(pdf)
}


# for UUU and UCU models
fabMix_UxU <- function(sameSigma = TRUE, dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle, gibbs_z = 1, warm_up_overfitting = 100, warm_up = 500, overfittingInitialization=TRUE, progressGraphs = FALSE, gwar = 0.05, lowerTriangular=TRUE){

	missingRowsIndex <- which(is.na(rowSums(rawData)) == TRUE)
	nMissingRows <- length( missingRowsIndex ) 
	if(nMissingRows > 0){
		stop("The data contains missing values. Use the `fabMix_missing_values()` function.")
	}
	if(progressGraphs == TRUE){
		dev.new(width=15, height=5)
	}
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
	}
	nChains <- length(dirPriorAlphas)
	registerDoParallel(cores = nChains)

	#get rid of other packages messages after printing logo
	nothingHere <- 0
	foreach(nothingHere=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		nothingHere + 1
	}
	mypal <- c(brewer.pal(9, "Set1"), "black") # up to 10 colours

	if( nChains > 1 ){
		if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
		if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
	dir.create(outDir)
	setwd(outDir)
	outputDirs <- paste0(getwd(), '/alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
	if( missing(thinning) ){thinning = 1}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat('\n')
#	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	if(sameSigma == TRUE){
		cat(paste0('-    Parameterization: UCU model'),'\n')
	}else{
		cat(paste0('-    Parameterization: UUU model'),'\n')
	}
	cat(paste0("-    Number of factors: q = ", q,"\n"))
#	cat(paste0('-    Using Nchains = ', nChains),'\n')
#	cat(paste0('-    Target posterior distribution corresponds to alpha = ', dirPriorAlphas[1]),'\n')
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = TRUE, scale = TRUE)
#		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		x_data <- rawData
#		cat('-    The sampler uses raw data (NOT GOOD PRACTICE).','\n')
	}
	kValues <- array(data = NA, dim = c(mCycles, nChains))
	mh_acceptance_rate <- 0
	dir.create('tmpDir')
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	kpq <- Kmax*p*q
	#	initialization
	iteration <- 1
	if(q == 0){warm_up_overfitting = 2*warm_up_overfitting}
	if(overfittingInitialization == TRUE){
		cat(paste('-    (1) Initializing from priors that lead to overfitting... '))
		d_per_cluster = 2*p + p*q + q*(q-1)/2
		if(q == 0){d_per_cluster = 10*p}
		initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_Sj(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}
	}else{
		cat(paste('-    (1) Initializing from random starting values (NOT A GOOD PRACTICE)... '))
		d_per_cluster = 2*p + p*q + q*(q-1)/2
		initialAlphas <- dirPriorAlphas
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_Sj(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_Sj(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
		}
	}
	cat(paste(' OK'),'\n')
	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	LambdaConnection_target <- vector('list',length = Kmax)
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	LambdaConnection_target <- file(paste0(getwd(),"/LambdaValues.txt"),open = "w") 
	swapRate <-  file(paste0(getwd(),"/swapRate.txt"),open = "w")
	ar <- 0 
	LambdaHeader <- character(Kmax*p*q)
	iter <- 0
	for(k in 1:Kmax){
		for(i in 1:p){
			for(j in 1:q){
				iter <- iter + 1
				LambdaHeader[iter] <- paste0("k",k,"_i",i,"_j",j)
			}
		}
	}
	cat(LambdaHeader, file = LambdaConnection_target, '\n', append = TRUE)

	cat(paste('-    (3) Running the sampler... '),'\n')
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	for( iteration in 2:mCycles ){
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_Sj(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}
		
		if(nChains > 1){
			chains <- sample(nChains - 1, 1)
			chains <- c(chains, chains + 1)
			weights[1, ] <- as.numeric(read.table( paste0(outputDirs[ chains[1] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))
			weights[2, ] <- as.numeric(read.table( paste0(outputDirs[ chains[2] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))

			mh_denom <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[2, ] )
			mh_nom   <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[2, ] )
			mh_ratio <- mh_nom - mh_denom
			if( log(runif(1)) < mh_ratio ){
				# dir1 to tmp
				file.copy( 
						from      = paste0(outputDirs[ chains[1] ], '/', file_names), 
						to        = 'tmpDir',
						overwrite = TRUE
					)
				# dir2 to dir1
				file.copy( 
						from      = paste0(outputDirs[ chains[2] ], '/', file_names), 
						to        = outputDirs[ chains[1] ],
						overwrite = TRUE
					)
				# tmp to dir2
				file.copy( 
						from      = paste0('tmpDir', '/', file_names), 
						to        = outputDirs[ chains[2] ],
						overwrite = TRUE
					)
				
				mh_acceptance_rate <- mh_acceptance_rate + 1
			}
		}
		for(myChain in 1:nChains){
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}

		z <- as.numeric(read.table('alpha_1/zValues.txt', colClasses = rep('numeric', n)))
		if( (iteration %% 10) == 0 ){
			if(progressGraphs == TRUE){
				par(mfrow = c(1,3))
				matplot(kValues[1:iteration, ], type = "l")
				points(1:iteration, kValues[1:iteration, 1], type = "b", col = 1)
				matplot(t(x_data), type = "l", col = mypal[as.numeric(as.factor(z))])
				matplot(t(originalX), type = "l", col = mypal[as.numeric(as.factor(z))])
			}
			ar <- round(100*mh_acceptance_rate/iteration, 3)
			cat("\r", paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.        '))
			#cat(paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.'), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				w        <- as.numeric(read.table('alpha_1/wValues.txt', colClasses = rep('numeric',Kmax))) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt', colClasses = rep('numeric',Kmax*p)))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
			if(q > 0){
				y        <- as.numeric(read.table('alpha_1/yValues.txt', colClasses = rep('numeric',n*q)))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt', colClasses = rep('numeric',q)))
				Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues.txt'), colClasses=rep('numeric',kpq) ) )
				cat(Lambda, file = LambdaConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
			}
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
				cat(ar	    , file = swapRate, '\n', append = TRUE)
			}
		}
	}
	cat('\n')
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	close(LambdaConnection_target)
	close(swapRate)
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	setwd("../")
	cat('-    DONE.','\n')
}


# new in version 3
# CUU and CCU models (sameSigma = TRUE => CCU, sameSigma = FALSE => CUU)
fabMix_CxU <- function(sameSigma = TRUE, dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle, gibbs_z = 1, warm_up_overfitting = 100, warm_up = 500, overfittingInitialization=TRUE, progressGraphs = FALSE, gwar = 0.05, lowerTriangular=TRUE){

	missingRowsIndex <- which(is.na(rowSums(rawData)) == TRUE)
	nMissingRows <- length( missingRowsIndex ) 
	if(nMissingRows > 0){
		stop("The data contains missing values. Use the `fabMix_missing_values()` function.")
	}
	if(progressGraphs == TRUE){
		dev.new(width=15, height=5)
	}
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
	}
	nChains <- length(dirPriorAlphas)
	registerDoParallel(cores = nChains)

	#get rid of other packages messages after printing logo
	nothingHere <- 0
	foreach(nothingHere=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		nothingHere + 1
	}
	mypal <- c(brewer.pal(9, "Set1"), "black") # up to 10 colours

	if( nChains > 1 ){
		if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
		if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
	dir.create(outDir)
	setwd(outDir)
	outputDirs <- paste0(getwd(), '/alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
	if( missing(thinning) ){thinning = 1}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat('\n')
#	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	if(sameSigma == TRUE){
		cat(paste0('-    Parameterization: CCU model'),'\n')
	}else{
		cat(paste0('-    Parameterization: CUU model'),'\n')
	}
	cat(paste0("-    Number of factors: q = ", q,"\n"))
#	cat(paste0('-    Using Nchains = ', nChains),'\n')
#	cat(paste0('-    Target posterior distribution corresponds to alpha = ', dirPriorAlphas[1]),'\n')
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = TRUE, scale = TRUE)
#		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		x_data <- rawData
#		cat('-    The sampler uses raw data (NOT GOOD PRACTICE).','\n')
	}
	kValues <- array(data = NA, dim = c(mCycles, nChains))
	mh_acceptance_rate <- 0
	dir.create('tmpDir')
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	kpq <- Kmax*p*q
	#	initialization
	iteration <- 1
	if(q == 0){warm_up_overfitting = 2*warm_up_overfitting}
	if(overfittingInitialization == TRUE){
		cat(paste('-    (1) Initializing from priors that lead to overfitting... '))
		d_per_cluster = 2*p + p*q + q*(q-1)/2
		if(q == 0){d_per_cluster = 10*p}
		initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CCU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CUU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}
	}else{
		cat(paste('-    (1) Initializing from random starting values (NOT A GOOD PRACTICE)... '))
		d_per_cluster = 2*p + p*q + q*(q-1)/2
		initialAlphas <- dirPriorAlphas
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CCU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CUU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_CCU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_CUU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular)
		}
	}
	cat(paste(' OK'),'\n')
	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	LambdaConnection_target <- vector('list',length = Kmax)
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	LambdaConnection_target <- file(paste0(getwd(),"/LambdaValues.txt"),open = "w")  
	LambdaHeader <- character(Kmax*p*q)
	swapRate <-  file(paste0(getwd(),"/swapRate.txt"),open = "w") 
	ar <- 0
	iter <- 0
	for(k in 1:Kmax){
		for(i in 1:p){
			for(j in 1:q){
				iter <- iter + 1
				LambdaHeader[iter] <- paste0("k",k,"_i",i,"_j",j)
			}
		}
	}
	cat(LambdaHeader, file = LambdaConnection_target, '\n', append = TRUE)



	cat(paste('-    (3) Running the sampler... '),'\n')
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	for( iteration in 2:mCycles ){
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CCU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CUU(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}
		if(nChains > 1){
			chains <- sample(nChains - 1, 1)
			chains <- c(chains, chains + 1)
			weights[1, ] <- as.numeric(read.table( paste0(outputDirs[ chains[1] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))
			weights[2, ] <- as.numeric(read.table( paste0(outputDirs[ chains[2] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))

			mh_denom <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[2, ] )
			mh_nom   <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[2, ] )
			mh_ratio <- mh_nom - mh_denom
			if( log(runif(1)) < mh_ratio ){
				# dir1 to tmp
				file.copy( 
						from      = paste0(outputDirs[ chains[1] ], '/', file_names), 
						to        = 'tmpDir',
						overwrite = TRUE
					)
				# dir2 to dir1
				file.copy( 
						from      = paste0(outputDirs[ chains[2] ], '/', file_names), 
						to        = outputDirs[ chains[1] ],
						overwrite = TRUE
					)
				# tmp to dir2
				file.copy( 
						from      = paste0('tmpDir', '/', file_names), 
						to        = outputDirs[ chains[2] ],
						overwrite = TRUE
					)
				
				mh_acceptance_rate <- mh_acceptance_rate + 1
			}
		}
		for(myChain in 1:nChains){
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}

		z <- as.numeric(read.table('alpha_1/zValues.txt', colClasses = rep('numeric', n)))
		if( (iteration %% 10) == 0 ){
			if(progressGraphs == TRUE){
				par(mfrow = c(1,3))
				matplot(kValues[1:iteration, ], type = "l")
				points(1:iteration, kValues[1:iteration, 1], type = "b", col = 1)
				matplot(t(x_data), type = "l", col = mypal[as.numeric(as.factor(z))])
				matplot(t(originalX), type = "l", col = mypal[as.numeric(as.factor(z))])
			}
			ar <- round(100*mh_acceptance_rate/iteration, 3)
			cat("\r", paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.        '))
			#cat(paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.'), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				w        <- as.numeric(read.table('alpha_1/wValues.txt', colClasses = rep('numeric',Kmax))) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt', colClasses = rep('numeric',Kmax*p)))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
			if(q > 0){
				y        <- as.numeric(read.table('alpha_1/yValues.txt', colClasses = rep('numeric',n*q)))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt', colClasses = rep('numeric',q)))
				Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues.txt'), colClasses=rep('numeric',kpq) ) )
				cat(Lambda, file = LambdaConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
			}
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
				cat(ar	    , file = swapRate, '\n', append = TRUE)
			}
		}
	}
	cat('\n')
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	close(LambdaConnection_target)
	close(swapRate)
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	setwd("../")
	cat('-    DONE.','\n')
}


# new in version 3
# CUC and CCC models (sameSigma = TRUE => CCC, sameSigma = FALSE => CUC)
fabMix_CxC <- function(sameSigma = TRUE, dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle, gibbs_z = 1, warm_up_overfitting = 100, warm_up = 500, overfittingInitialization=TRUE, progressGraphs = FALSE, gwar = 0.05, cccStart = FALSE, lowerTriangular=TRUE){

	missingRowsIndex <- which(is.na(rowSums(rawData)) == TRUE)
	nMissingRows <- length( missingRowsIndex ) 
	if(nMissingRows > 0){
		stop("The data contains missing values. Use the `fabMix_missing_values()` function.")
	}
	if(progressGraphs == TRUE){
		dev.new(width=15, height=5)
	}
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
	}
	nChains <- length(dirPriorAlphas)
	registerDoParallel(cores = nChains)

	#get rid of other packages messages after printing logo
	nothingHere <- 0
	foreach(nothingHere=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		nothingHere + 1
	}
	mypal <- c(brewer.pal(9, "Set1"), "black") # up to 10 colours

	if( nChains > 1 ){
		if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
		if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}	
	}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
	dir.create(outDir)
	setwd(outDir)
	outputDirs <- paste0(getwd(), '/alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
	if( missing(thinning) ){thinning = 1}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat('\n')
#	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	if(sameSigma == TRUE){
		cat(paste0('-    Parameterization: CCC model'),'\n')
	}else{
		cat(paste0('-    Parameterization: CUC model'),'\n')
	}
	cat(paste0("-    Number of factors: q = ", q,"\n"))
#	cat(paste0('-    Using Nchains = ', nChains),'\n')
#	cat(paste0('-    Target posterior distribution corresponds to alpha = ', dirPriorAlphas[1]),'\n')
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = TRUE, scale = TRUE)
#		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		x_data <- rawData
#		cat('-    The sampler uses raw data (NOT GOOD PRACTICE).','\n')
	}
	kValues <- array(data = NA, dim = c(mCycles, nChains))
	mh_acceptance_rate <- 0
	dir.create('tmpDir')
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	kpq <- Kmax*p*q
	#	initialization
	if(cccStart == FALSE){
		iteration <- 1
		if(q == 0){warm_up_overfitting = 2*warm_up_overfitting}
		if(overfittingInitialization == TRUE){
			cat(paste('-    (1) Initializing from priors that lead to overfitting... '))
			d_per_cluster = 2*p + p*q + q*(q-1)/2
			if(q == 0){d_per_cluster = 10*p}
			initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
			if(sameSigma == TRUE){
				foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
					overfittingMFA_CCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
						Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
						alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, zStart = zStart, lowerTriangular=lowerTriangular)
				}
			}else{
				foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
					overfittingMFA_CUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
						Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
						alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, zStart = zStart, lowerTriangular=lowerTriangular)
				}
			}
		}else{
			cat(paste('-    (1) Initializing from random starting values (NOT A GOOD PRACTICE)... '))
			d_per_cluster = 2*p + p*q + q*(q-1)/2
			initialAlphas <- dirPriorAlphas
			if(sameSigma == TRUE){
				foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
					overfittingMFA_CCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
						Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
						alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, zStart = zStart, lowerTriangular=lowerTriangular)
				}
			}else{
				foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
					overfittingMFA_CUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
						Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
						alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, zStart = zStart, lowerTriangular=lowerTriangular)
				}
			}
		}
		cat(paste(' OK'),'\n')
		cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
			}
		}
	}else{
		# starting from the ccc model
		iteration <- 1
		cat(paste('-    (1) Initializing from the CCC model with priors that lead to overfitting... '))
		d_per_cluster = 2*p + p*q + q*(q-1)/2
		if(q == 0){d_per_cluster = 10*p}
		initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_CCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, zStart = zStart, lowerTriangular=lowerTriangular)
		}
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_CCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
			if(sameSigma == FALSE){
				tmp <- read.table(paste0(outputDirs[myChain],"/sigmainvValues.txt"))
				tmp <- array(as.numeric(rep(tmp, Kmax)), dim = c(1, p*Kmax) )
				write.table(tmp, file = paste0(outputDirs[myChain],"/sigmainvValues.txt"), append = FALSE, col.names=F, row.names=F)
			}

		}
		cat(paste(' OK'),'\n')
		#cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
		#if(sameSigma == TRUE){
		#	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		#		overfittingMFA_CCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
		#			Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
		#			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z)
		#	}
		#}else{
		#	foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		#		overfittingMFA_CUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
		#			Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
		#			alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z)
		#	}
		#}
	}
	cat(paste(' OK'),'\n')
	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	LambdaConnection_target <- file(paste0(getwd(),"/LambdaValues.txt"),open = "w")  
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	LambdaHeader <- character(Kmax*p*q)
	swapRate <-  file(paste0(getwd(),"/swapRate.txt"),open = "w") 
	ar <- 0
	iter <- 0
	for(k in 1:Kmax){
		for(i in 1:p){
			for(j in 1:q){
				iter <- iter + 1
				LambdaHeader[iter] <- paste0("k",k,"_i",i,"_j",j)
			}
		}
	}
	cat(LambdaHeader, file = LambdaConnection_target, '\n', append = TRUE)




	cat(paste('-    (3) Running the sampler... '),'\n')
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	for( iteration in 2:mCycles ){
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_CUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}
		if(nChains > 1){
			chains <- sample(nChains - 1, 1)
			chains <- c(chains, chains + 1)
			weights[1, ] <- as.numeric(read.table( paste0(outputDirs[ chains[1] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))
			weights[2, ] <- as.numeric(read.table( paste0(outputDirs[ chains[2] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))

			mh_denom <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[2, ] )
			mh_nom   <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[2, ] )
			mh_ratio <- mh_nom - mh_denom
			if( log(runif(1)) < mh_ratio ){
				# dir1 to tmp
				file.copy( 
						from      = paste0(outputDirs[ chains[1] ], '/', file_names), 
						to        = 'tmpDir',
						overwrite = TRUE
					)
				# dir2 to dir1
				file.copy( 
						from      = paste0(outputDirs[ chains[2] ], '/', file_names), 
						to        = outputDirs[ chains[1] ],
						overwrite = TRUE
					)
				# tmp to dir2
				file.copy( 
						from      = paste0('tmpDir', '/', file_names), 
						to        = outputDirs[ chains[2] ],
						overwrite = TRUE
					)
				
				mh_acceptance_rate <- mh_acceptance_rate + 1
			}
		}
		for(myChain in 1:nChains){
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}

		z <- as.numeric(read.table('alpha_1/zValues.txt', colClasses = rep('numeric', n)))
		if( (iteration %% 10) == 0 ){
			if(progressGraphs == TRUE){
				par(mfrow = c(1,3))
				matplot(kValues[1:iteration, ], type = "l")
				points(1:iteration, kValues[1:iteration, 1], type = "b", col = 1)
				matplot(t(x_data), type = "l", col = mypal[as.numeric(as.factor(z))])
				matplot(t(originalX), type = "l", col = mypal[as.numeric(as.factor(z))])
			}
			ar <- round(100*mh_acceptance_rate/iteration, 3)
			cat("\r", paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.        '))
			#cat(paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.'), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				w        <- as.numeric(read.table('alpha_1/wValues.txt', colClasses = rep('numeric',Kmax))) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt', colClasses = rep('numeric',Kmax*p)))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
			if(q > 0){
				y        <- as.numeric(read.table('alpha_1/yValues.txt', colClasses = rep('numeric',n*q)))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt', colClasses = rep('numeric',q)))
				Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues.txt'), colClasses=rep('numeric',kpq) ) )
				cat(Lambda, file = LambdaConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
			}
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
				cat(ar	    , file = swapRate, '\n', append = TRUE)
			}
		}
	}
	cat('\n')
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	close(LambdaConnection_target)
	close(swapRate)
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	setwd("../")
	cat('-    DONE.','\n')
}


# new in version 3
# UUC and UCC models (sameSigma = TRUE => UCC, sameSigma = FALSE => UUC)
fabMix_UxC <- function(sameSigma = TRUE, dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle, gibbs_z = 1, warm_up_overfitting = 100, warm_up = 500, overfittingInitialization=TRUE, progressGraphs = FALSE, gwar = 0.05, lowerTriangular=TRUE){

	missingRowsIndex <- which(is.na(rowSums(rawData)) == TRUE)
	nMissingRows <- length( missingRowsIndex ) 
	if(nMissingRows > 0){
		stop("The data contains missing values. Use the `fabMix_missing_values()` function.")
	}
	if(progressGraphs == TRUE){
		dev.new(width=15, height=5)
	}
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
	}
	nChains <- length(dirPriorAlphas)
	registerDoParallel(cores = nChains)

	#get rid of other packages messages after printing logo
	nothingHere <- 0
	foreach(nothingHere=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
		nothingHere + 1
	}
	mypal <- c(brewer.pal(9, "Set1"), "black") # up to 10 colours
	if( nChains > 1 ){
		if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
		if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
	dir.create(outDir)
	setwd(outDir)
	outputDirs <- paste0(getwd(), '/alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
	if( missing(thinning) ){thinning = 1}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat('\n')
#	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	if(sameSigma == TRUE){
		cat(paste0('-    Parameterization: UCC model'),'\n')
	}else{
		cat(paste0('-    Parameterization: UUC model'),'\n')
	}
	cat(paste0("-    Number of factors: q = ", q,"\n"))
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = TRUE, scale = TRUE)
#		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		x_data <- rawData
#		cat('-    The sampler uses raw data (NOT GOOD PRACTICE).','\n')
	}
	kValues <- array(data = NA, dim = c(mCycles, nChains))
	mh_acceptance_rate <- 0
	dir.create('tmpDir')
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	kpq <- Kmax*p*q
	#	initialization
	iteration <- 1
	if(q == 0){warm_up_overfitting = 2*warm_up_overfitting}
	if(overfittingInitialization == TRUE){
		cat(paste('-    (1) Initializing from priors that lead to overfitting... '))
		d_per_cluster = 2*p + p*q + q*(q-1)/2
		if(q == 0){d_per_cluster = 10*p}
		initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_UCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_UUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = warm_up_overfitting, thinning = 1, burn = warm_up_overfitting - 1, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}
	}else{
		cat(paste('-    (1) Initializing from random starting values (NOT A GOOD PRACTICE)... '))
		d_per_cluster = 2*p + p*q + q*(q-1)/2
		initialAlphas <- dirPriorAlphas
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_UCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_UUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = 10, thinning = 1, burn = 9, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
			}
		}
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_UCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_UUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
		}
	}
	cat(paste(' OK'),'\n')
	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	LambdaConnection_target <- file(paste0(getwd(),"/LambdaValues.txt"),open = "w")  
	LambdaHeader <- character(Kmax*p*q)
	swapRate <-  file(paste0(getwd(),"/swapRate.txt"),open = "w") 
	ar <- 0
	iter <- 0
	for(k in 1:Kmax){
		for(i in 1:p){
			for(j in 1:q){
				iter <- iter + 1
				LambdaHeader[iter] <- paste0("k",k,"_i",i,"_j",j)
			}
		}
	}
	cat(LambdaHeader, file = LambdaConnection_target, '\n', append = TRUE)


	cat(paste('-    (3) Running the sampler... '),'\n')
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	for( iteration in 2:mCycles ){
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_UCC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_UUC(q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}
		if(nChains > 1){
			chains <- sample(nChains - 1, 1)
			chains <- c(chains, chains + 1)
			weights[1, ] <- as.numeric(read.table( paste0(outputDirs[ chains[1] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))
			weights[2, ] <- as.numeric(read.table( paste0(outputDirs[ chains[2] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))

			mh_denom <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[2, ] )
			mh_nom   <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[1, ] ) 
					+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[2, ] )
			mh_ratio <- mh_nom - mh_denom
			if( log(runif(1)) < mh_ratio ){
				# dir1 to tmp
				file.copy( 
						from      = paste0(outputDirs[ chains[1] ], '/', file_names), 
						to        = 'tmpDir',
						overwrite = TRUE
					)
				# dir2 to dir1
				file.copy( 
						from      = paste0(outputDirs[ chains[2] ], '/', file_names), 
						to        = outputDirs[ chains[1] ],
						overwrite = TRUE
					)
				# tmp to dir2
				file.copy( 
						from      = paste0('tmpDir', '/', file_names), 
						to        = outputDirs[ chains[2] ],
						overwrite = TRUE
					)
				
				mh_acceptance_rate <- mh_acceptance_rate + 1
			}
		}
		for(myChain in 1:nChains){
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}

		z <- as.numeric(read.table('alpha_1/zValues.txt', colClasses = rep('numeric', n)))
		if( (iteration %% 10) == 0 ){
			if(progressGraphs == TRUE){
				par(mfrow = c(1,3))
				matplot(kValues[1:iteration, ], type = "l")
				points(1:iteration, kValues[1:iteration, 1], type = "b", col = 1)
				matplot(t(x_data), type = "l", col = mypal[as.numeric(as.factor(z))])
				matplot(t(originalX), type = "l", col = mypal[as.numeric(as.factor(z))])
			}
			ar <- round(100*mh_acceptance_rate/iteration, 3)
			cat("\r", paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.        '))
			#cat(paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.'), '\n')
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				w        <- as.numeric(read.table('alpha_1/wValues.txt', colClasses = rep('numeric',Kmax))) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt', colClasses = rep('numeric',Kmax*p)))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
			if(q > 0){
				y        <- as.numeric(read.table('alpha_1/yValues.txt', colClasses = rep('numeric',n*q)))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt', colClasses = rep('numeric',q)))
				Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues.txt'), colClasses=rep('numeric',kpq) ) )
				cat(Lambda, file = LambdaConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
			}
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
				cat(ar	    , file = swapRate, '\n', append = TRUE)
			}
		}
	}
	cat('\n')
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	close(LambdaConnection_target)
	close(swapRate)

	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	setwd("../")
	cat('-    DONE.','\n')
}




fabMix_missing_values <- function(sameSigma = TRUE, dirPriorAlphas, rawData, outDir, Kmax, mCycles, burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize, thinning, zStart, nIterPerCycle, gibbs_z = 1, warm_up = 500, progressGraphs = FALSE, gwar = 0.05, lowerTriangular=TRUE){
	cat("         ____      __    __  ____     ", "\n")
	cat("        / __/___ _/ /_  /  |/  (_)  __", "\n")
	cat("       / /_/ __ `/ __ \\/ /|_/ / / |/_/", "\n")
	cat("      / __/ /_/ / /_/ / /  / / />  <  ", "\n")
	cat("     /_/  \\__,_/_.___/_/  /_/_/_/|_|  ", "\n")
	dev.new(width=15, height=5)
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if( missing(dirPriorAlphas) ){
		nChains <- 8
		dN <- 1
		dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
	}
	nChains <- length(dirPriorAlphas)
	if( nChains > 1 ){
		if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
		if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
	dir.create(outDir)
	setwd(outDir)
	registerDoParallel(cores = nChains)
	outputDirs <- paste0(getwd(), '/alpha_',1:nChains)
	originalX <- rawData
	x_data <- originalX
#	missing
	missingRowsIndex <- which(is.na(rowSums(x_data)) == TRUE)
	nMissingRows <- length( missingRowsIndex ) 
	missing_entries <- apply(x_data[missingRowsIndex, ], 1, function(y){which( is.na(y) == TRUE)})
	for(i in 1:nMissingRows){
		missing_entries[[i]] <- c(missingRowsIndex[i], missing_entries[[i]])
	}
	cat(paste0('-    Found ',nMissingRows,' rows containing missing values'),'\n')
#
	if( missing(thinning) ){thinning = 1}
	if( thinning < 1 ){ stop('thinning should be larger than or equal to 1.') }
	thinning <- floor(thinning)
	if( missing(normalize) ){normalize <- TRUE}
	cat(paste0("-    p = ", p, ", q = ", q, ", n = ",n,", g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	if(sameSigma == TRUE){
		cat(paste0('-    Parameterization: Same error variance per component'),'\n')
	}else{
		cat(paste0('-    Parameterization: Different error variance per component'),'\n')
	}
	cat(paste0('-    Using Nchains = ', nChains),'\n')
	cat(paste0('-    Target posterior distribution corresponds to alpha = ', dirPriorAlphas[1]),'\n')
	if( normalize == TRUE ){
		x_data <- scale(originalX, center = T, scale = apply(originalX, 2, sd, na.rm = TRUE))
		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		x_data <- rawData
		cat('-    The sampler uses raw data.','\n')
	}
	kValues <- array(data = NA, dim = c(mCycles, nChains))
	mh_acceptance_rate <- 0
	dir.create('tmpDir')
	v_r <- rep(q, p) #indicates the non-zero values of Lambdas
	if(lowerTriangular){
		for( r in 1:p ){
			v_r[r] <- min(r,q)
		}
	}

	kpq <- Kmax*p*q
	#	initialization
	iteration <- 1
	cat(paste('-    (1) Initializing from priors that lead to overfitting... '))
	d_per_cluster = 2*p + p*q + q*(q-1)/2
	initialAlphas <- seq(d_per_cluster/2, d_per_cluster, length = nChains)
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_Sj_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = 100, thinning = 1, burn = 99, alpha_prior= rep(initialAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = FALSE, gibbs_z = gwar, lowerTriangular=lowerTriangular)
		}
	}
	cat(paste(' OK'),'\n')
	cat(paste('-    (2) Initializing the actual model from the previously obtained values... '))
	if(sameSigma == TRUE){
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
		}
	}else{
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			overfittingMFA_Sj_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
				Kmax = Kmax, m = warm_up, thinning = 1, burn = warm_up - 1, alpha_prior= rep(dirPriorAlphas[myChain], Kmax), g = g, h = h, 
				alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, start_values = TRUE, gibbs_z = gibbs_z, lowerTriangular=lowerTriangular)
		}
	}
	cat(paste(' OK'),'\n')
	for(myChain in 1:nChains){
		kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
	}
	file_names <- list.files(outputDirs[1])	# the name of files are the same for each folder

	#	connections for saving the MCMC output corresponding to the target posterior distribution
	zConnection_target <- file(paste0(getwd(),"/zValues.txt"),open = "w")
	yConnection_target <- file(paste0(getwd(),"/yValues.txt"),open = "w")
	sigmainvConnection_target <- file(paste0(getwd(),"/sigmainvValues.txt"),open = "w")
	omegainvConnection_target <- file(paste0(getwd(),"/omegainvValues.txt"),open = "w")
	muConnection_target <- file(paste0(getwd(),"/muValues.txt"),open = "w")
	wConnection_target <- file(paste0(getwd(),"/wValues.txt"),open = "w")
	cllConnection_target <- file(paste0(getwd(),"/cllValues.txt"),open = "w")
	LambdaConnection_target <- file(paste0(getwd(),"/LambdaValues.txt"),open = "w")  
	LambdaHeader <- character(Kmax*p*q)
	iter <- 0
	for(k in 1:Kmax){
		for(i in 1:p){
			for(j in 1:q){
				iter <- iter + 1
				LambdaHeader[iter] <- paste0("k",k,"_i",i,"_j",j)
			}
		}
	}
	cat(LambdaHeader, file = LambdaConnection_target, '\n', append = TRUE)

	cat(paste('-    (3) Running the sampler... '),'\n')
	#	main loops
	weights <- array(data = NA, dim = c(2, Kmax))
	#par(mfrow = c(1, 2))
	bb <- nIterPerCycle - 1
	xMean <- array(data = 0, dim = c(n, p))
	for( iteration in 2:mCycles ){
		if(sameSigma == TRUE){
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}else{
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				overfittingMFA_Sj_missing_values(missing_entries = missing_entries, q = q, originalX = originalX, x_data = x_data, outputDirectory = outputDirs[myChain], 
					Kmax = Kmax, m = nIterPerCycle, thinning = 1, burn = bb, alpha_prior= rep( dirPriorAlphas[myChain], Kmax), g = g, h = h, 
					alpha_sigma = alpha_sigma, beta_sigma = beta_sigma,  start_values = TRUE, lowerTriangular=lowerTriangular)
				kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
			}
		}

		chains <- sample(nChains - 1, 1)
		chains <- c(chains, chains + 1)
		weights[1, ] <- as.numeric(read.table( paste0(outputDirs[ chains[1] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))
		weights[2, ] <- as.numeric(read.table( paste0(outputDirs[ chains[2] ],'/wValues.txt'), colClasses = rep('numeric',Kmax) ))

		mh_denom <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[1, ] ) 
				+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[2, ] )
		mh_nom   <- log_dirichlet_pdf( rep( dirPriorAlphas[ chains[2] ], Kmax ), weights[1, ] ) 
				+ log_dirichlet_pdf( rep( dirPriorAlphas[ chains[1] ], Kmax ), weights[2, ] )
		mh_ratio <- mh_nom - mh_denom
		if( log(runif(1)) < mh_ratio ){
#		if( 0 < 1 ){

			# dir1 to tmp
			file.copy( 
					from      = paste0(outputDirs[ chains[1] ], '/', file_names), 
					to        = 'tmpDir',
					overwrite = TRUE
				)
			# dir2 to dir1
			file.copy( 
					from      = paste0(outputDirs[ chains[2] ], '/', file_names), 
					to        = outputDirs[ chains[1] ],
					overwrite = TRUE
				)
			# tmp to dir2
			file.copy( 
					from      = paste0('tmpDir', '/', file_names), 
					to        = outputDirs[ chains[2] ],
					overwrite = TRUE
				)
			
			mh_acceptance_rate <- mh_acceptance_rate + 1
		}

		for(myChain in 1:nChains){
			kValues[iteration, myChain] <- read.table( paste0(outputDirs[myChain],'/k.and.logl.Values.txt') )[1,1]
		}

		z <- as.numeric(read.table('alpha_1/zValues.txt', colClasses = rep('numeric', n)))
		if( (iteration %% 10) == 0 ){
			if(progressGraphs == TRUE){
				par(mfrow = c(1,3))
				matplot(kValues[1:iteration, ], type = "l")
				points(1:iteration, kValues[1:iteration, 1], type = "b", col = 1)
				matplot(t(x_data), type = "l", col = as.numeric(as.factor(z)))
				matplot(t(originalX), type = "l", col = as.numeric(as.factor(z)))
			}

			ar <- round(100*mh_acceptance_rate/iteration, 3)
#			cat(paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.'), '\n')
			cat("\r", paste0('-        mcmc cycle: ',iteration,' (<=> iteration: ',iteration*nIterPerCycle,'). Swap acceptance rate: ', ar, '%.        '))
		}
		if(iteration %% thinning == 0){
			if(iteration > burnCycles){
				y        <- as.numeric(read.table('alpha_1/yValues.txt', colClasses = rep('numeric',n*q)))
				w        <- as.numeric(read.table('alpha_1/wValues.txt', colClasses = rep('numeric',Kmax))) 
				mu       <- as.numeric(read.table('alpha_1/muValues.txt', colClasses = rep('numeric',Kmax*p)))
				omegainv <- as.numeric(read.table('alpha_1/omegainvValues.txt', colClasses = rep('numeric',q)))
				sigmainv <- as.numeric(read.table('alpha_1/sigmainvValues.txt')) 
				cll      <- as.numeric(read.table('alpha_1/k.and.logl.Values.txt') )[2]
				xCurrent <- as.matrix(read.table('alpha_1/x_complete.txt', header = TRUE) )
				xMean    <- xMean + xCurrent
				Lambda <- as.numeric( read.table(paste0('alpha_1/LambdaValues.txt'), colClasses=rep('numeric',kpq) ) )
				cat(Lambda, file = LambdaConnection_target, '\n', append = TRUE)
				cat(z       , file = zConnection_target, '\n', append = TRUE)
				cat(y       , file = yConnection_target, '\n', append = TRUE)
				cat(w       , file = wConnection_target, '\n', append = TRUE)
				cat(mu      , file = muConnection_target, '\n', append = TRUE)
				cat(omegainv, file = omegainvConnection_target, '\n', append = TRUE)
				cat(sigmainv, file = sigmainvConnection_target, '\n', append = TRUE)
				cat(cll     , file = cllConnection_target, '\n', append = TRUE)
			}
		}
	}
	stopImplicitCluster()
	close(zConnection_target)
	close(yConnection_target)
	close(wConnection_target)
	close(muConnection_target)
	close(sigmainvConnection_target)
	close(omegainvConnection_target)
	close(cllConnection_target)
	close(LambdaConnection_target)
	keepedSeq <- seq(burnCycles + thinning, mCycles, by = thinning)
	if( burnCycles > 0){
		burnedSeq <- 1:burnCycles
		write.table(file = 'burn_in_period_k.txt', kValues[burnedSeq, ])
	}
	write.table(file = 'kValues.txt', kValues[keepedSeq, ], quote = FALSE, row.names = FALSE, col.names = dirPriorAlphas)
	write.table(file = 'xMeanEstimate.txt', xMean/(mCycles - burnCycles), quote = FALSE, row.names = FALSE)
	setwd("../")
	cat('-    DONE.','\n')
}




observed.log.likelihood0 <- function(x_data, w, mu, Lambda, Sigma, z){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	ct <- -(p/2)*log(2*pi)

        probs <- numeric(n)
        if( is.null(z) ){alive <- 1}else{
        alive <- as.numeric(names(table(z)))}
        newW <- numeric(length(w))
        newW[alive] <- w[alive]/sum(w[alive])
        loggedValues <- array(data = NA, dim = c(n, length(alive)))
        colnames(loggedValues) <- alive
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
                diag(x_var) <- diag(x_var) + Sigma
#                x_var <- try(solve(x_var), TRUE)
#                loggedValues[ ,as.character(k)] <- log(newW[k]) - 0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) + ct
                loggedValues[ ,as.character(k)] <- log(newW[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
        }
        lMax <- apply(loggedValues, 1, max)
        if( length(alive) == 1 ){
                logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
        }else{
                logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
        }
        return( logL )
}


observed.log.likelihood0_Sj <- function(x_data, w, mu, Lambda, Sigma, z){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	ct <- -(p/2)*log(2*pi)

        probs <- numeric(n)
        if( is.null(z) ){alive <- 1}else{
        alive <- as.numeric(names(table(z)))}
        newW <- numeric(length(w))
        newW[alive] <- w[alive]/sum(w[alive])
        loggedValues <- array(data = NA, dim = c(n, length(alive)))
        colnames(loggedValues) <- alive
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                x_var <- Lambda[k,,] %*% t(Lambda[k,,]) 
                diag(x_var) <- diag(x_var) + Sigma[k,]
#                x_var <- try(solve(x_var), TRUE)
#                loggedValues[ ,as.character(k)] <- log(newW[k]) - 0.5*apply(center_x,1,function(tmp){return( as.numeric(t(tmp) %*% x_var %*% tmp) )}) + 0.5*log(det(x_var)) + ct
                loggedValues[ ,as.character(k)] <- log(newW[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
        }
        lMax <- apply(loggedValues, 1, max)
        if( length(alive) == 1 ){
                logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
        }else{
                logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
        }
        return( logL )
}


observed.log.likelihood0_Sj_q0 <- function(x_data, w, mu, Sigma, z){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	ct <- -(p/2)*log(2*pi)

        probs <- numeric(n)
        if( is.null(z) ){alive <- 1}else{
        alive <- as.numeric(names(table(z)))}
        newW <- numeric(length(w))
        newW[alive] <- w[alive]/sum(w[alive])
        loggedValues <- array(data = NA, dim = c(n, length(alive)))
        colnames(loggedValues) <- alive
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                x_var <- array(data = 0, dim = c(p,p)) 
                diag(x_var) <- diag(x_var) + Sigma[k,]
                loggedValues[ ,as.character(k)] <- log(newW[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
        }
        lMax <- apply(loggedValues, 1, max)
        if( length(alive) == 1 ){
                logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
        }else{
                logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
        }
        return( logL )
}


observed.log.likelihood0_q0_sameSigma <- function(x_data, w, mu, Sigma, z){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	ct <- -(p/2)*log(2*pi)

        probs <- numeric(n)
        if( is.null(z) ){alive <- 1}else{
        alive <- as.numeric(names(table(z)))}
        newW <- numeric(length(w))
        newW[alive] <- w[alive]/sum(w[alive])
        loggedValues <- array(data = NA, dim = c(n, length(alive)))
        colnames(loggedValues) <- alive
        x_var <- array(data = 0, dim = c(p,p)) 
        diag(x_var) <- Sigma
        for(k in alive){
                center_x <- x_data - matrix(mu[k,], nrow = n, ncol = p, byrow=TRUE)
                loggedValues[ ,as.character(k)] <- log(newW[k]) + dmvnorm(center_x, mean = rep(0, p), sigma = x_var, log = TRUE)
        }
        lMax <- apply(loggedValues, 1, max)
        if( length(alive) == 1 ){
                logL <- sum(lMax + log( exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ) ) )
        }else{
                logL <- sum(lMax + log( rowSums(exp( apply(loggedValues, 2, function(y){return(y - lMax)}) ))))
        }
        return( logL )
}



getStuffForDIC <- function(sameSigma = TRUE, sameLambda = FALSE, isotropic  = FALSE, x_data, outputFolder, q, burn, Km, normalize, discardLower){
	cat(paste0('-    (4) Computing information criteria for q = ', q), '\n')
	if(missing(normalize)){normalize = TRUE}
	if(normalize){
		x_data <- scale(x_data, center = TRUE, scale = TRUE)
		cat('-    NOTE: using standardized data.','\n')
	}
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	if(missing(Km)){Km <- 20}
	if(missing(burn)){burn <- 0}
	setwd(outputFolder)
	cat(paste0('         - Entering directory: ', getwd()),'\n')
	z <- as.matrix(read.table("zValues.txt", colClasses = rep('numeric', n)))
	logl <- read.table("cllValues.txt")
        tmp  <- apply(z,1,function(y){length(table(y))})
        logl <- cbind(tmp, logl)
	if(burn > 0){
		z <- z[-(1:burn),]
		logl <- logl[-(1:burn),]
	}
	cat(paste0('            Nclusters:    ', paste(as.character(names(table(logl[,1]))), collapse="     ") ), '\n')
	cat(paste0('            Frequency:    ', paste(as.character(as.numeric(table(logl[,1]))), collapse="    ") ), '\n')
#	K <- Km
#	kSelected <- K
#	index <- 1:dim(z)[1]
#	Kindex <- index
        K <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
        kSelected <- K
        index <- which(logl[,1] == K)
        Kindex <- index
	m <- length(index)

	#this is artificial

	ECR <- matrix(1:Km, nrow = m, ncol = Km, byrow=T)
	permutations <- vector('list', length = 1)
	permutations[[1]] <- ECR
	names(permutations) <- "ECR"
	ls <- vector('list', length = 1)
	names(ls) <- "permutations"
	ls$permutations <- permutations	


	if(q > 0){
		l <- read.table("LambdaValues.txt", header=TRUE)
		if(burn > 0){
			l <- l[-(1:burn),]
		}
		J <- p*q
		mcmc <- array(data = NA, dim = c(m,Km,J))
		t <- 0
		i <- 1:p
		for(iter in index){
			t <- t + 1
			for(k in 1:Km){
#				for(i in 1:p){
					for(j in 1:q){
						mcmc[t, k, (i-1)*q + j] <- as.numeric(l[iter, paste0("k",k,"_i",i,"_j",j)])
					}
#				}
			}
		}
		lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
		for(k in 1:Km){
			lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
		}
	}
	mu <- read.table("muValues.txt", colClasses = rep('numeric',Km*p))# auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
	if(burn > 0){
		mu <- mu[-(1:burn),] 
	}
	mu <- mu[Kindex,]
	mu.mcmc <- array(data = NA, dim = c(m,Km,p))
	for(k in 1:Km){
		mu.mcmc[,k,] <- as.matrix(mu[,k + Km*((1:p)-1)])
	}
	mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
	#
	if(sameSigma == TRUE){
		SigmaINV <- as.matrix(read.table("sigmainvValues.txt"))# auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
		if(burn > 0){
			SigmaINV <- SigmaINV[-(1:burn),] 
		}
		SigmaINV <- SigmaINV[Kindex, ] 
		SigmaINV.mcmc <- SigmaINV
	}else{
		SigmaINV <- as.matrix(read.table("sigmainvValues.txt")) # auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
		if(burn > 0){
			SigmaINV <- SigmaINV[-(1:burn),] 
		}
		SigmaINV <- SigmaINV[Kindex,  ] 
		SigmaINV.mcmc <- array(data = NA, dim = c(m,Km,p))
		for(k in 1:Km){
		        SigmaINV.mcmc[,k,] <- as.matrix(SigmaINV[,((k-1)*p + 1):(k*p)])
		}
		SigmaINV.mcmc <- permute.mcmc(SigmaINV.mcmc, ls$permutations$ECR)$output

	}
	Sigma.mcmc <- 1/SigmaINV.mcmc

	#SigmaINV.mean <- as.numeric(apply(SigmaINV,2,mean))
	w.mcmc <- as.matrix(read.table("wValues.txt", colClasses = rep('numeric',Km)))
	w.mcmc <- array(w.mcmc, dim = c(dim(w.mcmc)[1], Km, 1))
	if(burn > 0){
		w.mcmc <- w.mcmc[-(1:burn),,]
		w.mcmc <- w.mcmc[Kindex,]
	}else{
		w.mcmc <- w.mcmc[Kindex,,]
	}
	w.mcmc <- array(w.mcmc[,1:Km],dim = c(length(Kindex),Km,1))
	w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$"ECR")$output

	lValues <- numeric(m)
	cll <- 0
	aic <- 0
	bic <- 0
	maxL <- 1
	i <- 1
	if(sameSigma == TRUE){
		Sigma.current <- Sigma.mcmc[i, ]
	}else{
	        Sigma.current <- Sigma.mcmc[i, , ]		
	}
	mu.current <- mu.mcmc[i,,]
	if( Km == 1 ){ 	
		Sigma.current <- array(Sigma.current, dim = c(1, p))
		mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
	}
	if(q > 0){lambda.current <- array(data = NA, dim = c(Km,p,q))}
	for(k in 1:Km){
		if(q > 0){
			ldraw <- lambda.perm.mcmc[i,k,]
			lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
		}
		if(sameSigma == TRUE){
			for(i1 in 1:p){
				if( Sigma.current[ i1] > 1000 ){  Sigma.current[i1] <- 1000 }
			}
		}else{
			for(i1 in 1:p){
				if( Sigma.current[k, i1] > 1000 ){ Sigma.current[k, i1] <- 1000 }
			}
		}
	}
	if(sameSigma == TRUE){
		if(q > 0){
			obsL <- observed.log.likelihood0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
		}else{
			obsL <- observed.log.likelihood0_q0_sameSigma(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
		}
	}else{
		if(q > 0){
			obsL <- observed.log.likelihood0_Sj(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
		}else{
			obsL <- observed.log.likelihood0_Sj_q0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
		}
	}
	lValues[i] <- obsL
	maxL <- obsL
	cll <- cll + obsL
	aic <- aic + obsL
	bic <- bic + obsL
	iterMax <- i
	if(sameSigma == TRUE){
		for(i in 2:m){
	#		cat(paste0("i  = ", i), "\n")
			lambda.current <- array(data = NA, dim = c(Km,p,q))
			Sigma.current <- Sigma.mcmc[i, ]
			mu.current <- mu.mcmc[i,,]
			if( Km == 1 ){ 	
				Sigma.current <- array(Sigma.current, dim = c(1, p))
				mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
			}
			for(k in 1:Km){
	#			cat(paste0("  k  = ", k), "\n")
				if(q > 0){
					ldraw <- lambda.perm.mcmc[i,k,]
					lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
				}
				for(i1 in 1:p){
					if( Sigma.current[ i1] > 1000 ){ 
#						cat(paste0('oops: ', i),'\n'); 
					Sigma.current[ i1] <- 1000 }
				}
			}
			if(q > 0){
				obsL <- observed.log.likelihood0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
			}else{
				obsL <- observed.log.likelihood0_q0_sameSigma(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
			}
			lValues[i] <- obsL
			if( obsL > maxL ){
				maxL <- obsL	
				iterMax <- i
			}
			cll <- cll + obsL
			aic <- aic + obsL
			bic <- bic + obsL
		}
	}else{
		for(i in 2:m){
	#		cat(paste0("i  = ", i), "\n")
			lambda.current <- array(data = NA, dim = c(Km,p,q))
			Sigma.current <- Sigma.mcmc[i, , ]
			mu.current <- mu.mcmc[i,,]
			if( Km == 1 ){ 	
				Sigma.current <- array(Sigma.current, dim = c(1, p))
				mu.current <- array(mu.mcmc[i,,], dim = c(1, p))
			}
			for(k in 1:Km){
	#			cat(paste0("  k  = ", k), "\n")
				if(q > 0){
					ldraw <- lambda.perm.mcmc[i,k,]
					lambda.current[k,,] <- matrix(ldraw,nrow = p, ncol = q, byrow=TRUE)
				}
				for(i1 in 1:p){
					if( Sigma.current[k, i1] > 1000 ){ Sigma.current[k, i1] <- 1000 }
				}
			}
			if(q > 0){
				obsL <- observed.log.likelihood0_Sj(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Lambda = lambda.current, Sigma = Sigma.current, z = z[i,])
			}else{
				obsL <- observed.log.likelihood0_Sj_q0(x_data = x_data, w = w.mcmc[i,,1], mu = mu.current, Sigma = Sigma.current, z = z[i,])
			}
			lValues[i] <- obsL
			if( obsL > maxL ){
				maxL <- obsL	
				iterMax <- i
			}
			cll <- cll + obsL
			aic <- aic + obsL
			bic <- bic + obsL
		}
	}
	if(missing(discardLower)){ discardLower <- 0.01 }
	if ( discardLower == FALSE){
		cll <- cll/m
	}else{
		cll <- mean( lValues[which( lValues > as.numeric(quantile(lValues, discardLower)) )] )
	}
	dic_classicMAP <- -4*cll + 2*maxL
	dic_starMAP <- -6*cll + 4*maxL
	dic_classic <- dic_classicMAP
	dic_star <- dic_starMAP

	if(sameSigma == TRUE){
		if(sameLambda == FALSE){
			if(isotropic == FALSE){
				#UCU
				aic <- -2*aic/m + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + p + kSelected - 1 )
			}else{
				#UCC
				aic <- -2*aic/m + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + 1 + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + 1 + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*( p+p*q - q*(q-1)/2 ) + 1 + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*( p+p*q - q*(q-1)/2 ) + 1 + kSelected - 1 )
			}
		}else{
			if(isotropic == FALSE){
				#CCU
				aic <- -2*aic/m + 2*(kSelected*p + p*q - q*(q-1)/2 + p + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*p + p*q - q*(q-1)/2 + p + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*p + p*q - q*(q-1)/2 + p + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*p + p*q - q*(q-1)/2 + p + kSelected - 1 )
			}else{
				#CCC
				aic <- -2*aic/m + 2*(kSelected*p + p*q - q*(q-1)/2 + 1 + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*p + p*q - q*(q-1)/2 + 1 + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*p + p*q - q*(q-1)/2 + 1 + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*p + p*q - q*(q-1)/2 + 1 + kSelected - 1 )
			}
		}
	}else{
		if(sameLambda == FALSE){
			if(isotropic == FALSE){
				#UUU
				aic <- -2*aic/m + 2*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*( 2*p+p*q - q*(q-1)/2 ) + kSelected - 1 )
			}else{
				#UUC
				aic <- -2*aic/m + 2*(kSelected*( p + 1 + p*q - q*(q-1)/2 ) + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*( p + 1 + p*q - q*(q-1)/2 ) + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*( p + 1 + p*q - q*(q-1)/2 ) + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*( p + 1 + p*q - q*(q-1)/2 ) + kSelected - 1 )
			}
		}else{
			if(isotropic == FALSE){
				#CUU
				aic <- -2*aic/m + 2*(kSelected*2*p + p*q - q*(q-1)/2 + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*2*p + p*q - q*(q-1)/2 + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*2*p + p*q - q*(q-1)/2 + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*2*p + p*q - q*(q-1)/2 + kSelected - 1 )
			}else{
				#CUC
				aic <- -2*aic/m + 2*(kSelected*(p + 1) + p*q - q*(q-1)/2 + kSelected - 1 ) 
				bic <- -2*bic/m + log(n)*(kSelected*(p + 1) + p*q - q*(q-1)/2 + kSelected - 1 )
				aic_MAX <- -2*maxL + 2*(kSelected*(p + 1) + p*q - q*(q-1)/2 + kSelected - 1 ) 
				bic_MAX <- -2*maxL + log(n)*(kSelected*(p + 1) + p*q - q*(q-1)/2 + kSelected - 1 )
			}
		}
	}

	dic <- c(aic, bic, dic_classic, dic_star, dic_classicMAP, dic_starMAP, aic_MAX, bic_MAX)
	names(dic) <- c('AIC', 'BIC', 'DIC1', 'DIC*2', 'DIC', 'DIC_2', 'AIC_map', 'BIC_map')	
	write.table(file = 'informationCriteria_map_model.txt', dic[c(5,6,7,8)], col.names = paste0('q_',q), quote = FALSE)
	write.table(file = 'lValues_map.txt', lValues, quote = FALSE)
	setwd("../")
	cat(paste0('         - Information criteria written to `', outputFolder,'/informationCriteria_map_model.txt`.'), '\n')
}



dealWithLabelSwitching <- function(sameSigma = TRUE, x_data, outputFolder, q, burn, z.true, compute_regularized_expression, Km){
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	if(missing(Km)){Km <- 20}
	if(missing(burn)){burn <- 0}
	if(missing(compute_regularized_expression)){ compute_regularized_expression = FALSE }
	cat(paste0('-    (5) Dealing with label switching for q = ', q), '\n')
	setwd(outputFolder)
	cat(paste0('         * Entering directory: ', getwd()),'\n')
	z <- as.matrix(read.table("zValues.txt", colClasses = rep('numeric', n)))
	logl <- read.table("kValues.txt", header=T)
	cllValues <- read.table("cllValues.txt")
	if(burn > 0){
		logl <- logl[-(1:burn), ]
		cllValues <- cllValues[-(1:burn), ]
		z <- z[-(1:burn),]
	}

	K <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
	kSelected <- K
	cat(paste0('         * Posterior mode corresponds to K = ', K),'\n')
	if(K == 1){
		cat(paste0('         *  no label switching algorithms are applied. Processing output...'),'\n')
		index <- which(logl[,1] == K)
		Kindex <- index
		logl <- logl[index,1]
		cllValues <- cllValues[index,1]
		z <- z[index,]
		mapIndex <- which(cllValues == max(cllValues))[1]
		m <- length(logl)

		K <- max(z)
		skiniko <- FALSE
		if(K < Km){
			skiniko <- TRUE
			kLAST <- K 
			newZ <- t(apply(z, 1, function(y){ myI <- which(y == K); y[myI] <- rep(Km, length(myI));return(y) } ) )
			z <- newZ
			K <- max(z)
		}

		write.table(matrix(1,nrow = m, ncol = n), file = 'reordered_allocations_ecr.txt')
		sbc <- matrix(1,nrow = 2, ncol = n)
		colnames(sbc) <- paste0('V',1:n)
		write.table(sbc, file = "singleBestClusterings.txt", quote = FALSE, row.names = FALSE)
		write.table(rep(1,n), file = "classificationProbabilities.txt")


		l <- read.table("LambdaValues.txt", header=TRUE)
		if(burn > 0){
			l <- l[-(1:burn),]
		}
		J <- p*q
		mcmc <- array(data = NA, dim = c(m,K,J))
		t <- 0
		i <- 1:p
		for(iter in index){
			t <- t + 1
			for(k in 1:K){
				for(j in 1:q){
					mcmc[t, k, (i-1)*q + j] <- as.numeric(l[iter, paste0("k",k,"_i",i,"_j",j)])
				}
			}
		}
		if(skiniko == TRUE){
			tmp1 <- mcmc[,kLAST,]
			tmp2 <-  mcmc[,Km,]	
			mcmc[,kLAST,] <- tmp2
			mcmc[,Km,] <- tmp1
		}
		lambda.perm.mcmc <- mcmc
		t <- 0
		for(iter in index){
			t <- t + 1
			j <- z[t, 1]
			lambda.perm.mcmc[t, 1, ] <- mcmc[t, j, ]
			lambda.perm.mcmc[t, j, ] <- mcmc[t, 1, ]
		}
		
		lCon <- file('reordered_lambda_ecr.txt', open = "w")
		cat(colnames(l), file = lCon, '\n', append = TRUE)
		for (i in 1:m){
			for(k in 1:K){
				cat(lambda.perm.mcmc[i,k, ], " ", file = lCon, append = TRUE)
			}
			cat("\n", file = lCon, append = TRUE)
		}
		close(lCon)
		lambda.map <- array(data = NA, dim = c(K,p,q))
		for(k in 1:K){
			lambda.map[k,,] <- matrix(lambda.perm.mcmc[mapIndex,k,],nrow = p, ncol = q, byrow=TRUE)
		}
		write.table(lambda.map, file = 'lambda_map.txt', col.names = paste('lambda',1:p, rep(1:q, each = p), sep = "_"))


		mu <- read.table("muValues.txt", colClasses = rep('numeric',Km*p))# auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
		if(burn > 0){
			mu <- mu[-(1:burn),] 
		}
		mu <- mu[Kindex,]

		mu.mcmc <- array(data = NA, dim = c(m,K,p))
		for(k in 1:K){
			mu.mcmc[,k,] <- as.matrix(mu[,k + K*((1:p)-1)])
		}
		if(skiniko == TRUE){
			tmp1 <- mu.mcmc[,kLAST,]
			tmp2 <-  mu.mcmc[,Km,]	
			mu.mcmc[,kLAST,] <- tmp2
			mu.mcmc[,Km,] <- tmp1
		}
		tmp1 <- tmp2 <- 0
		mu.perm.mcmc <- mu.mcmc
		t <- 0
		for(iter in index){
			t <- t + 1
			j <- z[t, 1]
			mu.perm.mcmc[t, 1, ] <- mu.mcmc[t, j, ]
			mu.perm.mcmc[t, j, ] <- mu.mcmc[t, 1, ]
		}
		write.table(mu.perm.mcmc, file = 'reordered_mu_ecr.txt')
		mu.mean <- array(data = NA, dim = c(K,p))
		mu.map <- array(data = NA, dim = c(K,p))
		for(k in 1:K){
			for(j in 1:p){
				mu.mean[k,j] <- mean(mu.perm.mcmc[,k,j])
				mu.map[k,j] <- mu.mcmc[mapIndex,k,j]
			}
		}
		write.table(mu.mean, file = 'mu_estimate_ecr.txt')
		if(sameSigma == TRUE){
			sigmaINV <- as.matrix(read.table("sigmainvValues.txt"))# auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
			if(burn > 0){
				sigmaINV <- sigmaINV[-(1:burn),] 
			}
			sigmaINV <- sigmaINV[Kindex, ] 
		        Sigma.mcmc <- 1/sigmaINV
		        write.table(Sigma.mcmc, file = 'reordered_sigma_ecr.txt')
		}else{	
		
		        SigmaINV <- as.matrix(read.table("sigmainvValues.txt"))
			if(burn > 0){
				SigmaINV <- SigmaINV[-(1:burn),] 
			}
			SigmaINV <- SigmaINV[Kindex,]

		        SigmaINV.mcmc <- array(data = NA, dim = c(m,K,p))
		        for(k in 1:K){
		                SigmaINV.mcmc[,k,] <- as.matrix(SigmaINV[,((k-1)*p + 1):(k*p)])
		        }
			if(skiniko == TRUE){
				tmp1 <- SigmaINV.mcmc[,kLAST,]
				tmp2 <-  SigmaINV.mcmc[,Km,]	
				SigmaINV.mcmc[,kLAST,] <- tmp2
				SigmaINV.mcmc[,Km,] <- tmp1
			}
			tmp1 <- tmp2 <- 0
		        SigmaINV.perm.mcmc <- SigmaINV.mcmc
		        Sigma.mcmc <- 1/SigmaINV.mcmc
			Sigma.perm.mcmc <- Sigma.mcmc
			t <- 0
			for(iter in index){
				t <- t + 1
				j <- z[t, 1]
				Sigma.perm.mcmc[t, 1, ] <- Sigma.mcmc[t, j, ]
				Sigma.perm.mcmc[t, j, ] <- Sigma.mcmc[t, 1, ]
			}
		        write.table(Sigma.perm.mcmc, file = 'reordered_sigma_ecr.txt')
		        sigma.mean <- array(data = NA, dim = c(K,p))
		        sigma.map <- array(data = NA, dim = c(K,p))
		        for(k in 1:K){
		                for(j in 1:p){
		                        sigma.mean[k,j] <- mean(Sigma.perm.mcmc[,k,j])
		                        sigma.map[k,j] <- Sigma.mcmc[mapIndex,k,j]
		                }
		        }
		        write.table(sigma.mean, file = 'sigma_estimate_ecr.txt')
		}

		w.mcmc <- as.matrix(read.table("wValues.txt", colClasses = rep('numeric',Km)))
		w.mcmc <- array(w.mcmc, dim = c(dim(w.mcmc)[1], K, 1))
		if(burn > 0){
			w.mcmc <- w.mcmc[-(1:burn),,]
			w.mcmc <- w.mcmc[Kindex,]
		}else{
			w.mcmc <- w.mcmc[Kindex,,]
		}
		w.mcmc <- array(w.mcmc[,1:K],dim = c(length(Kindex),K,1))
		w.mcmc_raw <- w.mcmc
		if(skiniko == TRUE){
			tmp1 <- w.mcmc[,kLAST,]
			tmp2 <-  w.mcmc[,Km,]	
			w.mcmc[,kLAST,] <- tmp2
			w.mcmc[,Km,] <- tmp1
		}
		w.perm.mcmc <- w.mcmc
		t <- 0
		for(iter in index){
			t <- t + 1
			j <- z[t, 1]
			w.perm.mcmc[t, 1, ] <- w.mcmc[t, j, ]
			w.perm.mcmc[t, j, ] <- w.mcmc[t, 1, ]
		}

		w.mean <- numeric(K)
		w.map <- numeric(K)
		for(k in 1:K){
			w.mean[k] <- mean(w.perm.mcmc[,k,1])
			w.map[k] <- w.mcmc[mapIndex,k,1]
		}
		write.table(w.perm.mcmc, file = 'reordered_weights_ecr.txt')
		write.table(w.mean, file = 'weights_estimate_ecr.txt')

		aliveClusters <- 1
		covmat <- array(data = 0, dim = c( length(aliveClusters), p, p ))
		rownames(covmat) <- as.character(aliveClusters)
		if(sameSigma == TRUE){
			for(k in aliveClusters){
				for(iter in 1:m){
					lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
					covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
				}
				covmat[as.character(k), , ] <- covmat[as.character(k), , ]/m
				diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(apply(1/sigmaINV, 2, median))
			}
		}else{
			for(k in aliveClusters){
				for(iter in 1:m){
					lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
					covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
				}
				covmat[as.character(k), , ] <- covmat[as.character(k), , ]/m
			
				diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(apply( Sigma.perm.mcmc[ ,k, ], 2, median ))
			}
		}
		for(k in 1:length(aliveClusters)){
			write.table(covmat[k, , ], file = paste0("estimated_cov_cluster_",k,".txt"))
			write.table(lambda.map[k, , ], file = paste0('lambda_map_',k,'.txt'))
		}



		cat(paste0('-    Done.'), '\n')
	}else{
		index <- which(logl[,1] == K)
		Kindex <- index
		logl <- logl[index,1]
		cllValues <- cllValues[index,1]
		z <- z[index,]
		mapIndex <- which(cllValues == max(cllValues))[1]
		zPivot <- as.numeric(z[mapIndex,])
		K <- max(z)
		skiniko <- FALSE
		if(K < Km){
			skiniko <- TRUE
			kLAST <- K 
			newZ <- t(apply(z, 1, function(y){ myI <- which(y == K); y[myI] <- rep(Km, length(myI));return(y) } ) )
			z <- newZ
			newZ <- zPivot
			myI <- which(newZ == K)
			zPivot[myI] <- rep(Km, length(myI))
			K <- max(z)
		}
		m <- length(logl)
		if(missing(z.true)){
			ls <- label.switching(method = c("ECR", "ECR-ITERATIVE-1"), zpivot = zPivot, z = z, K = K)}else{
			ls <- label.switching(method = c("ECR", "ECR-ITERATIVE-1"), zpivot = zPivot, z = z, K = K, groundTruth=z.true)
		}

		lsMethod <- "ECR"
		if(  kSelected != length(table(ls$clusters[1,])) ){
			cat(paste0('         * WARNING: `ECR-ITERATIVE-1` selected as label switching method.'),'\n')
			lsMethod <- "ECR-ITERATIVE-1"
		}

		oldLabels <- 1:K
		index <- Kindex
		allocationsECR <- z
		for (i in 1:m){
			myPerm <- order(ls$permutations[[lsMethod]][i,])
			allocationsECR[i,] <- myPerm[z[i,]]
		}
		write.table(allocationsECR, file = 'reordered_allocations_ecr.txt')

		if(q > 0){
			l <- read.table("LambdaValues.txt", header=TRUE)
			if(burn > 0){
				l <- l[-(1:burn),]
			}
			J <- p*q
			mcmc <- array(data = NA, dim = c(m,K,J))
			t <- 0
			i <- 1:p
			for(iter in index){
				t <- t + 1
				for(k in 1:K){
					for(j in 1:q){
						mcmc[t, k, (i-1)*q + j] <- as.numeric(l[iter, paste0("k",k,"_i",i,"_j",j)])
					}
				}
			}

			if(skiniko == TRUE){
				tmp1 <- mcmc[,kLAST,]
				tmp2 <-  mcmc[,Km,]	
				mcmc[,kLAST,] <- tmp2
				mcmc[,Km,] <- tmp1
			}
			lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations[[lsMethod]])$output
			lCon <- file('reordered_lambda_ecr.txt', open = "w")
			cat(colnames(l), file = lCon, '\n', append = TRUE)
			for (i in 1:m){
				for(k in 1:K){
					cat(lambda.perm.mcmc[i,k, ], " ", file = lCon, append = TRUE)
				}
				cat("\n", file = lCon, append = TRUE)
			}
			close(lCon)

			lambda.mean <- array(data = NA, dim = c(K,p,q))
			lambda.map <- array(data = NA, dim = c(K,p,q))
			for(k in 1:K){
				lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
				lambda.mean[k,,] <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
				lambda.map[k,,] <- matrix(lambda.perm.mcmc[mapIndex,k,],nrow = p, ncol = q, byrow=TRUE)
			}
			write.table(lambda.map, file = 'lambda_map.txt', col.names = paste('lambda',1:p, rep(1:q, each = p), sep = "_"))
		}
		mu <- read.table("muValues.txt", colClasses = rep('numeric',Km*p))# auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
		if(burn > 0){
			mu <- mu[-(1:burn),] 
		}
		mu <- mu[Kindex,]

		mu.mcmc <- array(data = NA, dim = c(m,K,p))
		for(k in 1:K){
			mu.mcmc[,k,] <- as.matrix(mu[,k + K*((1:p)-1)])
		}
		if(skiniko == TRUE){
			tmp1 <- mu.mcmc[,kLAST,]
			tmp2 <-  mu.mcmc[,Km,]	
			mu.mcmc[,kLAST,] <- tmp2
			mu.mcmc[,Km,] <- tmp1
		}
		mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations[[lsMethod]])$output
		write.table(mu.mcmc, file = 'reordered_mu_ecr.txt')
		mu.mean <- array(data = NA, dim = c(K,p))
		mu.map <- array(data = NA, dim = c(K,p))
		for(k in 1:K){
			for(j in 1:p){
				mu.mean[k,j] <- mean(mu.mcmc[,k,j])
				mu.map[k,j] <- mu.mcmc[mapIndex,k,j]
			}
		}
		write.table(mu.mean, file = 'mu_estimate_ecr.txt')
		if(sameSigma == TRUE){
			sigmaINV <- as.matrix(read.table("sigmainvValues.txt"))# auto to grafei ws (s_{11},...,s_{p1}),....,(s_{1k},...,s_{pk}),....,(s_{1K},...,s_{pK})
			if(burn > 0){
				sigmaINV <- sigmaINV[-(1:burn),] 
			}
			sigmaINV <- sigmaINV[Kindex, ] 
		        Sigma.mcmc <- 1/sigmaINV
		        write.table(Sigma.mcmc, file = 'reordered_sigma_ecr.txt')
		}else{	
		
		        SigmaINV <- as.matrix(read.table("sigmainvValues.txt"))
			if(burn > 0){
				SigmaINV <- SigmaINV[-(1:burn),] 
			}
			SigmaINV <- SigmaINV[Kindex,]

		        SigmaINV.mcmc <- array(data = NA, dim = c(m,K,p))
		        for(k in 1:K){
		                SigmaINV.mcmc[,k,] <- as.matrix(SigmaINV[,((k-1)*p + 1):(k*p)])
		        }
			if(skiniko == TRUE){
				tmp1 <- SigmaINV.mcmc[,kLAST,]
				tmp2 <-  SigmaINV.mcmc[,Km,]	
				SigmaINV.mcmc[,kLAST,] <- tmp2
				SigmaINV.mcmc[,Km,] <- tmp1
			}
		        SigmaINV.mcmc <- permute.mcmc(SigmaINV.mcmc, ls$permutations[[lsMethod]])$output
		        Sigma.mcmc <- 1/SigmaINV.mcmc
		        write.table(Sigma.mcmc, file = 'reordered_sigma_ecr.txt')
		        sigma.mean <- array(data = NA, dim = c(K,p))
		        sigma.map <- array(data = NA, dim = c(K,p))
		        for(k in 1:K){
		                for(j in 1:p){
		                        sigma.mean[k,j] <- mean(Sigma.mcmc[,k,j])
		                        sigma.map[k,j] <- Sigma.mcmc[mapIndex,k,j]
		                }
		        }
		        write.table(sigma.mean, file = 'sigma_estimate_ecr.txt')
		}

		w.mcmc <- as.matrix(read.table("wValues.txt", colClasses = rep('numeric',Km)))
		w.mcmc <- array(w.mcmc, dim = c(dim(w.mcmc)[1], K, 1))
		if(burn > 0){
			w.mcmc <- w.mcmc[-(1:burn),,]
			w.mcmc <- w.mcmc[Kindex,]
		}else{
			w.mcmc <- w.mcmc[Kindex,,]
		}
		w.mcmc <- array(w.mcmc[,1:K],dim = c(length(Kindex),K,1))
		w.mcmc_raw <- w.mcmc
		if(skiniko == TRUE){
			tmp1 <- w.mcmc[,kLAST,]
			tmp2 <-  w.mcmc[,Km,]	
			w.mcmc[,kLAST,] <- tmp2
			w.mcmc[,Km,] <- tmp1
		}
		w.mcmc <- permute.mcmc(w.mcmc, ls$permutations[[lsMethod]])$output
		w.mean <- numeric(K)
		w.map <- numeric(K)
		for(k in 1:K){
			w.mean[k] <- mean(w.mcmc[,k,1])
			w.map[k] <- w.mcmc[mapIndex,k,1]
		}
		write.table(w.mcmc, file = 'reordered_weights_ecr.txt')
		write.table(w.mean, file = 'weights_estimate_ecr.txt')
		write.table(ls$clusters, file = "singleBestClusterings.txt", quote = FALSE, row.names = FALSE)
		if(  kSelected != length(table(ls$clusters[1,])) ){
			write.table(ls$clusters[c(2,1),], file = "singleBestClusterings.txt", quote = FALSE, row.names = FALSE)
		}
		zMAP <- ls$clusters[lsMethod,]

		mapAllocationsPosteriorProbs <- numeric(n)
		for(i in 1:n){
			mapAllocationsPosteriorProbs[i] <- length(which(allocationsECR[,i] == zMAP[i]))
		}
		mapAllocationsPosteriorProbs <- mapAllocationsPosteriorProbs/dim(allocationsECR)[1]
		write.table(mapAllocationsPosteriorProbs, file = "classificationProbabilities.txt")

		aliveClusters <- as.numeric(names(table(zMAP)))
		if(q > 0){
#			covmat <- array(data = 0, dim = c( length(aliveClusters), p, p ))
#			rownames(covmat) <- as.character(aliveClusters)
#			if(sameSigma == TRUE){
#				for(iter in 1:m){
#					for(k in aliveClusters){
#						lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
#						covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
#						diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(1/sigmaINV[iter, ])
#					}
#				}
#			}else{
#				for(iter in 1:m){
#					for(k in aliveClusters){
#						lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
#						covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
#						diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(Sigma.mcmc[iter,k,])
#					}
#				}
#			}
#			covmat <- covmat/m
######################################################################################################################################################################3
#			computing covmat but with median of Sigma instead ergodic mean
			covmat <- array(data = 0, dim = c( length(aliveClusters), p, p ))
			rownames(covmat) <- as.character(aliveClusters)
			if(sameSigma == TRUE){
				for(k in aliveClusters){
					for(iter in 1:m){
						lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
						covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
					}
					covmat[as.character(k), , ] <- covmat[as.character(k), , ]/m
					diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(apply(1/sigmaINV, 2, median))
				}
			}else{
				for(k in aliveClusters){
					for(iter in 1:m){
						lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
						covmat[as.character(k), , ] <- covmat[as.character(k), , ] + lambda_k %*% t(lambda_k) 
					}
					covmat[as.character(k), , ] <- covmat[as.character(k), , ]/m
				
					diag(covmat[as.character(k), , ]) <- diag(covmat[as.character(k), , ]) + as.numeric(apply( Sigma.mcmc[ ,k, ], 2, median ))
				}
			}
######################################################################################################################################################################3
			for(k in 1:length(aliveClusters)){
				write.table(covmat[k, , ], file = paste0("estimated_cov_cluster_",k,".txt"))
#				cat(paste0('         * write file: `estimated_cov_cluster_',k,'.txt`'),'\n')
				write.table(lambda.map[k, , ], file = paste0('lambda_map_',k,'.txt'))
#				cat(paste0('         * write file: `lambda_map_',k,'.txt`'),'\n')
			}
		}



		if( compute_regularized_expression == TRUE ){
#			cat(paste("-    computing regularized expressions..."),'\n')
			yValues <- read.table("yValues.txt", colClasses = rep('numeric',n*q))
			if(burn > 0){
				yValues <- yValues[-(1:burn), ]
			}
			yValues <- yValues[Kindex, ]
			regularizedExpression <- array(data = 0, dim = c(length(aliveClusters), p))
			regularizedExpression2 <- array(data = 0, dim = c(length(aliveClusters), p, q))			
			rownames(regularizedExpression) <- rownames(regularizedExpression2) <- as.character(aliveClusters)
			for(iter in 1:m){
				for(k in aliveClusters){
					lambda_k <- matrix( lambda.perm.mcmc[iter, k, ], nrow = p, ncol = q, byrow=T)
					yMean <- array(data = 0, dim = c(q,1))
					index_k <- as.numeric(which( allocationsECR[iter, ] == k ))
					for(i in index_k ){
						yMean <- yMean + t(yValues[iter , (1:q - 1)*n + i])
					}
					yMean <- yMean/length(index_k)
					tmp <- t(apply(lambda_k,1, function(y){ y*yMean }))
					regularizedExpression2[as.character(k), , ] <- regularizedExpression2[as.character(k), , ] + tmp
					regularizedExpression[as.character(k), ] <- regularizedExpression[as.character(k), ] + rowSums(tmp)
				}
#				if(iter %% 50 == 0){cat(paste('iter = ', iter),'\n')}
			}
			regularizedExpression <- regularizedExpression/m
			regularizedExpression2 <- regularizedExpression2/m
			for(j in 1:q){
				write.table(
				regularizedExpression2[,,j], 
				file = paste0("estimated_regularized_expression_per_cluster_",j,".txt"))
#				cat(paste0('         * write file: `estimated_regularized_expression_per_cluster_',j,'.txt`'),'\n')

			}
			write.table(regularizedExpression, file = "estimated_regularized_expression_per_cluster.txt")
#			cat(paste0('         * write file: `estimated_regularized_expression_per_cluster.txt`'),'\n')
		}


		cat(paste0('-    Done.'), '\n')
	}
	setwd("../")

}




#new in version 3
# overall main function
fabMix <- function(model = c("UUU", "CUU", "UCU", "CCU", "UCC", "UUC", "CUC", "CCC"), 
			nChains = NULL,	dirPriorAlphas, rawData, outDir, Kmax, mCycles, 
			burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize = TRUE, 
			thinning = 1, zStart, nIterPerCycle, gibbs_z = 1, 
			warm_up_overfitting = 500, warm_up = 5000,  overfittingInitialization=TRUE, 
			progressGraphs = FALSE, gwar = 0.05, rmDir = TRUE, parallelModels = NULL, lowerTriangular=TRUE			
			){

	cat("         ____      __    __  ____     ", "\n")
	cat("        / __/___ _/ /_  /  |/  (_)  __", "\n")
	cat("       / /_/ __ `/ __ \\/ /|_/ / / |/_/", "\n")
	cat("      / __/ /_/ / /_/ / /  / / />  <  ", "\n")
	cat("     /_/  \\__,_/_.___/_/  /_/_/_/|_|  version 5.0", "\n\n")

	model = intersect(model, c("UUU", "CUU", "UCU", "CCU", "UCC", "UUC", "CUC", "CCC"))
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if(is.null(nChains)){
		if( missing(dirPriorAlphas) ){
			nChains <- 8
			dN <- 1
			dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
		}else{
			nChains <- length(dirPriorAlphas)
		}
	}else{
		if( missing(dirPriorAlphas) ){
			if(nChains > 1){
				dN <- 1
				dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
			}else{
				dirPriorAlphas <- 1/Kmax
			}
		}else{
			if( length(dirPriorAlphas) != nChains ){
				stop('dirPriorAlphas should have length equal to nChains.')
			}
		}
	}


	if( nChains > 1 ){
		if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
		if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	}
	if(burnCycles < 1){ burnCycles <- 1;  mCycles <- mCycles + 1}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}

	p <- dim(rawData)[2]
	n <- dim(rawData)[1]
	cat(paste0("-    Data consists of p = ", p, " variables and n = ",n," observations","\n"))
	cat(paste0("-    Prior parameters: g = ", g, ", h = ", h, ", alpha_sigma = ", alpha_sigma, ", beta_sigma = ", beta_sigma,"\n"))
	cat(paste0('-         using Nchains = ', nChains),'\n')
	cat(paste0('-         target posterior distribution corresponds to alpha = ', dirPriorAlphas[1]),'\n')
	if( normalize == TRUE ){
		cat('-    The sampler uses standardized data.','\n')
	}
	if( normalize == FALSE ){
		cat('-    The sampler uses raw data (NOT GOOD PRACTICE).','\n')
	}




	ledermannBound <- ( 2*p + 1 - sqrt(8*p + 1) ) / 2
	if( max(q) >  ledermannBound){
		cat(paste0('-    WARNING: I will consider only the number of factors that do not exceed the ledermann bound.'),'\n')
		cat(paste0('-         so: q <= ', ledermannBound),'\n')
		q <- q[q <= ledermannBound]
	}



	if( is.numeric(parallelModels) ){
		fpr <- fabMix_parallelModels(model = model, 
			nChains = nChains, dirPriorAlphas = dirPriorAlphas, rawData = rawData, outDir = outDir, Kmax = Kmax, mCycles = mCycles, 
			burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, q = q, normalize = normalize, 
			thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, gibbs_z = gibbs_z, 
			warm_up_overfitting = warm_up_overfitting, warm_up = warm_up,  overfittingInitialization=overfittingInitialization, 
			progressGraphs = FALSE, gwar = 0.05, rmDir = rmDir, parallelModels = parallelModels, lowerTriangular=lowerTriangular)
			return(fpr)	#exit
	}


	# define output objects
	bic <- array(data = NA, dim = c(length(q), length(model)))
	colnames(bic) <- model
	rownames(bic) <- q
	nClusters <- array(data = NA, dim = c(length(q), length(model)))
	colnames(nClusters) <- model
	rownames(nClusters) <- q
	Kmap_prob <- array(data = NA, dim = c(length(q), length(model)))
	swap_rate <- array(data = NA, dim = c(length(q), length(model)))
	colnames(Kmap_prob) <- colnames(swap_rate) <- model
	rownames(Kmap_prob) <- rownames(swap_rate) <- q

	if(dir.exists(outDir)){
		stop(paste0('Directory `',outDir,'` exists, please supply another name.'))
	}
	dir.create(outDir)
	setwd(outDir)
	rememberOutDir <- outDir
	check_if_at_least_one_model <- 0
	if("UUU" %in% model){
		for(nFactors in q){
			myDir <- paste0("UUU_",nFactors)
			fabMix_UxU(sameSigma = FALSE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = FALSE, sameLambda = FALSE, isotropic = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "UUU"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "UUU"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "UUU"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "UUU"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if("UCU" %in% model){
		for(nFactors in q){
			myDir <- paste0("UCU_",nFactors)
			fabMix_UxU(sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = TRUE, sameLambda = FALSE, isotropic = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "UCU"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "UCU"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "UCU"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "UCU"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if("CUU" %in% model){
		for(nFactors in q){
			myDir <- paste0("CUU_",nFactors)
			fabMix_CxU(sameSigma = FALSE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = FALSE, sameLambda = TRUE, isotropic = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "CUU"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "CUU"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "CUU"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "CUU"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if("CCU" %in% model){
		for(nFactors in q){
			myDir <- paste0("CCU_",nFactors)
			fabMix_CxU(sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = TRUE, sameLambda = TRUE, isotropic = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "CCU"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "CCU"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "CCU"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "CCU"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if("CUC" %in% model){
		for(nFactors in q){
			myDir <- paste0("CUC_",nFactors)
			fabMix_CxC(sameSigma = FALSE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = FALSE, sameLambda = TRUE, isotropic = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "CUC"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "CUC"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "CUC"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "CUC"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if("CCC" %in% model){
		for(nFactors in q){
			myDir <- paste0("CCC_",nFactors)
			fabMix_CxC(sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = TRUE, sameLambda = TRUE, isotropic = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "CCC"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "CCC"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "CCC"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "CCC"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if("UUC" %in% model){
		for(nFactors in q){
			myDir <- paste0("UUC_",nFactors)
			fabMix_UxC(sameSigma = FALSE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = FALSE, sameLambda = FALSE, isotropic = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = FALSE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "UUC"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "UUC"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "UUC"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "UUC"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if("UCC" %in% model){
		for(nFactors in q){
			myDir <- paste0("UCC_",nFactors)
			fabMix_UxC(sameSigma = TRUE, dirPriorAlphas = dirPriorAlphas, rawData = rawData, 
				outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
				burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, 
				beta_sigma = beta_sigma, q = nFactors, normalize = normalize, 
				thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, 
				gibbs_z = gibbs_z, warm_up_overfitting = warm_up_overfitting, warm_up = warm_up, 
				overfittingInitialization=overfittingInitialization, progressGraphs = progressGraphs, gwar = gwar, lowerTriangular=lowerTriangular)
			if(progressGraphs==TRUE){dev.off()}
			getStuffForDIC(sameSigma = TRUE, sameLambda = FALSE, isotropic = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, Km = Kmax)
#			dealWithLabelSwitching(sameSigma = TRUE, x_data = rawData, outputFolder = myDir, q = nFactors, compute_regularized_expression = FALSE, Km = Kmax)
			bic[as.character(nFactors), "UCC"] <- read.table(paste0(myDir,"/informationCriteria_map_model.txt"))[4, ]
			logl <- read.table(paste0(myDir,"/kValues.txt"), header=T)
			nClusters[as.character(nFactors), "UCC"] <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
			r.freq <- 100*as.numeric(table(logl[,1])/length(logl[,1])) 
			Kmap_prob[as.character(nFactors), "UCC"] <- max(r.freq)
			sr <- read.table(paste0(myDir,"/swapRate.txt"))
			swap_rate[as.character(nFactors), "UCC"] <- as.numeric(tail(sr, 1))
			check_if_at_least_one_model <- 	check_if_at_least_one_model + 1
		}
	}

	if(check_if_at_least_one_model < 1){ 
		stop("Please specify at least one valid model.")
	}

	#selected model
	model_selected <- colnames(bic)[which.min(apply(bic,2,min))]
	q_selected <- rownames(bic)[which.min(apply(bic,1,min))]
	outDir <- paste0(model_selected,"_",q_selected)
	cat(paste0('-    The label.switching package says hello.'), '\n')
	if( strsplit(model_selected,split="")[[1]][2] == "C" ){
		dealWithLabelSwitching(sameSigma = TRUE, x_data = rawData, outputFolder = outDir, q = as.numeric(q_selected), compute_regularized_expression = TRUE, Km = Kmax)
	}else{
		dealWithLabelSwitching(sameSigma = FALSE, x_data = rawData, outputFolder = outDir, q = as.numeric(q_selected), compute_regularized_expression = TRUE, Km = Kmax)
	}
#	if ( nClusters[q_selected, model_selected] > 1){
		z <- as.numeric(read.table(paste0(outDir,"/singleBestClusterings.txt"), header=TRUE)[1,])
		if( nClusters[q_selected, model_selected] != length(table(z)) ){
			nClusters[q_selected, model_selected] <- length(table(z))
			cat(paste0('-    WARNING: the number of classes in the reordered single best clustering does not agree with the raw output of non-empty components.'), '\n')
		}
		classification_probabilities_MAP <- read.table(paste0(outDir,"/classificationProbabilities.txt"))
		covariance_matrix_MAP <- vector("list", length = nClusters[q_selected, model_selected])
#		for (k in 1:nClusters[q_selected, model_selected]){ # this produces error when K_ecr < K, so change to the next line
		for (k in 1:length(names(table(z)))){
			covariance_matrix_MAP[[k]] <- as.matrix(read.table(paste0(outDir,"/estimated_cov_cluster_",k,".txt")))
		}
		names(covariance_matrix_MAP) <- names(table(z))
		mu <- read.table(paste0(outDir,"/mu_estimate_ecr.txt"), header=TRUE)
		mu <- t(mu[as.numeric(names(table(z))),])
		colnames(mu) <- names(table(z))
		w <- read.table(paste0(outDir,"/weights_estimate_ecr.txt"))[as.numeric(names(table(z))),]
		w <- w/sum(w)	
		names(w) <- names(table(z))

		MCMC <- vector("list", length = 8)
		names(MCMC) <- c("y", "w", "Lambda","mu","z","Sigma","K_all_chains", "lambda_map")
		lFile <- read.table(paste0(outDir,"/reordered_lambda_ecr.txt"), header = TRUE)
		if( strsplit(model_selected, split = "")[[1]][1] == "U" ){
			MCMC$Lambda <- vector("list", length = nClusters[q_selected, model_selected])
			names(MCMC$Lambda) <- names(table(z))
			for(k in names(table(z))){
				MCMC$Lambda[[k]] <- array(data = NA, dim = c(dim(lFile)[1], p*as.numeric(q_selected)) )
				colnames(MCMC$Lambda[[k]]) <- paste0('V',rep(1:p, each = as.numeric(q_selected)),'_F',rep(1:as.numeric(q_selected), p))
				f <- 0
				for(i in 1:p){
					for(j in 1:as.numeric(q_selected)){
						f <- f + 1
						MCMC$Lambda[[k]][, f] <- lFile[, paste0("k",k,"_i",i,"_j",j)]
					}
				}
				MCMC$Lambda[[k]] <- mcmc(MCMC$Lambda[[k]], 
							start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
							thin = nIterPerCycle)
			}
		}else{
			k <- names(table(z))[1]
			MCMC$Lambda <- array(data = NA, dim = c(dim(lFile)[1], p*as.numeric(q_selected)) )
			colnames(MCMC$Lambda) <- paste0('V',rep(1:p, each = as.numeric(q_selected)),'_F',rep(1:as.numeric(q_selected), p))
			f <- 0
			for(i in 1:p){
				for(j in 1:as.numeric(q_selected)){
					f <- f + 1
					MCMC$Lambda[, f] <- lFile[, paste0("k",k,"_i",i,"_j",j)]
				}
			}
			MCMC$Lambda <- mcmc(MCMC$Lambda, 
						start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
						thin = nIterPerCycle)
		}
		muValues <- read.table(paste0(outDir,"/reordered_mu_ecr.txt"))
		MCMC$mu <- vector("list", length = nClusters[q_selected, model_selected])
		names(MCMC$mu) <- names(table(z))
		for(k in names(table(z))){
			MCMC$mu[[k]] <- as.matrix(muValues[,as.numeric(k) + Kmax*((1:p)-1)])
			colnames(MCMC$mu[[k]]) <- paste0("V",1:p)
			MCMC$mu[[k]] <- mcmc(MCMC$mu[[k]], 
						start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
						thin = nIterPerCycle)
		}
		MCMC$w <- read.table(paste0(outDir,"/reordered_weights_ecr.txt"))[,as.numeric(names(table(z)))]
		if( nClusters[q_selected, model_selected] == 1 ){
			MCMC$w <- as.matrix(MCMC$w)
		}
		colnames(MCMC$w) <- names(table(z))
		MCMC$w <- mcmc(MCMC$w,
				start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
				thin = nIterPerCycle)
		MCMC$y <- read.table(paste0(outDir,"/yValues.txt"), colClasses = rep('numeric',n*as.numeric(q_selected)))
		colnames(MCMC$y) <- paste0(rep( paste0('OBS',1:n), as.numeric(q_selected)), rep( paste0('_F',1:as.numeric(q_selected)), each = n))
		MCMC$y <- mcmc(MCMC$y,
				start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
				thin = nIterPerCycle)
		MCMC$z <- read.table(paste0(outDir,"/reordered_allocations_ecr.txt"))
		if( strsplit(model_selected, split = "")[[1]][2] == "U" ){
			sMATRIX <- read.table(paste0(outDir,'/reordered_sigma_ecr.txt'))
			MCMC$Sigma <- vector("list", length = nClusters[q_selected, model_selected])
			names(MCMC$Sigma) <- names(table(z))
			if( strsplit(model_selected, split = "")[[1]][3] == "U" ){
				for(k in names(table(z))){
					myColumns <- as.numeric(Kmax)*(1:p) - (Kmax - as.numeric(k))
					MCMC$Sigma[[k]] <- as.matrix(sMATRIX[,myColumns])
					MCMC$Sigma[[k]] <- mcmc(MCMC$Sigma[[k]], 
								start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle,  
								thin = nIterPerCycle)
					colnames(MCMC$Sigma[[k]]) <- paste0('V',1:p)
				}
			}else{
				for(k in names(table(z))){
					MCMC$Sigma[[k]] <- as.matrix(sMATRIX[,as.numeric(k)])
					MCMC$Sigma[[k]] <- mcmc(MCMC$Sigma[[k]], 
								start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
								thin = nIterPerCycle)
				}
			}		
		}else{
#			MCMC$Sigma <- 1/read.table(paste0(outDir, '/sigmainvValues.txt'))
			MCMC$Sigma <- read.table(paste0(outDir,'/reordered_sigma_ecr.txt'))
			if( strsplit(model_selected, split = "")[[1]][3] == "U" ){
				MCMC$Sigma <- MCMC$Sigma
				colnames(MCMC$Sigma) <- paste0('V',1:p)
			}else{
				MCMC$Sigma <- MCMC$Sigma[,1]
			}		
			MCMC$Sigma <- mcmc(MCMC$Sigma, 
						start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
						thin = nIterPerCycle)

		}
		MCMC$K_all_chains <- read.table(paste0(outDir,"/kValues.txt"), header=TRUE)
		MCMC$K_all_chains <- mcmc(MCMC$K_all_chains, 
					start = warm_up_overfitting + warm_up + burnCycles*nIterPerCycle, 
					thin = nIterPerCycle)
		MCMC$lambda_map <- read.table(paste0(outDir, '/lambda_map.txt'))
		rr <- vector("list", length = as.numeric(q_selected))
		names(rr) <- paste0('F',1:as.numeric(q_selected))
		if ( nClusters[q_selected, model_selected] > 1){
			for(j in 1:as.numeric(q_selected)){
				rr[[j]] <- read.table(paste0(outDir,"/estimated_regularized_expression_per_cluster_",j,".txt"))
				rownames(rr[[j]]) <-  names(table(z))
			}
		}
#	}


	setwd("../")
	if(rmDir == TRUE){
		cat(paste0('-    Cleaning: deleting directory `', rememberOutDir, '` ...'))
		unlink(rememberOutDir, recursive=TRUE)
		cat(' done.', '\n')
	}


	cat(paste0("\n","Given the specified range of models, factors, maximum number of clusters and Prior parameters,","\n", "the best model corresponds to the ", model_selected, " parameterization with q = ", q_selected, " factors and K = ",nClusters[q_selected, model_selected]," clusters. ","\n","The BIC for this model equals ", round(min(bic),3), "."),"\n")
	best_model <- data.frame(parameterization = model_selected, num_Clusters = nClusters[q_selected, model_selected], num_Factors = as.numeric(q_selected))

	results <- vector("list", length = 13)
	results[[1]] <- bic
	results[[3]] <- nClusters
	results[[8]] <- best_model
	results[[10]] <- rawData
	results[[12]] <- Kmap_prob
	results[[13]] <- swap_rate
#	if ( nClusters[q_selected, model_selected] > 1){
		results[[4]] <- classification_probabilities_MAP
		results[[2]] <- z
		results[[5]] <- covariance_matrix_MAP
		results[[6]] <- mu
		results[[7]] <- w
		results[[9]] <- MCMC
		results[[11]] <- rr
#	}else{
#		results[[2]] <- rep(1,n)
#		results[[4]] <- rep(1,n)
#		results[[7]] <- 1
#	}

	names(results) <- c(	"bic", 
				"class", 
				"n_Clusters_per_model", 
				"posterior_probability", 
				"covariance_matrix", 
				"mu", 
				"weights", 
				"selected_model", 
				"mcmc",
				"data",
				"regularizedExpression",
				"Kmap_prob",
				"swap_rate"
			)
	class(results) <- c('list', 'fabMix.object')
	return(results)
}





fabMix_parallelModels <- function(model = c("UUU", "CUU", "UCU", "CCU", "UCC", "UUC", "CUC", "CCC"), 
			nChains = NULL,dirPriorAlphas, rawData, outDir, Kmax, mCycles, 
			burnCycles, g, h, alpha_sigma, beta_sigma, q, normalize =TRUE, 
			thinning=1, zStart, nIterPerCycle, gibbs_z = 1, 
			warm_up_overfitting, warm_up,  overfittingInitialization = TRUE, 
			progressGraphs = FALSE, gwar = 0.05, rmDir = TRUE, parallelModels, lowerTriangular=TRUE){

#	copy default values
	model = intersect(model, c("UUU", "CUU", "UCU", "CCU", "UCC", "UUC", "CUC", "CCC"))
	if(missing(Kmax)){Kmax <- 20}
	if(missing(nIterPerCycle)){nIterPerCycle = 10}
	if(missing(zStart)){zStart = FALSE}
	if(is.null(nChains)){
		if( missing(dirPriorAlphas) ){
			nChains <- 8
			dN <- 1
			dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
		}else{
			nChains <- length(dirPriorAlphas)
		}
	}else{
		if( missing(dirPriorAlphas) ){
			if(nChains > 1){
				dN <- 1
				dirPriorAlphas <- c(1, 1 + dN*(2:nChains - 1))/Kmax
			}else{
				dirPriorAlphas <- 1/Kmax
			}
		}else{
			if( length(dirPriorAlphas) != nChains ){
				stop('dirPriorAlphas should have length equal to nChains.')
			}
		}
	}


	if( nChains > 1 ){
		if( range(diff(order(dirPriorAlphas)))[1] != 1){stop('dirPriorAlphas should be in increasing order.')}
		if( range(diff(order(dirPriorAlphas)))[2] != 1){stop('dirPriorAlphas should be in increasing order.')}
	}
	if(mCycles < burnCycles + 1){ stop('`burnCycles` should be less than `mCycles`.') } 
	if(missing(g)){g <- 0.5}
	if(missing(h)){h <- 0.5}
	if(missing(alpha_sigma)){alpha_sigma <- 0.5}
	if(missing(beta_sigma)){beta_sigma <- 0.5}
###################################################
	dir.create(outDir)
	setwd(outDir)
	nModels <- length(model)
	if(parallelModels > nModels){parallelModels <- nModels}
	subRuns <- split(model, factor(sort(rank(model)%%parallelModels)))
	registerDoParallel(parallelModels)
	cat(paste0('-    Fitting ',nModels,' models using ',parallelModels,' parallel threads.'), '\n')
	for(i in 1:parallelModels){
		cat(paste0('-        [NOTE] Screen output from parallel run ',i,' is redirected to `',outDir ,'/print_parallelRun_',i,'.txt`.'), '\n')
	}
	i <-0
	gwd <- getwd()	# for windows compatibility (otherwise, the foreach loops changes the working directory)
	myList <- foreach(i=1:parallelModels, .export=ls(globalenv())) %dopar% {
		setwd(gwd)
		sink(paste0("print_parallelRun_",i,".txt"))
		myDir <- paste0('parallelRun_',i)
		fam <- fabMix(model = subRuns[[i]], 
			nChains = nChains, dirPriorAlphas = dirPriorAlphas, rawData = rawData, outDir = myDir, Kmax = Kmax, mCycles = mCycles, 
			burnCycles = burnCycles, g = g, h = h, alpha_sigma = alpha_sigma, beta_sigma = beta_sigma, q = q, normalize = normalize, 
			thinning = thinning, zStart = zStart, nIterPerCycle = nIterPerCycle, gibbs_z = gibbs_z, 
			warm_up_overfitting = warm_up_overfitting, warm_up = warm_up,  overfittingInitialization=overfittingInitialization, 
			progressGraphs = FALSE, gwar = gwar, rmDir = rmDir, lowerTriangular=lowerTriangular)
		return(fam)
		sink()		
	}
	stopImplicitCluster()
	setwd('../')

	# merging all results together and return a fabMix object
	mergedResult <- vector('list', length = 13)	
	names(mergedResult) <- c("bic", 
				"class", 
				"n_Clusters_per_model", 
				"posterior_probability", 
				"covariance_matrix", 
				"mu", 
				"weights", 
				"selected_model", 
				"mcmc",
				"data",
				"regularizedExpression",
				"Kmap_prob",
				"swap_rate"
			)
	class(mergedResult) <- c('list', 'fabMix.object')

	mergedResult$bic <- matrix(unlist(lapply(myList, function(y)y$bic)), nrow = length(q), ncol = length(model), dimnames = list(q, model))
	mergedResult$Kmap_prob <- matrix(unlist(lapply(myList, function(y)y$Kmap_prob)), nrow = length(q), ncol = length(model), dimnames = list(q, model))
	mergedResult$swap_rate <- matrix(unlist(lapply(myList, function(y)y$swap_rate)), nrow = length(q), ncol = length(model), dimnames = list(q, model))
	mergedResult$data <- myList[[1]]$data
	mergedResult$n_Clusters_per_model <- matrix(unlist(lapply(myList, function(y)y$n_Clusters_per_model)), nrow = length(q), ncol = length(model), dimnames = list(q, model)) 
	model_selected <- model[which.min(apply(mergedResult$bic,2,min))]
	q_selected <- q[which.min(apply(mergedResult$bic,1,min))]
#	model_pointer <- which(model == model_selected)
	model_pointer <- as.numeric(which(unlist(lapply(subRuns, function(y)model_selected %in% y)) == TRUE))
	mergedResult$class <- myList[[model_pointer]][['class']]
	mergedResult$posterior_probability <- myList[[model_pointer]][['posterior_probability']]
	mergedResult$covariance_matrix <- myList[[model_pointer]][['covariance_matrix']]
	mergedResult$mu <- myList[[model_pointer]][['mu']]
	mergedResult$weights <- myList[[model_pointer]][['weights']]
	mergedResult$selected_model <- myList[[model_pointer]][['selected_model']]
	mergedResult$mcmc <- myList[[model_pointer]][['mcmc']]
	mergedResult$regularizedExpression <- myList[[model_pointer]][['regularizedExpression']]

	cat(paste0("\n","Given the specified range of models, factors, maximum number of clusters and Prior parameters,","\n", "the best model corresponds to the ", model_selected, " parameterization with q = ", q_selected, " factors and K = ",mergedResult$selected_model$num_Clusters," clusters. ","\n","The BIC for this model equals ", round(min(mergedResult$bic[as.character(q_selected),model_selected]),3), "."),"\n")

	if(rmDir == TRUE){
		unlink(outDir, recursive=TRUE)
	}


	return(mergedResult)
}



simData <- function(sameSigma = TRUE, sameLambda = FALSE, p, q, K.true, n, loading_means, loading_sd, sINV_values){
	if(missing(p)){p = 40}
	if(missing(q)){q = 4}
	if(missing(K.true)){K.true = 6}
	if(missing(n)){n = 200}
	if(missing(sINV_values)){
		if(sameSigma){ 
			sINV_values = rgamma(p, shape = 1, rate = 1) 
		}else{
			sINV_values = matrix(rgamma(K.true*p, shape = 1, rate = 1), nrow = K.true, ncol = p )
		}
	}
	if( missing(loading_means) ){ loading_means <- c(-30,-20,-10,10, 20, 30) }
	if( missing(loading_sd) ){ loading_sd <- rep(2, length(loading_means)) }
	if ( length(loading_means) !=  length(loading_sd) ){
		stop("`loading_means` and `loading_sd` should have same length.")
	}
	cat(paste0("Simulation parameters:"),'\n')
	if(q >= p){stop("q should not be greater than p")}
	cat(paste0("   n = ", n),'\n')
	cat(paste0("   p = ", p),'\n')
	cat(paste0("   q = ", q),'\n')
	cat(paste0("   K = ", K.true),'\n')
	w.true <- myDirichlet(rep(10,K.true))
	z.true <- sample(K.true,n,replace = TRUE, prob = w.true)
	if(sameLambda == FALSE){
		Lambda.true <- array(data = rnorm(K.true*p*q, mean = 0, sd = 1), dim = c(K.true,p,q))
	}else{
		Lambda.true <- array(data = rnorm(K.true*p*q, mean = 0, sd = 1), dim = c(K.true,p,q))
		for(k in 2:K.true){
			Lambda.true[k,,] <- Lambda.true[1,,]
		}
	}
	mu.true <- array(data = 0, dim = c(K.true,p))
	for(k in 1:K.true){
		u <- runif(1)
		subROW <- floor(p/q)
		for(j in 1:q){
			meanCOL <- rep(0,p)
			pickAnIndex <- sample(length(loading_means), 1)
			meanCOL[ (j-1)*subROW + 1:subROW] <- loading_means[pickAnIndex]
			if(sameLambda == FALSE){
				Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = loading_sd[pickAnIndex] )
				if(j > 1)  {
					Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)
				}
			}else{
				Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = loading_sd[pickAnIndex] )
				if(j > 1)  {
					Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)
				}

				if(k > 1){
					Lambda.true[k, , j] <- Lambda.true[1, , j]
				}
			}
		}
		u <- runif(1)
		if(u < 1/3){
			mu.true[k, ] <- 20*sin(seq(0,k*pi, length = p))
		}else{
			if(u < 2/3){
				mu.true[k, ] <- 20*cos(seq(0,k*pi, length = p))
			}else{
				mu.true[k, ] <- 40*(sin(seq(0,k*pi, length = p)))^2 - 40*(cos(seq(0,k*pi, length = p)))^2
			}
		}
	}


	if(sameSigma == TRUE){
		SigmaINV.true <- array(data = 0, dim = c(p,p))
		diag(SigmaINV.true) <- sINV_values
		Sigma.true <- SigmaINV.true
		diag(Sigma.true) <- 1/diag(SigmaINV.true)

	}else{
		SigmaINV.true <- array(data = 0, dim = c(K.true, p,p))
		Sigma.true <- SigmaINV.true
		for(k in 1:K.true){
			diag(SigmaINV.true[k,,]) <- sINV_values[k,]
			diag(Sigma.true[k,,]) <- 1/diag(SigmaINV.true[k,,])
		}
	}
	y.true <- array(data = 0, dim = c(n,q))
	x_data <- array(data = 0, dim = c(n,p))
	ly <- q
	for(i in 1:n){
		y.true[i,] <- rnorm(ly,mean = 0,sd = 1)
		j <- z.true[i]
		if(q == 1){
			x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% array(y.true[i, ], dim = c(q,1))
		}else{
			x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% y.true[i, ]
		}
		if(sameSigma){
			x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true)
		}else{
			x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true[j,,])
		}
	}
	matplot(t(x_data), type = "l", col = z.true, lty = 1)
	legend("bottomleft", paste0("cluster ",1:K.true, ": ",as.character(as.numeric(table(z.true)))), col = 1:K.true, lty = 1)
	results <- vector('list', length = 7)
	results[[1]] <- x_data
	results[[2]] <- z.true
	results[[3]] <- Lambda.true
	results[[4]] <- mu.true 
	results[[5]] <- Sigma.true
	results[[6]] <- y.true
	results[[7]] <- w.true
	names(results) <- c("data", "class", "factorLoadings", "means", "variance","factors","weights")
	return(results)
}




simData2 <- function(sameSigma = TRUE, p, q, K.true, n, loading_means, loading_sd, sINV_values){
	if(missing(p)){p = 40}
	if(missing(q)){q = 4}
	if(missing(K.true)){K.true = 6}
	if(missing(n)){n = 200}
	if(missing(sINV_values)){
		if(sameSigma){ 
			sINV_values = rgamma(p, shape = 1, rate = 1) 
		}else{
			sINV_values = matrix(rgamma(K.true*p, shape = 1, rate = 1), nrow = K.true, ncol = p )
		}
	}
	if( missing(loading_means) ){ loading_means <- c(-30,-20,-10,10, 20, 30) }
	if( missing(loading_sd) ){ loading_sd <- rep(2, length(loading_means)) }

	if ( length(loading_means) !=  length(loading_sd) ){
		stop("`loading_means` and `loading_sd` should have same length.")
	}
	cat(paste0("Simulation parameters:"),'\n')
	if(q >= p){stop("q should not be greater than p")}
	cat(paste0("   n = ", n),'\n')
	cat(paste0("   p = ", p),'\n')
	cat(paste0("   q = ", q),'\n')
	cat(paste0("   K = ", K.true),'\n')
	w.true <- myDirichlet(rep(10,K.true))
	Lambda.true <- array(data = rnorm(K.true*p*q, mean = 0, sd = 1), dim = c(K.true,p,q))
	mu.true <- array(data = 0, dim = c(K.true,p))
	for(k in 1:K.true){
		if(runif(1) < 0.5){
			pickAnIndex <- sample(length(loading_means), 1)
			u <- runif(1)
			subROW <- floor(p/q)
			mySeq <- seq(from = 1, to = loading_means[pickAnIndex], length = p)
			mySign2 = 1
			for(j in 1:q){
				if(j%%2 == 0) { myIndex = p:1}else{myIndex = 1:p}
				meanCOL <- mySeq
				mySign = 1*mySign2
				for(fa in 1:q){
					meanCOL[ (fa-1)*subROW + 1:subROW] <- mySign*meanCOL[(fa-1)*subROW + 1:subROW]
					mySign = -1*mySign
				}
				Lambda.true[k, , j] <- rnorm(p, mean = meanCOL[myIndex], sd = loading_sd[pickAnIndex] )
				mySign2 = -1*mySign2
			}
		}else{
			u <- runif(1)
			subROW <- floor(p/q)
			for (j in 1:q) {
			    meanCOL <- rep(0, p)
			    pickAnIndex <- sample(length(loading_means), 1)
			    meanCOL[(j - 1) * subROW + 1:subROW] <- loading_means[pickAnIndex]
				Lambda.true[k, , j] <- rnorm(p, mean = meanCOL,  sd = loading_sd[pickAnIndex])
				if (j > 1) {
					Lambda.true[k, 1:(j - 1), j] <- rep(0, j - 1)
				}
			}

		}
		u <- runif(1)
		if(u < 1/3){
			mu.true[k, ] <- 20*sin(seq(0,k*pi, length = p))
		}else{
			if(u < 2/3){
				mu.true[k, ] <- 20*cos(seq(0,k*pi, length = p))
			}else{
				mu.true[k, ] <- 40*(sin(seq(0,k*pi, length = p)))^2 - 40*(cos(seq(0,k*pi, length = p)))^2
			}
		}
	}


	if(sameSigma == TRUE){
		SigmaINV.true <- array(data = 0, dim = c(p,p))
		diag(SigmaINV.true) <- sINV_values
		Sigma.true <- SigmaINV.true
		diag(Sigma.true) <- 1/diag(SigmaINV.true)

	}else{
		SigmaINV.true <- array(data = 0, dim = c(K.true, p,p))
		Sigma.true <- SigmaINV.true
		for(k in 1:K.true){
			diag(SigmaINV.true[k,,]) <- sINV_values[k,]
			diag(Sigma.true[k,,]) <- 1/diag(SigmaINV.true[k,,])
		}
	}
	z.true <- sample(K.true,n,replace = TRUE, prob = w.true)
	y.true <- array(data = 0, dim = c(n,q))
	x_data <- array(data = 0, dim = c(n,p))
	ly <- q
	for(i in 1:n){
		y.true[i,] <- rnorm(ly,mean = 0,sd = 1)
		j <- z.true[i]
		if(q == 1){
			x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% array(y.true[i, ], dim = c(q,1))
		}else{
			x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% y.true[i, ]
		}
		if(sameSigma){
			x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true)
		}else{
			x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true[j,,])
		}
	}
	matplot(t(x_data), type = "l", col = z.true, lty = 1)
	legend("bottomleft", paste0("cluster ",1:K.true, ": ",as.character(as.numeric(table(z.true)))), col = 1:K.true, lty = 1)
	results <- vector('list', length = 7)
	results[[1]] <- x_data
	results[[2]] <- z.true
	results[[3]] <- Lambda.true
	results[[4]] <- mu.true 
	results[[5]] <- Sigma.true
	results[[6]] <- y.true
	results[[7]] <- w.true
	names(results) <- c("data", "class", "factorLoadings", "means", "variance","factors","weights")
	return(results)
}





#' @export
print.fabMix.object <- function(x, ...){
       if( 'fabMix.object' %in% class(x) ){
                cat("\n")
                cat(paste0("* Run information:"),"\n")
                cat(paste0("      Number of fitted models: (", dim(x$bic)[1]," factor levels) x (",dim(x$bic)[2]," parameterizations) = ",prod(dim(x$bic))," models.","\n"))
                cat(paste0("      Selected model: ", as.character(x$selected_model$parameterization)," model with K = ", x$selected_model$num_Clusters, " clusters and q = ", x$selected_model$num_Factors ," factors.","\n"))
                cat("\n")
                cat(paste0("* Maximum A Posteriori (MAP) number of ``alive'' clusters and selected number of factors (BIC) per model:"),"\n")
		if(dim(x$bic)[1] > 1){
			l1 <- as.numeric(as.numeric(rownames(x$bic)[apply(x$bic, 2, order)[1,]]))
			l2 <- as.numeric(diag(x$n_Clusters_per_model[ rownames(x$bic)[apply(x$bic, 2, order)[1,]],]))
			l3 <- round(as.numeric(diag(x$Kmap_prob[ rownames(x$bic)[apply(x$bic, 2, order)[1,]],]))/100,2)
			l4 <- paste0(round(as.numeric(diag(x$swap_rate[ rownames(x$bic)[apply(x$bic, 2, order)[1,]],])),2), '%')
		}else{
			l1 <- rep(rownames(x$bic), dim(x$bic)[2])
			l2 <- as.numeric(x$n_Clusters_per_model)
			l3 <- as.numeric(round(x$Kmap_prob/100, 2))
			l4 <- paste0(as.numeric(round(x$swap_rate, 2)), '%')
		}
		model.info <- data.frame(model = as.factor(colnames(x$bic)), 
					K_MAP = l2, 
					K_MAP_prob =  l3,
					q = l1, 
					BIC_q = as.numeric(as.numeric(round(apply(x$bic, 2, min),1))),
					chain_swap_rate = l4
				)		
		print(model.info)
		cat('\n')
                cat(paste0("* Estimated number of observations per cluster (selected model):"),'\n')
                print(table(x$class, dnn = c('label')))
		cat('\n')

        }else{
                cat(paste("    The input is not in class `fabMix.object`"),'\n')
        }
}


#' @export
summary.fabMix.object <- function(object, quantile_probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...){
	if( 'fabMix.object' %in% class(object) ){
		results <- vector("list", length = 3)
		names(results) <- c("alive_cluster_labels", "posterior_means", "quantiles")
		results$posterior_means <- vector("list", length = 3)
		names(results$posterior_means) <- c("weight", "mean", "covariance_matrix")
		cat(paste0("* ``Alive'' cluster labels:"),'\n')
		print(names(table(object$class)))
		results$alive_cluster_labels <- names(table(object$class))
		cat('\n')		
		m <- round(object$weights, 2)
		cat(paste0('* Posterior mean of the mixing proportions:'),'\n')
		print(m)
		results$posterior_means$weight <- m
		cat('\n')
		m <- round(object$mu, 2)
		names(dimnames(m)) <- list('Variable', 'Cluster label')
		cat(paste0('* Posterior mean of the marginal means per cluster:'),'\n')
		print(m)
		results$posterior_means$mean <- m
		cat('\n')
		cat(paste0('* Posterior mean of the covariance matrix per cluster:'),'\n')		
		for (k in names(table(object$class)) ){
			cat('\n')
			cat(paste0("   Covariance matrix for cluster ``",k,"'':"),'\n')
			m <- round(object$covariance_matrix[[k]], 2)
			rownames(m) <- colnames(m)
			print(m)
		}
		results$posterior_means$covariance_matrix <- object$covariance_matrix
		cat('\n')
		
		if(is.null(quantile_probs)==FALSE){
			cat(paste0('Quantiles for each parameter:'),'\n')
			rp <- range(quantile_probs)
			p <- dim(object$mu)[1]
			if(rp[1] < 0){stop("`quantile_probs` should be between 0 and 1.")}
			if(rp[2] > 1){stop("`quantile_probs` should be between 0 and 1.")}
			K <- length(table(object$class))
			class_names <- names(table(object$class))
			nVariables <- length(object$mu) + length(object$weights) + K*(p^2 - p*(p-1)/2)
			quants <- matrix(nrow = nVariables, ncol = length(quantile_probs))
			colnames(quants) <- paste0(100*quantile_probs, '%')
			covnames <- character(K*(p^2 - p*(p-1)/2))
			cvm <- CovMat_mcmc_summary(object, quantile_probs = quantile_probs)
			t <- 0
			for(i in 1:p){
				for(j in i:p){
					kay <- 0
					for(k in class_names){
						kay <- kay + 1
						t <- t + 1
						covnames[t] <- paste0('cov_',k,'_V',i,'_V',j)
						quants[K*p + K + t, ] <- unlist(lapply(cvm, function(x){return(c(x[i,j,kay]))}))
					}
				}
			}
			rownames(quants) <- c(	paste0('weight_',class_names),
						paste0(paste0('mean_',class_names),'_V',rep(1:p, each = K)), 
						covnames
					)
			names(dimnames(quants)) <- list('parameter', 'quantile')
			# weights
			if(K > 1){
				w <- as.mcmc(t(apply(object$mcmc$w,1,function(x){x/sum(x)})))
			}else{
				w <- object$mcmc$w/object$mcmc$w
			}
			summcmc <- summary(w, quantiles = quantile_probs)$quantiles
			quants[1:K, ] <- as.matrix(summcmc)
			#means
			for(k in class_names){
				summcmc <- summary(object$mcmc$mu[[k]], quantiles = quantile_probs)$quantiles
				myRows <- paste0('mean_',k,'_V',1:p)
				quants[myRows, ] <- as.matrix(summcmc)
			}
			results$quantiles <- quants
			print(round(quants,3))
		}
		return(results)
	}else{
		cat(paste("    The input is not in class `fabMix.object`"),'\n')
	}
}





#' @export
plot.fabMix.object <- function(x, what, variableSubset, class_mfrow = NULL, sig_correlation = NULL, confidence = 0.95, ...){
        if( 'fabMix.object' %in% class(x) ){
	K <- as.numeric(x$selected_model['num_Clusters'])
	cMeans <- colMeans(x$data)
	sdevs <- sqrt(apply(x$data, 2, var))
	p <- dim(x$data)[2]
	v <- 99
	oldpar <- par(no.readonly = TRUE)
	mclustColors <- c("dodgerblue2","red3","green3","slateblue","darkorange","skyblue1",
			"violetred4","forestgreen","steelblue4","slategrey","brown",
			"black","darkseagreen","darkgoldenrod3","olivedrab", "royalblue", 
			"tomato4","cyan2","springgreen2")

	if(missing(what)){what = "menu"}
	if(missing(variableSubset)){
		variableSubset = 1:p
	}
	if(length(variableSubset) < 2){stop("variableSubset should contain at least 2 variables.")}
	if(min(confidence) < 0.5){stop('confidence should be between 0.5 and 1.')}
	if(max(confidence) > 1){stop('confidence should be between 0.5 and 1.')}
	while(v > 0){
		if(what=="menu"){
			cat(paste0("fabMix plots:"),"\n\n")
			cat(paste0("0: Exit"),"\n")
			cat(paste0("1: BIC"),"\n")
			cat(paste0("2: Classification (matplot view)"),"\n")
			cat(paste0("3: Classification (scatter view)"),"\n")
			cat(paste0("4: Correlation plot"),"\n")
			cat(paste0("5: Factor loadings (MAP estimate)"),"\n")
			cat("\n")
			v <- readline(prompt="Selection:")
		}else{
			v = 0
		}
		if((v == 1)||(what == "BIC")){
			on.exit(par())
			tmp.mar <- par()$mar
			# 1. Plot BIC values
#			par(mfrow=c(1,1),mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
			par(mar = c(4,4,1.5,1))
			myMat <- matrix(1:2, nrow = 1, ncol = 2, byrow=TRUE) 
			layout(myMat, widths = c(6.5,1))			
			myCols <- brewer.pal(8, "Set1")
			matplot(x$bic, type = "n", xlab = "number of factors", ylab = "BIC", cex = 0, col = myCols, lty = 1, xaxt = "n")
			rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgrey")
			matplot(x$bic, type = "b", add = TRUE, xlab = "number of factors", ylab = "BIC", cex = 0, col = myCols, lty = 1, xaxt = "n")
			axis(1, at = 1:dim(x$bic)[1], labels = as.numeric(rownames(x$bic)))
			for (i in 1:dim(x$bic)[2]){
				text(labels = x$n_Clusters_per_model[,i], y = x$bic[,i], x = 1:dim(x$bic)[1], col = myCols[i])
			}
			par(mar = c(0,0,0,0))
			plot(1, type = 'n', ylab = '', xlab = '', axes=FALSE)
			legend("center", legend=colnames(x$bic), col = myCols, lty = 1, title="Model", bty = 'n')
			par(mar = tmp.mar)
#			legend("topright", inset=c(-0.3,0.3), legend=colnames(x$bic), col = myCols, lty = 1, title="Model")
		}
		
		if((v == 2)||(what == "classification_matplot")){
			on.exit(par())
			ZNORM <- qnorm(p = (1 - confidence)/2)
			# 2. Matplot per cluster
			tmp.mar <- par()$mar
			par(mar = c(3,3,3,1))
			if(is.null( class_mfrow ) == FALSE){
				if(prod(class_mfrow) < length(table(x$class))){
					stop('class_mfrow not correct.')
				}
				remaining <- prod(class_mfrow) -  length(table(x$class))
				#par(mfrow = class_mfrow)
				myMat <- matrix(1:prod(class_mfrow), nrow = class_mfrow[1], ncol = class_mfrow[2], byrow=TRUE) 
				myMat <- rbind(myMat, rep(prod(class_mfrow) + 1,class_mfrow[2]))
				layout(myMat, heights = c(rep(9,class_mfrow[1]),1))
				for(k in 1:K){
					#ind <- which(x$class == as.numeric(names(table(x$class)))[k])
					ind <- which(x$class == names(table(x$class))[k])
					# apply the suitable transformation:
					sNew <- array(data = NA, dim = c(p,p))
					for(i in 1:p){
						for(j in 1:p){
							sNew[i,j] <- sdevs[i]*sdevs[j]*x$covariance_matrix[[names(table(x$class))[k]]][i,j]
						}
					}
					matplot(t(x$data[ind,]), col = "gray", type = "l", lty = 1, main = paste0("cluster `", names(table(x$class))[k],"' (n = ", as.numeric(table(x$class))[k] ,")"), xlab = "", ylab = "")
#					points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k], type = "l", col =  mclustColors[k], lwd = 2)
					points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k], type = "l", col =  1, lwd = 2)
					mylty <- 1
					for(z_norm in ZNORM){
						mylty = mylty + 1
#						points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k] + z_norm*sqrt(diag(sNew)), type = "l", col =  mclustColors[k], lty = mylty, lwd = 2)
#						points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k] - z_norm*sqrt(diag(sNew)), type = "l", col =  mclustColors[k], lty = mylty, lwd = 2)
						points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k] + z_norm*sqrt(diag(sNew)), type = "l", col =  1, lty = mylty, lwd = 2)
						points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k] - z_norm*sqrt(diag(sNew)), type = "l", col =  1, lty = mylty, lwd = 2)
					}
#					legend("bottomleft", col = c( mclustColors[k], rep(mclustColors[k], length(ZNORM)), "gray"), c("estimated mean", paste0(round(confidence*100,2),"% HDI"), "assigned data"), lty = c(1,2:mylty,1), lwd = 2)
				}
				if(remaining > 0){
					for (i in 1:remaining){
						plot(1, type = 'n', ylab = '', xlab = '', axes=FALSE)					
					}
				}
				par(mar = c(0,0,0,0))
				plot(1, type = 'n', ylab = '', xlab = '', axes=FALSE)
				legend("center",  col = c(rep(1,length(ZNORM) + 1),'gray'), c("mean", paste0(round(confidence*100,2),"% HDI"), "assigned data"), lty = c(1,2:mylty,1), lwd = 2, ncol = 2 + length(ZNORM), bty = 'n')
				par(mar = tmp.mar)
			}else{
				k = 0
				while(k < K ){
					k = k + 1
					#ind <- which(x$class == as.numeric(names(table(x$class)))[k])
					ind <- which(x$class == names(table(x$class))[k])
					# apply the suitable transformation:
					sNew <- array(data = NA, dim = c(p,p))
					for(i in 1:p){
						for(j in 1:p){
							sNew[i,j] <- sdevs[i]*sdevs[j]*x$covariance_matrix[[names(table(x$class))[k]]][i,j]
						}
					}
					matplot(t(x$data[ind,]), col = "gray", type = "l", lty = 1, main = paste0("cluster `", names(table(x$class))[k],"' (n = ", as.numeric(table(x$class))[k] ,")"), xlab = "", ylab = "")
					points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k], type = "l", col =  mclustColors[k], lwd = 2)
					mylty <- 1
					for(z_norm in ZNORM){
						mylty = mylty + 1
						points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k] + z_norm*sqrt(diag(sNew)), type = "l", col =  mclustColors[k], lty = mylty, lwd = 2)
						points(apply(x$mu, 2, function(y){cMeans + sdevs*y})[,k] - z_norm*sqrt(diag(sNew)), type = "l", col =  mclustColors[k], lty = mylty, lwd = 2)
					}
					legend("bottomleft", col = c( mclustColors[k], rep(mclustColors[k], length(ZNORM)), "gray"), 
						c("estimated mean", paste0(round(confidence*100,2),"% HDI"), "assigned data"), lty = c(1,2:mylty,1), lwd = 2)
					if(k < K){
						readline(prompt=paste0("Press ENTER to see next plot... (",k + 1,"/",K,")"))
					}
				}
			}
		}

		if((v == 3)||(what == "classification_pairs")){
			on.exit(par())
			tt <- vector("list", length = 2)
			tt[[1]] <- x$mu
			variance <- vector("list", length = 2)
			variance[[1]] <- "dada"
			variance[[2]] <- array(data = NA, dim = c(p,p,K))
			for(k in 1:K){
				variance[[2]][,,k] <- x$covariance_matrix[[k]]
			}
			names(variance) <- c("dada", "sigma")
			tt[[2]] <- variance
			names(tt) <- c("mean", "variance")
			# 2. Pairwise scatterplots (hacking coordProj() from mclust).
			dimens <- variableSubset
			d <- length(dimens)
			par(mfrow = c(d, d), mar = rep(c(0.3, 0.3/2), each = 2), oma = c(4, 4, 4, 4))
			for (i in seq(d)) {
			  for (j in seq(d)) {
			    if (i == j) {
			      plot(x$data[, c(j, i)], type = "n", xlab = "", ylab = "", axes = FALSE)
			      text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), labels = colnames(x$data[, dimens])[i], cex = 1.5, adj = 0.5)
			      box()
			    }
			    else {
			      coordProj(data = scale(x$data), what = "classification", 
				parameters = tt, classification = x$class, 
				addEllipses = TRUE, dimens = dimens[c(j, i)], 
				main = FALSE, xaxt = "n", yaxt = "n", ...)
			    }
			    if (i == 1 && (!(j%%2))) 
			      axis(3)
			    if (i == d && (j%%2)) 
			      axis(1)
			    if (j == 1 && (!(i%%2))) 
			      axis(2)
			    if (j == d && (i%%2)) 
			      axis(4)
			  }
			}
		}

		if((v == 4)||(what == "correlation")){
			on.exit(par())
			if(is.null(sig_correlation)){sig_correlation = 10}
			if(sig_correlation < 1){
				cm <- CorMat_mcmc_summary(x, quantile_probs = c(0.025, 0.975))
			}


			# 4. Correlation plot per cluster
			if(is.null( class_mfrow ) == FALSE){
				par(mfrow = class_mfrow)
				for(k in 1:K){
					sNew <- array(data = NA, dim = c(p,p))
					for(i in 1:p){
						for(j in 1:p){
							sNew[i,j] <- sdevs[i]*sdevs[j]*x$covariance_matrix[[names(table(x$class))[k]]][i,j]
						}
					}
					if(sig_correlation < 1){
						corrplot(cov2cor(sNew),method = "ellipse", title = paste0("cluster `", names(table(x$class))[k],"'"),  mar=c(0,0,2,0), 
							p.mat = cm$p_matrix[,,k], sig.level = sig_correlation, insig = "pch")
					}else{
						corrplot(cov2cor(sNew),method = "ellipse", title = paste0("cluster `", names(table(x$class))[k],"'"),  mar=c(0,0,2,0))
					}
				}

			}else{
				par(mfrow = c(1,1))
				k = 0
				while(k < K ){
					k = k + 1
					sNew <- array(data = NA, dim = c(p,p))
					for(i in 1:p){
						for(j in 1:p){
							sNew[i,j] <- sdevs[i]*sdevs[j]*x$covariance_matrix[[names(table(x$class))[k]]][i,j]
						}
					}
					if(sig_correlation < 1){
						corrplot(cov2cor(sNew),method = "ellipse", title = paste0("cluster `", names(table(x$class))[k],"'"),  mar=c(0,0,2,0), 
							p.mat = cm$p_matrix[,,k], sig.level = sig_correlation, insig = "pch")
					}else{
						corrplot(cov2cor(sNew),method = "ellipse", title = paste0("cluster `", names(table(x$class))[k],"'"),  mar=c(0,0,2,0))
					}
					if(k < K){
						readline(prompt=paste0("Press ENTER to see next plot... (",k + 1,"/",K,")"))
					}
				}
			}			
		}

		if((v==666)||(what == "regularized_expression")){
			on.exit(par())
			q <- as.numeric(x$selected_model['num_Factors'])
			par(mfrow = c(q, q), mar = rep(c(0.3, 0.3/2), each = 2), oma = c(4, 4, 4, 4))
			for(j in seq(q)){
				for(i in seq(q)){
					if(i == j){
						par(mgp=c(0,-1.4, 0)) 
						matplot(t(x$regularizedExpression[[i]]), type = "n", col = mclustColors, xlab = "", ylab = "", xaxt = 'n', yaxt='n')
						abline(h = 0, lty = 2, col = "gray")
						matplot(t(x$regularizedExpression[[i]]), type = "l", col = mclustColors, xlab = "", ylab = "", xaxt = 'n', yaxt='n', add= TRUE)
						title(paste0("Factor ",i), line = -1)
						#legend("bottomright", paste0("cluster label ", names(table(x$class))), col = mclustColors[1:K], pch = 16)
						if(q > 1){
							axis(1, cex.axis=.7, font=1, tck=.01)
						}
						par(mgp = c(3, 1, 0)) #reset mpg values
					}else{
						plot(range(x$regularizedExpression[[i]]), range(x$regularizedExpression[[j]]), type = "n", xlab = "", ylab = "", xaxt = 'n', yaxt='n')
						for(k in 1:K){
							abline(h = 0, lty = 2, col = "gray")
							abline(v = 0, lty = 2, col = "gray")
							points(as.numeric(x$regularizedExpression[[i]][k,]),as.numeric(x$regularizedExpression[[j]][k,]), pch = k, col = mclustColors[k])
						}
					}
					if (j == 1 && (!(i%%2))) 
					axis(3)
					if (j == q && (i%%2)) 
					axis(1)
					if (i == 1 && (!(j%%2))) 
					axis(2)
					if (i == q && (j%%2)) 
					axis(4)
				}
			}
		}

		if((v == 5)||(what == "factor_loadings")){
			on.exit(par())
			# 4. Factor loadings per cluster
		 	K <- x$selected_model$num_Clusters
			aliveClusters <-  as.numeric(names(table(x$class)))
		        par(mfrow = c(1, 1))
		        firstLetter <- strsplit(as.character(x$selected_model$parameterization), split = "")[[1]][1]
		        q <- as.numeric(x$selected_model["num_Factors"])
		        k = 0
		        while ((k < K)) {
		          k = k + 1
		          Factor = as.factor(rep(paste0("Factor ", 1:q), 
		            each = p))
		          Variable = as.factor(rep(1:p, q))
		          lambda = as.numeric(x$mcmc$lambda_map[aliveClusters[k], ])
		          df <- data.frame(Factor, Variable, lambda)
		          if (firstLetter == "C") {
		            myTitle <- "Factor loadings (MAP)"
		            mySubTitle <- "Common for all clusters"
		          }
		          else {
		            myTitle <- "Factor loadings (MAP)"
		            mySubTitle <- paste0("Cluster `", aliveClusters[k],"'")
		          }
		          p1 <- ggplot(df, aes(x = as.factor(Factor), 
		            y = lambda, fill = Variable)) + geom_bar(position = position_dodge(), 
		            stat = "identity", colour = "black") + coord_flip() + 
		            labs(y = "loading", x = "factor") + labs(title = myTitle, 
		            subtitle = mySubTitle)
			     print(p1)
		          if ((k < K) & ((firstLetter == "U"))) {
		            readline(prompt = paste0("Press ENTER to see next plot... (", 
		              k + 1, "/", K, ")"))
		          }
		          else {
		            k = K
		          }
		        }
		}
		par(mfrow=c(1,1))

	}

        }else{
                cat(paste("    The input is not in class `fabMix.object`"),'\n')
        }
}




CorMat_mcmc_summary <- function(x, quantile_probs){
	if( 'fabMix.object' %in% class(x) ){
		rp <- range(quantile_probs)
		
		if(rp[1] < 0){stop("`quantile_probs` should be between 0 and 1.")}
		if(rp[2] > 1){stop("`quantile_probs` should be between 0 and 1.")}

		firstLetter <- strsplit(as.character(x$selected_model$parameterization), split = "")[[1]][1]
		secondLetter <- strsplit(as.character(x$selected_model$parameterization), split = "")[[1]][2]
		thirdLetter <- strsplit(as.character(x$selected_model$parameterization), split = "")[[1]][3]
		K <- x$selected_model$num_Clusters
		p <- dim(x$data)[2]
		q <- x$selected_model$num_Factors
		quantileList <- vector("list", length = length(quantile_probs))
		names(quantileList) <- as.character(quantile_probs)
		for(j in 1:length(quantile_probs)){
			quantileList[[j]] <- array(0, dim = c(p, p, K))
		}
		Prob_larger_than_zero <- p_matrix <- array(0, dim = c(p, p, K))
		for(j in 1:length(quantile_probs)){
			for(kay in 1:K){
				diag(quantileList[[j]][,,kay]) <- rep(1,p)
			}
		}
		
		for(k in 1:K){
			if(firstLetter == 'C'){
				l <- x$mcmc$Lambda
			}else{
				l <- x$mcmc$Lambda[[k]]
			}
			m <- dim(l)[1]

			if(secondLetter == "C"){
				if(thirdLetter == "C"){
				# xCC
					sigma = matrix(x$mcmc$Sigma, nrow = m, ncol = p )
				}else{
				# xCU
					sigma = as.matrix(x$mcmc$Sigma)
				}
			}else{
				if(thirdLetter == "C"){
				#xUC
					sigma = matrix(x$mcmc$Sigma[[k]], nrow = m, ncol = p )
				}else{
				#xUU
					sigma = as.matrix(x$mcmc$Sigma[[k]])
				}
			}
			# in all cases sigma is an m x p matrix. 

			diag(Prob_larger_than_zero[,,k]) <- rep(1,p)
			seqQ <- 1:q
			for(i in 1:(p-1)){
				for(j in (i+1):p){
					y <- numeric(m)
					for(iter in 1:m){
						x1 <- l[iter,paste0('V',i,'_F',seqQ )]
						x2 <- l[iter,paste0('V',j,'_F',seqQ )]
						y[iter] <- sum(x1*x2)/sqrt((sum(x1^2) + sigma[iter, i])*(sum(x2^2) + sigma[iter, j]))
						Prob_larger_than_zero[i,j,k] <- Prob_larger_than_zero[i, j, k] + as.numeric(y[iter] > 0)
					}
					Prob_larger_than_zero[i,j, k] <- Prob_larger_than_zero[i,j, k]/m
					Prob_larger_than_zero[j, i, k] <- Prob_larger_than_zero[i, j, k]
					quants <- as.numeric(quantile(y, probs = quantile_probs))
					for(h in 1:length(quantile_probs)){
						quantileList[[h]][i,j,k] <- quantileList[[h]][j,i,k] <- quants[h]
					}
				}
			}
			p_matrix[,,k] <- 1 - 2*abs(0.5 - Prob_larger_than_zero[,,k])
		}
		results <- vector("list", length = 2)
		results[[1]] <- quantileList
		results[[2]] <- p_matrix
		names(results) <- c("quantiles", "p_matrix")
		return(results)
        }else{
                cat(paste("    The input is not in class `fabMix.object`"),'\n')
        }

}


CovMat_mcmc_summary <- function(x, quantile_probs){
	if( 'fabMix.object' %in% class(x) ){
		rp <- range(quantile_probs)
		
		if(rp[1] < 0){stop("`quantile_probs` should be between 0 and 1.")}
		if(rp[2] > 1){stop("`quantile_probs` should be between 0 and 1.")}

		firstLetter <- strsplit(as.character(x$selected_model$parameterization), split = "")[[1]][1]
		secondLetter <- strsplit(as.character(x$selected_model$parameterization), split = "")[[1]][2]
		thirdLetter <- strsplit(as.character(x$selected_model$parameterization), split = "")[[1]][3]
		K <- x$selected_model$num_Clusters
		p <- dim(x$data)[2]
		q <- x$selected_model$num_Factors
		quantileList <- vector("list", length = length(quantile_probs))
		names(quantileList) <- as.character(quantile_probs)
		for(j in 1:length(quantile_probs)){
			quantileList[[j]] <- array(0, dim = c(p, p, K))
		}

		for(k in 1:K){
			if(firstLetter == 'C'){
				l <- x$mcmc$Lambda
			}else{
				l <- x$mcmc$Lambda[[k]]
			}
			m <- dim(l)[1]

			if(secondLetter == "C"){
				if(thirdLetter == "C"){
				# xCC
					sigma = matrix(x$mcmc$Sigma, nrow = m, ncol = p )
				}else{
				# xCU
					sigma = as.matrix(x$mcmc$Sigma)
				}
			}else{
				if(thirdLetter == "C"){
				#xUC
					sigma = matrix(x$mcmc$Sigma[[k]], nrow = m, ncol = p )
				}else{
				#xUU
					sigma = as.matrix(x$mcmc$Sigma[[k]])
				}
			}
			# in all cases sigma is an m x p matrix. 


			seqQ <- 1:q
			for(i in 1:p){
				for(j in i:p){
					y <- numeric(m)
					for(iter in 1:m){
						x1 <- l[iter,paste0('V',i,'_F',seqQ )]
						x2 <- l[iter,paste0('V',j,'_F',seqQ )]
						if(j > i){ 
							y[iter] <- sum(x1*x2) 
						}else{
							y[iter] <- sum(x1*x2) + sigma[iter, i]
						}
					}
					quants <- as.numeric(quantile(y, probs = quantile_probs))
					for(h in 1:length(quantile_probs)){
						quantileList[[h]][i,j,k] <- quantileList[[h]][j,i,k] <- quants[h]
					}
				}
			}
		}
		return(quantileList)
        }else{
                cat(paste("    The input is not in class `fabMix.object`"),'\n')
        }

}

