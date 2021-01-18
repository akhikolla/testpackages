blsm_colors = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


proc_crr=function(Z,Z0){
  #' @title Procrustean corresponding positions
  #' @description Given a set of starting coordinates, the function returns the Procrustean Transform of the initial points that minimizes 
  #' the sum of squared positional difference from a set of reference coordinates. The (Euclidean) distances between a candidate 
  #' configuration and the reference are evaluated by considering the couples of corresponding points. 
  #' 
  #' The reference configuration must be centered at the origin. 
  #' 
  #' @param Z set of initial coordinates to be transformed
  #' @param Z0 set of reference coordinates centered at the origin
  #' 
  #' @return Set of coordinates minimizing the distance between the initial configuration and the reference one
  #' @examples 
  #' # Create configuration and center it at the origin
  #' pos_ref = matrix(runif(20), ncol=2)
  #' pos_ref = t(t(pos_ref)-colMeans(pos_ref))
  #' 
  #' # Create a new configuration by adding a perturbation to the previous one
  #' pos = pos_ref + matrix(rnorm(20, mean=1, sd=0.1), ncol=2)
  #' 
  #' # Compute the Procrustean Transform and inspect the results
  #' proc_pos = proc_crr(pos, pos_ref)
  #' plot(pos_ref, col="blue", pch=20, xlim=c(-1,3), ylim=c(-1,3))
  #' points(pos, col="red", pch=20)
  #' points(proc_pos, col="purple", pch=20)
  #' @export
  
  Z=t(t(Z)-colMeans(Z))
  
  A=t(Z)%*%(Z0%*%t(Z0))%*%Z
  eA=eigen(A,symmetric=T)
  Ahalf=eA$vec%*%diag(sqrt(eA$val))%*%t(eA$vec)
  
  t(t(Z0)%*%Z%*%solve(Ahalf)%*%t(Z))  
}


estimate_latent_positions = function (Y,W,
                                      procrustean = TRUE, 
                                      k=2,
                                      alpha=2,
                                      nscan=8*10^5, burn_in=5*10^5, odens=10^3,
                                      zdelta=1, z_norm_prior_mu=0, z_norm_prior_sd=10,
                                      adelta=.3, a_exp_prior_a=1, a_exp_prior_b=1,
                                      dynamic_plot = FALSE, dynamic_circles = FALSE,
                                      ...){
  #' @title BLSM simulation
  #' @description Core function of the BLSM package: run a simulation to obtain the positions of the network nodes 
  #' in the latent space for each sampled iteration.
  #' 
  #' The positions are simulated accordingly to the model assumptions, please refer to \link[BLSM]{BLSM} for further information. 
  #' The output of the function can be used to retrieve and compare specific iterations, observe their evolution or simply compute
  #' the average positions (more details in the descriptions and examples below). 
  #' 
  #' @param Y Adjacency matrix of the network
  #' @param W (Optional) BLSM Weight matrix of the network
  #' @param k Space dimensionality
  #' @param procrustean Boolean to include/exclude (\code{TRUE/FALSE}) the Procrustean Transform step in the algorithm. Set \code{TRUE} by default.
  #' @param alpha Starting value of the \eqn{\alpha} variable
  #' @param nscan Number of iterations
  #' @param burn_in Burn-in value (starting iterations to be discarded)
  #' @param odens Thinning: only 1 iteration every \code{odens} will be sampled and stored in the output
  #' @param zdelta Standard deviation of the Gaussian proposal for latent positions
  #' @param z_norm_prior_mu Mean of the Gaussian prior distribution for latent positions 
  #' @param z_norm_prior_sd Standard deviation of the Gaussian prior distribution for latent positions
  #' @param adelta The uniform proposal for \eqn{\alpha} is defined on the \eqn{[-adelta,+adelta]} interval
  #' @param a_exp_prior_a Shape parameter of the Gamma prior distribution for \eqn{\alpha}. As the value is usually set to 1 the prior is an exponential distribution.
  #' @param a_exp_prior_b Rate parameter of the Gamma prior distribution for \eqn{\alpha}. 
  #' @param dynamic_plot Boolean to plot dynamically the simulated positions (one update every \code{odens} iterations)
  #' @param dynamic_circles Boolean to add circles of radius \eqn{\alpha} to the dynamic plots
  #' @param \dots Additional parameters that can be passed to \link[BLSM]{plot_latent_positions}
  #' 
  #' @return Returns a "BLSM object" (\code{blsm_obj}), i.e. a list containing:
  #' \itemize{
  #' \item \code{Alpha }{\eqn{\alpha} values from the sampled iterations}
  #' \item \code{Likelihood }{Log-likelihood values from the sampled iterations}
  #' \item \code{Iterations }{Latent space coordinates from the sampled iterations. Latent positions are stored in a
  #' 3D array whose dimensions are given by (1: number of nodes, 2: space dimensionality, 3: number of iterations).
  #' In the non-Procrustean framework the latent distances are given instead of the positions: another 3D array is returned, whose dimensions
  #' are given by (1: number of nodes, 2: number of nodes, 3: number of iterations). The command needed in order to get the average values over the iterations for
  #' either the positions or the distances is \code{rowMeans(blsm_obj$Iterations, dims=2)} (see example below).}
  #' \item \code{StartingPositions }{Latent space coordinates right after the initialization step. In the non-Procrustean framework starting distances are given instead.}
  #' \item \code{Matrix }{Original matrices of the network (adjacency and BLSM weights)}
  #' \item \code{Parameters }{List of parameters specified during the call to \link[BLSM]{estimate_latent_positions}}
  #' }
  #' 
  #' @examples 
  #' \dontshow{
  #' blsm_obj_test_1 = estimate_latent_positions(example_adjacency_matrix, burn_in = 10^3, nscan = 3*10^3, odens=100)
  #' blsm_obj_test_2 = estimate_latent_positions(example_adjacency_matrix, procrustean=FALSE, burn_in = 10^3, nscan = 3*10^3, odens=100)
  #' }
  #' \dontrun{
  #'  # Procrustean version followed by clustering
  #'  blsm_obj = estimate_latent_positions(example_adjacency_matrix,  
  #'                           burn_in = 3*10^4, nscan = 10^5, dynamic_plot = TRUE)
  #'                           
  #'  avg_latent_positions = rowMeans(blsm_obj$Iterations, dims=2)                   
  #'  h_cl = hclust(dist(avg_latent_positions), method="complete")
  #'  n = 3
  #'  latent_space_clusters = cutree(h_cl, k=n)
  #'  print(latent_space_clusters)
  #'  plot(avg_latent_positions, col=rainbow(n)[latent_space_clusters], pch=20)
  #'  
  #'  # Non-Procrustean version followed by clustering                    
  #'  blsm_obj_2 = estimate_latent_positions(example_adjacency_matrix, procrustean=FALSE,
  #'                           burn_in = 3*10^4, nscan = 10^5)
  #'  avg_latent_distances = rowMeans(blsm_obj_2$Iterations, dims=2)
  #'  h_cl = hclust(as.dist(avg_latent_distances), method="complete")
  #'  n = 3
  #'  latent_space_clusters_2 = cutree(h_cl, k=n)
  #'  print(latent_space_clusters_2)
  #'                            
  #'  # Weighted network 
  #'  blsm_obj_3 = estimate_latent_positions(example_adjacency_matrix, example_weights_matrix, 
  #'                           burn_in = 10^5, nscan = 2*10^5, dynamic_plot = TRUE)
  #' }
  #' @export
  
  if (missing(W)) {
    W=Y*0+1
  }
  
  if (k==3 & dynamic_plot){
    if (!requireNamespace("rgl", quietly = TRUE)) {
      message("rgl package needed for the 3D plot. Please install it or set dynamic_plot=FALSE.")
      return(NULL)      
    }
  }
  
  params = list(procrustean=procrustean, 
                k=k,
                alpha=alpha,
                nscan=nscan, burn_in=burn_in, odens=odens,
                zdelta=zdelta, z_norm_prior_mu=z_norm_prior_mu, z_norm_prior_sd=z_norm_prior_sd,
                adelta=adelta, a_exp_prior_a=a_exp_prior_a, a_exp_prior_b=a_exp_prior_b)
  
  nscan = nscan + burn_in
  it_cont = 1
  create_window_flag=dynamic_plot
  
  rem = which(rowMeans(Y)==0)
  if(length(rem)>0) {
    Y=Y[-rem,-rem]
    W=W[-rem,-rem]
  }
  n=dim(Y)[1]   
  my_colors=blsm_colors(n)
  cc=(Y>0)+0
  D=dst(cc)
  Z=cmdscale(D, k)
  Z=t(t(Z)-colMeans(Z))
  
  tmp_opt=c(alpha,c(Z))
  tmp_opt=optim(tmp_opt,mlpY,Y=Y,W=W, method="SANN")$par 
  tmp_opt=optim(tmp_opt,mlpY,Y=Y,W=W,method="Nelder-Mead")$par
  
  tmp_opt=tmp_opt*2/(tmp_opt[1])  
  alpha=tmp_opt[1]
  
  Z_Proc = matrix(tmp_opt[-1],nrow=n,ncol=k)
  Z_Proc=t(t(Z_Proc)-colMeans(Z_Proc))
  Z = Z_Proc
  lpz = lpz_dist(Z)
  
  acc_a=0
  acc_z=0 
  Alpha=alpha       
  Lik=lpY(Y,lpz,alpha, W) 
  
  inputs = list(Adjacency = Y, Weight = W)
  
  if (procrustean){
    blsm_obj=list(Alpha=rep(NA, (nscan-burn_in)/odens), 
                  Likelihood=rep(NA, (nscan-burn_in)/odens),
                  Iterations=array(NA,dim=c(n,k,(nscan-burn_in)/odens)), 
                  StartingPositions=Z,
                  Matrix=inputs, 
                  Parameters=params)
    
    avg_Z_est = Z
    for(ns in 1:nscan){
      tmp = Z_up(Y,Z,W,alpha,zdelta, z_norm_prior_mu, z_norm_prior_sd)
      if(any(tmp!=Z)){
        acc_z=acc_z+sum(tmp!=Z)/(2*n*odens)
        tryCatch({
          Z = proc_crr(tmp,Z_Proc)
        }, error = function(e) {
          message("The matrix used to compute the Procrustean transformation is singular. \nIf changing the parameters doesn't solve the issue, please try to lower the space dimensionality.")
          return(blsm_obj)
          graphics.off()
        }
        )
      }
      lpz = lpz_dist(Z)
      tmp = alpha_up(Y,lpz,W,alpha,adelta,a_exp_prior_a,a_exp_prior_b)
      if(tmp!=alpha) {
        acc_a=acc_a+1/odens
        alpha=tmp
      }
      if (ns%%odens==0){
        
        if(ns>burn_in){
          if (create_window_flag){
            if (k==2) {
              dev.new(noRStudioGD = TRUE)
            } else if (k==3) {
              rgl::par3d(windowRect = c(50, 50, 800, 800))
            } else {
              message("Error: plot cannot be displayed since space dimensionality is bigger than 3.")
              dynamic_plot=FALSE
            }
            create_window_flag = FALSE
          }
          blsm_obj$Alpha[it_cont]= alpha
          lik = lpY(Y, lpz, alpha, W)
          blsm_obj$Likelihood[it_cont] = lik
          cat(ns-burn_in,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          blsm_obj$Iterations[,,it_cont] = Z
          avg_Z_est= ((it_cont-1)*avg_Z_est+Z)/it_cont 
          it_cont = it_cont+1 
          
          if (dynamic_plot){
            plot_latent_positions(blsm_obj, circles_2D = dynamic_circles, ...)
          }
        } else {
          if (ns==odens){
            cat("\nBeginning burn-in period...\n\n")
          }
          lik = lpY(Y, lpz, alpha, W)
          cat(ns,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          if (ns==burn_in){
            cat("\nBurn-in period ended.\n\nBeginning simulation...\n\n")
          }
        }
      }
    }
  } else {
    blsm_obj=list(Alpha=rep(NA, (nscan-burn_in)/odens), 
                  Likelihood=rep(NA, (nscan-burn_in)/odens),
                  Iterations=array(NA,dim=c(n,n,(nscan-burn_in)/odens)),
                  StartingDistances=-lpz_dist(Z),
                  Matrix=inputs, 
                  Parameters=params)
    
    for(ns in 1:nscan){
      Z_tmp=Z_up(Y,Z,W,alpha,zdelta, z_norm_prior_mu, z_norm_prior_sd)  
      if(any(Z_tmp!=Z)){
        acc_z=acc_z+sum(Z_tmp!=Z)/(2*n*odens)
        Z=Z_tmp
      }
      lpz = lpz_dist(Z)
      tmp=alpha_up(Y,lpz,W,alpha,adelta,a_exp_prior_a,a_exp_prior_b)
      if(tmp!=alpha) {
        acc_a=acc_a+1/odens
        alpha=tmp
      }
      if (ns%%odens==0){
        if( ns > burn_in){
          blsm_obj$Alpha[it_cont]= alpha
          lik = lpY(Y, lpz, alpha, W)
          blsm_obj$Likelihood[it_cont] = lik
          cat(ns-burn_in,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          blsm_obj$Iterations[,,it_cont] = -lpz
          it_cont = it_cont + 1
        } else {
          if (ns==odens){
            cat("\nBeginning burn-in period...\n\n")
          }
          lik = lpY(Y, lpz, alpha, W)
          cat(ns,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          if (ns==burn_in){
            cat("\nBurn-in period ended.\n\nBeginning simulation...\n\n")
          }
        }
      }
    }
    
  }
  cat("\nThe simulation ended successfully.")
  return(blsm_obj)
}


plot_traceplots_acf = function (blsm_obj, chosen_node=1,  coordinate=1, chosen_pair=c(1,2)){
  #' @title BLSM traceplots and ACF
  #' @description Traceplots and autocorrelation functions for the \eqn{\alpha} variable and a selected node (or pair of nodes in the non-Procrustean framework).
  #' 
  #' @param blsm_obj BLSM object obtained through \link[BLSM]{estimate_latent_positions}
  #' @param chosen_node Specified node for traceplot and autocorrelation function (Procrustean framework)
  #' @param coordinate Specified coordinate dimension from the n-dimensional latent space
  #' @param chosen_pair Specified pair of nodes for traceplot and autocorrelation function (non-Procrustean framework)
  #' 
  #' @examples 
  #' plot_traceplots_acf(example_blsm_obj, chosen_node=3, coordinate=1)
  #' 
  #'\dontrun{
  #'  # Run the simulation without Procrustean step
  #'  blsm_obj = estimate_latent_positions(example_adjacency_matrix, procrustean = FALSE, 
  #'                           burn_in = 3*10^4, nscan = 10^5)
  #'  
  #'  # Plot 
  #'  plot_traceplots_acf(blsm_obj, chosen_pair=c(2,5))
  #'}
  #' @export
  
  par(mfrow=c(2,2))
  if (blsm_obj$Parameters$procrustean==TRUE){
    plot(blsm_obj$Alpha,type="l",main="Alpha")
    plot(blsm_obj$Iterations[chosen_node,coordinate,],type="l", ylab=paste0("Coordinate ",coordinate), main = paste0("Node ",chosen_node))
    plot(acf(blsm_obj$Alpha,plot=F),main="Alpha")
    plot(acf(blsm_obj$Iterations[chosen_node,coordinate,],plot=F),main= paste0("Node ",chosen_node))
  } else {
    plot(blsm_obj$Alpha,type="l",main="Alpha")
    plot(blsm_obj$Iterations[chosen_pair[1],chosen_pair[2],],type="l", ylab="Euclidean distance", main = paste0("Nodes ",chosen_pair[1], " and ", chosen_pair[2]))
    plot(acf(blsm_obj$Alpha,plot=F),main="Alpha")
    plot(acf(blsm_obj$Iterations[chosen_pair[1],chosen_pair[2],],plot=F),main= paste0("Nodes ",chosen_pair[1], " and ", chosen_pair[2]))
  }
}


plot_latent_positions = function(blsm_obj, colors, points_size=0.1, labels_point_size=5, labels_point_color="yellow", 
                                 labels_text_size=1, labels_text_color="blue", circles_2D = FALSE){
  #' @title Base BLSM plot function
  #' @description Plot latent positions from a Procrustean simulation.
  #' 
  #' @param blsm_obj BLSM object obtained through \link[BLSM]{estimate_latent_positions}
  #' @param colors (Optional) Colors of the simulated coordinate points in the latent space. Internal default colors are used if the argument is missing. 
  #' @param points_size Size of the coordinate points
  #' @param labels_point_size Size of the label points
  #' @param labels_point_color Color of the label points
  #' @param labels_text_size Text size in the label points
  #' @param labels_text_color Text color in the label points
  #' @param circles_2D Plot circles of radius \eqn{\alpha} (see the model's main variables) centered around the label points
  #' 
  #' @examples 
  #' plot_latent_positions(example_blsm_obj, labels_point_color = "black", labels_text_color = "black")
  #' 
  #' plot_latent_positions(example_blsm_obj, circles_2D = TRUE)
  #' @export
  
  n=dim(blsm_obj$Iterations)[1]
  if (missing(colors)) {
    colors=blsm_colors(n)
  }
  
  if (blsm_obj$Parameters$procrustean==TRUE){
    if (dim(blsm_obj$Iterations)[2]==2){
      par(mfrow=c(1,1))
      par(mar=c(1,1,1,1))
      par(mgp=c(2,1,0))
      dev.hold()
      plot(blsm_obj$Iterations[,1,],blsm_obj$Iterations[,2,],pch=20,cex=points_size, col=colors, xlab="", ylab="", xaxt="n", yaxt="n")
      
      avg_Z_est=rowMeans(blsm_obj$Iterations, dims=2, na.rm = TRUE)
      
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          lines(avg_Z_est[c(i,j),1], avg_Z_est[c(i,j),2] ,lty=blsm_obj$Matrix$Adjacency[i,j], col="blue", lwd =2)
        }
      }
      
      points(avg_Z_est[,1],avg_Z_est[,2],xaxt="n",yaxt="n",xlab="",ylab="", col=labels_point_color,pch=20,cex=labels_point_size)
      text(avg_Z_est[,1],avg_Z_est[,2],labels(blsm_obj$Matrix$Adjacency)[[1]],col=labels_text_color,cex=labels_text_size)
      
      if (circles_2D){symbols(avg_Z_est, circles = rep(mean(blsm_obj$Alpha),n), add=TRUE, fg = colors, inches=F)}
      dev.flush()
    } else if (dim(blsm_obj$Iterations)[2]==3){
      if (!requireNamespace("rgl", quietly = TRUE)) {
        message("rgl package needed for the 3D plot. Please install it.")
        return(NULL)      
      }
      par(mfrow=c(1,1))
      par(mar=c(3,3,1,1))
      par(mgp=c(2,1,0))
      
      rgl::plot3d(blsm_obj$Iterations[,1,],blsm_obj$Iterations[,2,],blsm_obj$Iterations[,3,], size=points_size*10, 
             col=colors, xlab = "", ylab="", zlab="")
      rgl::rgl.viewpoint(theta=30, phi=10, fov=30)
      
      if (missing(avg_Z_est)){
        avg_Z_est=rowMeans(blsm_obj$Iterations, dims=2, na.rm = TRUE)
      }
      for(i in 1:(n-1)){ 
        for(j in (i+1):n){
          if (blsm_obj$Matrix$Adjacency[i,j]){
            rgl::lines3d( avg_Z_est[c(i,j),], col="blue", lwd =2)
          }
        }
      }
      
      rgl::points3d(avg_Z_est[,1],avg_Z_est[,2],avg_Z_est[,3],xaxt="n",yaxt="n",xlab="",ylab="", col=labels_point_color,size=labels_point_size*10)
      
      rgl::text3d(avg_Z_est[,1],avg_Z_est[,2],avg_Z_est[,3],labels(blsm_obj$Matrix$Adjacency)[[1]],col=labels_text_color,cex=labels_text_size)
    } else {
      message("Error: plot cannot be displayed since space dimensionality is bigger than 3.")
    }
  } else {
    message("Sampled latent positions are not available in the non-Procrustean framework.\n")
  }
}
