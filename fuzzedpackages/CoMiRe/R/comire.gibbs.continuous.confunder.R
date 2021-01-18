# @name comire.gibbs.continuous.confunder
#
# @title Gibbs sampler for CoMiRe model with continuous response and one confounder variable.
# 
# @description Posterior inference via Gibbs sampler for CoMiRe model with continuous response and one confounder variable.
# 
# @param y numeric vector for the response.
# @param x numeric vector for the covariate relative to the dose of exposure.
# @param z numeric vector for the confunders; a vector if there is only one confounder or a matrix for two or more confunders.
# @param grid a list giving the parameters for plotting the posterior mean density and the posterior mean \eqn{\beta(x)} over finite grids.
# @param mcmc a list giving the MCMC parameters. It must include the following integers: \code{nb} giving the number of burn-in iterations, \code{nrep} giving the total number of iterations, \code{thin} giving the thinning interval, \code{ndisplay} giving the multiple of iterations to be displayed on screen while the algorithm is running (a message will be printed every \code{ndisplay} iterations).
# @param prior a list containing the values of the hyperparameters. 
# It must include the following values: 
# \itemize{
# \item \code{mu.theta}, the prior mean \eqn{\mu_\theta} for each location parameter \eqn{\theta_{0h}}{\theta_0h} and \eqn{\theta_1}, 
# \item \code{k.theta}, the prior variance \eqn{k_\theta} for each location paramter \eqn{\theta_{0h}}{\theta_0h} and \eqn{\theta_1}, 
# \item \code{mu.gamma} (if \code{p} confounding covariates are included in the model) a \code{p}-dimentional vector of prior means \eqn{\mu_\gamma}{\mu_gamma} of the parameters \eqn{\gamma} corresponding to the confounders,
# \item \code{k.gamma}, the prior variance \eqn{k_\gamma}{k_gamma} for parameter corresponding to the confounding covariate (if \code{p=1}) or \code{sigma.gamma} (if \code{p>1}), that is the covariance matrix \eqn{\Sigma_\gamma}{\Sigma_gamma} for the parameters corresponding to the \code{p} confounding covariates; this must be a symmetric positive definite matrix.
# \item \code{eta}, numeric vector of size \code{J} for the Dirichlet prior on the beta basis weights, 
# \item \code{alpha}, prior for the mixture weights,
# \item \code{a} and \code{b}, prior scale and shape parameter for the gamma distribution of each precision parameter, 
# \item \code{J}, parameter controlling the number of elements of the I-spline basis,
# \item \code{H}, total number of components in the mixture at \eqn{x_0}.
# }
# @param state if \code{family="continuous"}, a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis or if we want to start the MCMC algorithm from some particular value of the parameters.
# @param seed seed for random initialization.
# @param max.x maximum value allowed for \code{x}.
# @param z.val optional numeric vector containing a fixed value of interest for each of the confounding covariates to be used for the plots. Default value is \code{mean(z)} for numeric covariates or the mode for factorial covariates.
# @param verbose logical, if \code{TRUE} a message on the status of the MCMC algorithm is printed to the console. Default is \code{TRUE}.
#
#' @importFrom stats rgamma dnorm
#' @importFrom truncnorm rtruncnorm 

.comire.gibbs.continuous.confunder <-function(y, x, z, grid=NULL, mcmc, prior, state=NULL, seed, 
                             max.x=ceiling(max(x)), z.val=NULL, verbose = TRUE){
  
    if(is.null(z.val)){
      if(!is.factor(z)){
        z.val <- mean(z)
        }
      else{
        z.val <- mode(z)
      }
    }
    
    # internal working variables
    n <- length(y)
    H <- prior$H 
    J <- prior$J
    print_now <- c(mcmc$nb + 1:mcmc$ndisplay*(mcmc$nrep)/mcmc$ndisplay)
    
    # create the objects to store the MCMC samples
    w <- matrix(NA, mcmc$nrep+mcmc$nb, J)
    nu0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)         
    th0 <- tau0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)
    nu1 <- th1 <- tau1 <- rep(NA, mcmc$nrep+mcmc$nb)
    ga <- rep(NA, mcmc$nrep+mcmc$nb)
    
    # initialize each quantity
    
    ## parameters of the model
    if(is.null(state))
    {
      set.seed(seed)
      w[1,] <- prior$eta/sum(prior$eta)
      nu0[1,] <- rep(1/H, H)
      tau0[1,] <- stats::rgamma(H, prior$a, prior$b)
      th0[1,] <- rep(0, H)
      nu1[1] <- 1
      tau1[1] <- stats::rgamma(1, prior$a, prior$b)
      th1[1] <- truncnorm::rtruncnorm(1, a=-Inf, b=min(th0[1,]), prior$mu.theta, sqrt(prior$k.theta))
      ga[1]<- 0

    }
    else
    {
      w[1,] <- state$w
      nu0[1,] <- state$nu0
      nu1[1] <- 1
      th0[1,] <- state$th0
      tau0[1,] <- state$tau0
      th1[1] <- state$th1
      tau1[1] <- state$tau1
      ga[1] <- state$ga
    }
    
    ## quantity of interest
    if(is.null(grid$grids))
    {
      x.grid <- seq(0, max.x, length=100)
      y.grid <- seq(min(y)-sqrt(var(y)), max(y)+sqrt(var(y)), length = 100)	
      beta_x <- matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))	
      f0 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
      f1 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))

    }
    else
    {
      x.grid <- grid$xgrid
      y.grid <- grid$ygrid
      beta_x <- matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))
    }
    
    ## basis expansion
    knots <- seq(min(x)+1, max.x, length=prior$J-3)
    basisX <- function(x) iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
    phiX <- basisX(x)
    phi.grid <- basisX(x.grid)
    
    ## beta_i is the interpolating function evaluated at x_i
    beta_i <- as.double(phiX %*% w[1, ])
    
    ## f0i and f1i
    f0i <- sapply(1:n, .mixdensity_uni, y=y, z=z, nu=nu0[1,], theta=th0[1,], tau=tau0[1,], ga=ga[1])
    f1i <- stats::dnorm(y, (th1[1]+z*ga[1]) , sqrt(1/tau1[1]))
    
    # start the MCMC simulation 
    set.seed(seed)
    for(ite in 2:(mcmc$nrep+mcmc$nb))
    {
      # 0. Print the iteration
      if(verbose)
      {
        if(ite==mcmc$nb) cat("Burn in done\n")
        if(ite %in% print_now) cat(ite, "iterations over",
                                   mcmc$nrep+mcmc$nb, "\n")
      }
      
      # 1. Update d_i marginalising out b_i from
      d = rbinom(n, 1, prob=(beta_i*f1i)/((1-beta_i)*f0i + beta_i*f1i))
      
      # 2. Update b_i from the multinomial 
      b = sapply(1:n, .labelling_b_uni, w[ite-1,], phi=phiX, f0i=f0i, f1i=f1i)
      
      # 3. Update c_i, marginalizing over b_i and d_i, from the multinomial 
      ind0 <- c(1:n)[d==0]
      ind1 <- c(1:n)[d==1]
      c <- sapply(ind0, .labelling_c_uni, y=y, z=z, nu=nu0[ite-1,], theta=th0[ite-1,], tau=tau0[ite-1,], ga=ga[ite-1])
      
      # 4. Update the mixture weights sampling from Dirichlet
      n_h <- table(factor(c, levels=1:H))
      nu0[ite,] <- as.double(rdirichlet(1, as.double(prior$alpha + n_h)))
      
      # 5. Update w from the Dirichlet and obtain an update function beta_i
      n_j <- table(factor(b, levels=1:J))  
      w[ite, ] <- as.double(rdirichlet(1, as.double(prior$eta + n_j)))
      
      beta_i <- as.numeric(phiX %*% w[ite, ])
      beta_i[beta_i>1] <- 1
      beta_i[beta_i<0] <- 0

      # 6. Update gamma
      n_0h <- n_h 
      n_1h <- length(ind1)
      
      z_mean_0h <- tapply(z[ind0], factor(c, levels=1:H), mean); z_mean_0h[is.na(z_mean_0h)]<- 0
      z_mean_1h <- mean(z[ind1]); z_mean_1h[is.na(z_mean_1h)]<- 0
      
      zy_sum_0h <- sapply(1:H, .psdp_uni, y=y[ind0], z=z[ind0], cluster=c)
      zy_sum_1h <- sum(y[ind1]*z[ind1])
      
      z2_0h <- tapply(z[ind0], factor(c, levels=1:H), function(x) sum(x^2)); z2_0h[is.na(z2_0h)]<- 0
      z2_1h <- sum(z[ind1]^2); z2_1h[is.na(z2_1h)]<- 0
      
      post_var <- (1/prior$k.gamma + sum(tau0[ite-1,]*z2_0h) + tau1[ite-1]*z2_1h )^-1 
      
      post_mean <- post_var * ( prior$mu.gamma/prior$k.gamma + ( 
        sum( tau0[ite-1,]*(zy_sum_0h - th0[ite-1,]*n_0h*z_mean_0h ) ) 
        + tau1[ite-1]*(zy_sum_1h - th1[ite-1]*n_1h*z_mean_1h ) ) ) 
      
      ga[ite] <- rnorm(1, mean=post_mean, sd=sqrt(post_var))

      # 7. Update theta and tau 
      
      ## in cluster 0: m=0
      hat_a <- prior$a + n_0h/2
      y_mean_0h <- tapply(y[ind0], factor(c, levels=1:H), mean); y_mean_0h[n_0h==0]  = 0
      deviance_0h <- sapply(1:H, .pssq_uni, y=y[ind0], z=z[ind0], cluster=c, theta=th0[ite-1,], gamma=ga[ite] )
      hat_b <- prior$b + 0.5*(deviance_0h)
      tau0[ite, ] <- stats::rgamma(H, hat_a, hat_b)
      
      hat_kappa <- (1/prior$k.theta + n_0h*tau0[ite,])^-1
      hat_mu <- hat_kappa * (prior$mu.theta/prior$k.theta+n_0h*tau0[ite, ]*(y_mean_0h-ga[ite]*z_mean_0h) )
      th0[ite, ] <- truncnorm::rtruncnorm(H, a=th1[ite-1], b=Inf, hat_mu, sqrt(hat_kappa))

      ## in cluster 1: m=1
      hat_a <- prior$a + n_1h/2
      y_mean_1h <- mean(y[ind1])
      y_mean_1h[n_1h==0] <- 0
      deviance_1h <- sum( (y[ind1] - (th1[ite-1]+z[ind1]*ga[ite]))^2 )
      hat_b <- prior$b + 0.5*(deviance_1h)
      tau1[ite] <- stats::rgamma(1, hat_a, hat_b)
      
      hat_kappa <- (1/prior$k.theta + n_1h*tau1[ite])^-1
      hat_mu <- hat_kappa * (prior$mu.theta/prior$k.theta + n_1h*tau1[ite]*(y_mean_1h-ga[ite]*z_mean_1h) )
      th1[ite] <- truncnorm::rtruncnorm(1, a=-Inf, b=min(th0[ite, ]), hat_mu, sqrt(hat_kappa))

      # update the values of the densities in the observed points
      f0i <- sapply(1:n, .mixdensity_uni, y=y, z=z, nu=nu0[ite,], theta=th0[ite,], 
                   tau=tau0[ite,], ga=ga[ite])
      f1i <- stats::dnorm(y, (th1[ite]+z*ga[ite]), sqrt(1/tau1[ite]))
      
      # 7. compute some posterior quantities of interest
      beta_x[ite, ] <- phi.grid %*% w[ite, ]
      f0[ite,] <- sapply(1:length(y.grid), .mixdensity_uni, y=y.grid, z=rep(z.val,length(y.grid)), nu=nu0[ite,], theta=th0[ite,], tau=tau0[ite,], ga=ga[ite])
      f1[ite,] <- stats::dnorm(y.grid, (th1[ite]+rep(z.val,length(y.grid))*ga[ite]) , sqrt(1/tau1[ite]))

      
    }
    
    post.mean.beta <- colMeans(beta_x[-c(1:mcmc$nb),])
    post.mean.w <- colMeans(w[-c(1:mcmc$nb),])
    post.mean.th0 <- colMeans(th0[-c(1:mcmc$nb),])
    post.mean.tau0 <- colMeans(tau0[-c(1:mcmc$nb),])
    post.mean.th1 <- mean(th1[-c(1:mcmc$nb)])
    post.mean.tau1 <- mean(tau1[-c(1:mcmc$nb)])
    post.mean.f0 <- colMeans(f0[-c(1:mcmc$nb),])
    post.mean.f1 <- colMeans(f1[-c(1:mcmc$nb),])
    post.mean.nu0 <- colMeans(nu0[-c(1:mcmc$nb),])
    post.mean.nu1 <- mean(nu1[-c(1:mcmc$nb)])
    post.mean.ga <- mean(ga[-c(1:mcmc$nb)])

    ci.beta <- apply(beta_x[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
    ci.w<- apply(w[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
    ci.th0 <- apply(th0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
    ci.tau0 <- apply(tau0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
    ci.th1 <- quantile(th1[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
    ci.tau1 <- quantile(tau1[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
    ci.f0 <- apply(f0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
    ci.f1 <- apply(f1[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
    ci.nu0 <- apply(nu0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
    #ci.nu1 <- quantile(nu1[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
    ci.ga <- quantile(ga[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
    
    # output
    output <- list(
      post.means = list(beta=post.mean.beta, w=post.mean.w, 
                       th0=post.mean.th0, tau0=post.mean.tau0, 
                       th1=post.mean.th1, tau1=post.mean.tau1, 
                       f0=post.mean.f0, f1=post.mean.f1,
                       nu0=post.mean.nu0, nu1=post.mean.nu1,
                       ga=post.mean.ga),
      ci = list(beta=ci.beta, w=ci.w, th0=ci.th0, tau0=ci.tau0, 
                th1=ci.th1, tau1=ci.tau1, f0=ci.f0, f1=ci.f1, nu0=ci.nu0, #nu1=ci.nu1, 
                ga= ci.ga),
      mcmc = list(beta=beta_x, w=w, th0=th0, tau0=tau0, th1=th1, tau1=tau1,
                  f0=f0, f1=f1, nu0=nu0, nu1=nu1, ga=ga)    
      )
  list(out = output, z.val = z.val)
}
