# @name comire.gibbs.continuous
#
# @title Gibbs sampler for CoMiRe model with continuous response.
# 
# @description Posterior inference via Gibbs sampler for CoMiRe model with continuous response.
# 
# @param y numeric vector for the response.
# @param x numeric vector for the covariate relative to the dose of exposure.
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
# @param verbose logical, if \code{TRUE} a message on the status of the MCMC algorithm is printed to the console. Default is \code{TRUE}.
# 
#' @importFrom stats rgamma dnorm rnorm
#' @importFrom truncnorm rtruncnorm

.comire.gibbs.continuous <-
  function(y, x, grid=NULL, mcmc, prior, state=NULL, seed, max.x=max(x), verbose = TRUE){
    # prior: mu.theta, k.theta, eta(Jx1), alpha(Hx1), a, b, H, J)

    # internal working variables
    n <- length(y)
    H <- prior$H
    J <- prior$J 
    print_now <- c(mcmc$nb + 1:mcmc$ndisplay*(mcmc$nrep)/mcmc$ndisplay)

    # create the objects to store the MCMC samples
    nu0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)
    th0 <- tau0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)
    th1 <- tau1 <- nu1 <- rep(NA, mcmc$nrep+mcmc$nb)
    w <- matrix(NA, mcmc$nrep+mcmc$nb, length(prior$eta)) # J?

    # initialize each quantity

    ## parameters of the model
    if(is.null(state))
    {
      set.seed(seed)
      nu0[1,] <- rep(1/H, H)
      nu1[1] <- 1
      tau0[1,] <- stats::rgamma(H, prior$a, prior$b)
      th0[1,] <- stats::rnorm(H, prior$mu.theta, sqrt(prior$k.theta))
      tau1[1] <- stats::rgamma(1, prior$a, prior$b)
      th1[1] <- truncnorm::rtruncnorm(1, a=-Inf, b=min(th0[1,]), prior$mu.theta, sqrt(prior$k.theta))
      w[1,] <- prior$eta/sum(prior$eta) # ?
    }
    else
    {
      nu0[1,] <- state$nu0
      nu1[1] <- 1
      th0[1,] <- state$th0
      tau0[1,] <- state$tau0
      th1[1] <- state$th1
      tau1[1] <- state$tau1
      w[1,] <- state$w
    }

    ## quantity of interest
    if(is.null(grid$grids))
    {
      x.grid <- seq(0, max.x, length=100)
      y.grid <- seq(min(y)-sqrt(var(y)), max(y) + sqrt(var(y)), length = 100)
      beta_x <- matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))
      f0 <- matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
      f1 <- matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
    }
    else
    {
      x.grid = grid$xgrid
      y.grid = grid$ygrid
      beta_x = matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))
      f0 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
      f1 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
    }

    ## basis expansion
    knots <- seq(min(x)+1, max.x, length=prior$J-3)
    basisX <- function(x) iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
    
    phiX <- basisX(x)
    phi.grid <- basisX(x.grid)

    ## beta_i is the interpolating function evaluated at x_i
    beta_i =  as.double(phiX %*% w[1, ])

    ## f0i and f1i
    f0i = .mixdensity_C(y=y, pi=nu0[1,], mu=th0[1,], tau=tau0[1,])
    #sapply(1:n, mixdensity, y=y, nu=nu0[1,], mu=th0[1,], tau=tau0[1,])
    f1i = stats::dnorm(y, th1[1], sqrt(1/tau1[1]))

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
      d <- rbinom(n, 1, prob=(beta_i*f1i)/((1-beta_i)*f0i + beta_i*f1i))

      # 2. Update b_i from the multinomial
      b <- .labelling_b_C(w=w[ite-1,], phi=phiX, f0i=f0i, f1i=f1i)
      
      # 3. Update c_i, marginalizing over b_i and d_i, from the multinomial
      ind0 <- c(1:n)[d==0]
      ind1 <- c(1:n)[d==1]
      c0 <- .labelling_c_C(y=y[ind0], logpi=log(nu0[ite-1,]), mu=th0[ite-1,], tau=tau0[ite-1,])
      
      # 4. Update the mixture weights sampling from dirichlet
      n_0h <- table(factor(c0, levels=1:H))
      nu0[ite,] <- as.double(rdirichlet(1, as.double(prior$alpha + n_0h)))
      nu1[ite] <- 1

      # 5. Update w from the Dirichlet and obtain an updated function beta_i
      eta.post <- as.double(prior$eta + table(factor(b, levels=1:length(w[ite-1,]))))
      w[ite, ] <- as.double(rdirichlet(1, eta.post))
      beta_i <- as.numeric(phiX %*% w[ite, ])
      beta_i[beta_i>1] <- 1
      beta_i[beta_i<0] <- 0

      # 6. Updated mu and tau from the usual normal inverse gamma

      ## in cluster 0
      n_h <- table(factor(c0, levels=1:H))
      hat_a <- prior$a + n_h/2
      mean_h <- tapply(y[ind0], factor(c0, levels=1:H), mean)
      mean_h[n_h==0]  <- 0
      deviance_h <- sapply(1:H, .pssq_gaussian, data=y[ind0], cluster = c0, locations = th0[ite-1,])
      hat_b <- prior$b + 0.5*(deviance_h)
      tau0[ite, ] <- stats::rgamma(H, hat_a, hat_b)
      hat_k.theta <- 1/(1/prior$k.theta + n_h*tau0[ite,])
      hat_mu <- hat_k.theta*(1/prior$k.theta*prior$mu.theta + n_h*mean_h*tau0[ite, ])
      th0[ite, ] <- truncnorm::rtruncnorm(H, a=th1[ite-1], b=Inf, hat_mu, sqrt(hat_k.theta))

      ## in cluster 1
      n_h <- sum(d==1)
      hat_a <- prior$a + n_h/2
      mean_h <- mean(y[ind1])
      mean_h[n_h==0] <- 0
      deviance_h <- sum((y[ind1] -  th1[ite-1])^2)
      hat_b <- prior$b + 0.5*(deviance_h)
      tau1[ite] <- stats::rgamma(1, hat_a, hat_b)
      hat_k.theta <- 1/(1/prior$k.theta + n_h*tau1[ite])
      hat_mu <- hat_k.theta*(1/prior$k.theta*prior$mu.theta + n_h*mean_h*tau1[ite])
      th1[ite] <- truncnorm::rtruncnorm(1, a=-Inf, b=min(th0[ite, ]), hat_mu, sqrt(hat_k.theta))

      # update the values of the densities in the observed points
      f0i <- .mixdensity_C(y=y, pi=nu0[ite,], mu=th0[ite,], tau=tau0[ite,])
      #sapply(1:n, mixdensity, y=y, nu=nu0[ite,], mu=th0[ite,], tau=tau0[ite,])
      f1i <- stats::dnorm(y, th1[ite], sqrt(1/tau1[ite]))

      # 7. compute some posterior quanties of interest
      beta_x[ite, ] <- phi.grid %*% w[ite, ]
      f0[ite, ] <- .mixdensity_C(y=y.grid, pi=nu0[ite,], mu=th0[ite,], tau=tau0[ite,])
      #sapply(1:length(y.grid), mixdensity, y=y.grid, nu=nu0[ite,], mu=th0[ite,], tau=tau0[ite,])
      f1[ite, ] <- stats::dnorm(y.grid, th1[ite], sqrt(1/tau1[ite]))

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

    # output
    output <- list(
      post.means = list(beta=post.mean.beta, w=post.mean.w,
                       th0=post.mean.th0, tau0=post.mean.tau0,
                       th1=post.mean.th1, tau1=post.mean.tau1,
                       f0=post.mean.f0, f1=post.mean.f1, nu0=post.mean.nu0, nu1=post.mean.nu1),
      ci = list( beta=ci.beta, w=ci.w,
                 th0=ci.th0, tau0=ci.tau0,
                 th1=ci.th1, tau1=ci.tau1,
                 f0=ci.f0,f1=ci.f1, nu0=ci.nu0#, nu1=ci.nu1
                 ),
      mcmc = list(beta=beta_x, w=w, th0=th0, tau0=tau0, th1=th1, tau1=tau1, f0=f0, f1=f1, nu0=nu0, nu1=nu1))
      output
  }
