# create a comire object class

#' @name as.classCoMiRe
#' 
#' @export 
#' 
#' @usage as.classCoMiRe(call = NULL, out = NULL, z = NULL, z.val = NULL, f0 = NULL, f1 = NULL, 
#' nrep, nb, bin = FALSE, univariate = TRUE)
#'
#' @title classCoMiRe class constructor
#' 
#' @description A constructor for the \code{classCoMiRe} class. The class \code{classCoMiRe} is a named list containing
#' the output of the posterior estimation of CoMiRe model implemented in \code{comire.gibbs}
#'  
#' @param call a formula for \code{comire.gibbs}.
#' @param out an output of \code{comire.gibbs}.
#' @param z optional numeric vector or matrix for the confounding covariates.
#' @param z.val optional numeric vector containing a fixed value of interest for each of the confounding covariates to be used for the plots. Default value is \code{mean(z)} for numeric covariates or the mode for factorial covariates.
#' @param f0,f1 optional matrices containing simulated values of the mixture densities at low and high dose exposure; default values are simulated with \code{comire.gibbs}. It is possible to change these for differente fixed values of \code{z}: see \code{predict_new_z} function.
#' @param nrep integer giving the total number of iterations used in \code{comire.gibbs}.
#' @param nb integer giving the number of burn-in iterations used in \code{comire.gibbs}.
#' @param bin logical. It is \code{TRUE} if y is drawn for a binomial distribution.
#' @param univariate logical. It is \code{TRUE} if the model is univariate.
#' @author Antonio Canale, Arianna Falcioni
#' 


as.classCoMiRe <- function(
  call = NULL,
  out = NULL,
  z = NULL,
  z.val = NULL,
  f0 = NULL,
  f1 = NULL,
  nrep= NULL, 
  nb= NULL,
  bin = FALSE,
  univariate = TRUE
){
  value <- list(
    call = call, 
    post.means = out$post.means,
                ci = out$ci,
                mcmc = out$mcmc,
                z = z,
                z.val = z.val,
                f0 = if(is.null(f0)) out$mcmc$f0 else f0,
                f1 = if(is.null(f1)) out$mcmc$f1 else f1,
                nrep= nrep, 
                nb= nb,
                bin = bin,
                univariate = univariate)
  attr(value, "class") <- "classCoMiRe"
  value
}

# Function for f0 and f1 with non-default values of z.val
#' @name predict_new_z
#' 
#' @importFrom stats dnorm var
#' @export 
#' 
#' @usage predict_new_z(fit, y, z.val)
#'
#' @title comire.gibbs for different fixed values of z
#' 
#' @description This function computes the predicted values of the density al low 
#' dose \code{f_0} and of the density at high dose \code{f_{inf}}, for fixed values 
#' of the confounders \code{z}. 
#' 
#' @param fit the output of \code{comire.gibbs} opportunely trasformed in \code{classCoMiRe} class
#' @param y numeric vector for the response used in \code{comire.gibbs}
#' @param z.val optional numeric vector containing a fixed value of interest for each of the confounding covariates to be used for the plots. Default value is \code{mean(z)} for numeric covariates or the mode for factorial covariates.
#' 
#' @return An object of class \code{classCoMiRe}.
#' 
#' @author Antonio Canale, Arianna Falcioni
#' 
#' 
predict_new_z <- function(fit, y, z.val){
  if(fit$univariate){
    y.grid <- seq(min(y)-sqrt(stats::var(y)), max(y)+sqrt(stats::var(y)), length = 100)	
    f0 = matrix(NA, fit$nrep+fit$nb , length(y.grid))
    f1 = matrix(NA, fit$nrep+fit$nb , length(y.grid))
    for(ite in 2:(fit$nrep+fit$nb)){
      f0[ite,] <- sapply(1:length(y.grid), .mixdensity_uni, y=y.grid, z=rep(z.val,length(y.grid)), nu=fit$mcmc$nu0[ite,], theta=fit$mcmc$th0[ite,], tau=fit$mcmc$tau0[ite,], ga=fit$mcmc$ga[ite])
      f1[ite,] <- stats::dnorm(y.grid, (fit$mcmc$th1[ite]+rep(z.val,length(y.grid))*fit$mcmc$ga[ite]) , sqrt(1/fit$mcmc$tau1[ite]))
    }
    as.classCoMiRe(
      call = fit$call,
      out=list(post.means=fit$post.means, ci=fit$ci, mcmc=fit$mcmc), 
                 z=fit$z, 
                 z.val=z.val, 
                 f0=f0, 
                 f1=f1, 
                 nrep=fit$nrep, 
                 nb=fit$nb, 
                 bin=fit$bin, 
                 univariate = fit$univariate)
  }
  else{
    y.grid <- seq(min(y)-sqrt(stats::var(y)), max(y)+sqrt(stats::var(y)), length = 100)	
    z_val<- matrix(rep(z.val, length(y.grid)), ncol=length(z.val), byrow=T) #?
    f0 = matrix(NA, fit$nrep+fit$nb, length(y.grid))
    f1 = matrix(NA, fit$nrep+fit$nb, length(y.grid))
    for(ite in 2:(fit$nrep+fit$nb)){
      f0[ite,] <- sapply(1:length(y.grid), .mixdensity_multi, y=y.grid, z=z_val, 
                         nu=fit$mcmc$nu0[ite,], theta=fit$mcmc$th0[ite,], 
                         tau=fit$mcmc$tau0[ite,], ga=fit$mcmc$ga[ite,])
      f1[ite,] <- stats::dnorm(y.grid,(fit$mcmc$th1[ite]+z_val%*%fit$mcmc$ga[ite,]) , sqrt(1/fit$mcmc$tau1[ite]))
    }
    as.classCoMiRe(call = fit$call,
                   out=list(post.means=fit$post.means, ci=fit$ci, mcmc=fit$mcmc), 
                 z=fit$z, 
                 z.val=z.val, 
                 f0=f0, 
                 f1=f1, 
                 nrep=fit$nrep, 
                 nb=fit$nb, 
                 bin=fit$bin, 
                 univariate = fit$univariate)
  }
}



# -----------------------------------------------------------------------
# PLOT 1: postetrior predictive check
# -----------------------------------------------------------------------

#' @rdname post.pred.check
#' 
#' @importFrom ggplot2 theme_set theme_bw ggplot aes geom_line labs geom_point coord_cartesian
#' @importFrom splines2 iSpline
#' @importFrom stats pnorm rbinom
#' @importFrom KernSmooth locpoly
#' @importFrom NonpModelCheck localpoly.reg
#' @export 
#' 
#' @usage post.pred.check(y, x, z, fit, mcmc, J=10, H=10, a, max.x=max(x), 
#' xlim=c(0, max(x)), bandwidth = 20, oneevery = 20)
#'
#' @title Posterior predictive check plot
#' 
#' @description A plot for an object of \code{classCoMiRe} class. The plot is a goodness-of-fit assessment of CoMiRe model. 
#' If \code{family = 'continuous'}, a smoothed empirical estimate of F(a|x,z) = pr(y < a | x,z) is computed from the observed data (black line) 
#' and from some of the data sets simulated from the posterior predictive distribution in  the \code{fit} object (grey lines).
#' If \code{family = 'binary'}, a smoothed empirical estimate of the proportion of events (black line) and of the smoothed empirical 
#' proportion of data simulated from the posterior predictive distribution in the \code{fit} object (grey lines). 
#' In the x axis are reported the observed exposures.
#'  
#' @param y numeric vector for the response used in \code{comire.gibbs}
#' @param x numeric vector for the covariate relative to the dose of exposure used in \code{comire.gibbs}
#' @param z optional numeric vector or matrix for the confounding covariates.
#' @param fit the output of \code{comire.gibbs} opportunely trasformed in \code{classCoMiRe} class
#' @param mcmc a list giving the MCMC parameters
#' @param J parameter controlling the number of elements of the I-spline basis
#' @param H total number of components in the mixture at \eqn{x_0}
#' @param a threshold of clinical interest to compute the F(a|x,z)
#' @param max.x maximum value allowed for x
#' @param xlim numeric vectors of length 2, giving the x coordinates ranges for the plot
#' @param bandwidth the kernel bandwidth smoothing parameter
#' @param oneevery integer number representing how many MCMC draws to plot in the posterior predictive check. It draws one sample every \code{oneevery}.
#' 
#' @author Antonio Canale, Arianna Falcioni
#' 
#' @examples{
#' data(CPP)
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' premature <- as.numeric(gestage<=37)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## too few iterations to be meaningful. see below for safer and more comprehensive results
#' 
#' mcmc <- list(nrep=10, nb=2, thin=1, ndisplay=4) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit.dummy <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=1, max.x=180)
#'                      
#' post.pred.check(y = gestage, x = dde, fit = fit.dummy, mcmc = mcmc, J = 10, H = 10, a = 37, 
#'                 max.x = max(dde), xlim = c(0,150), oneevery = 4)
#'                 
#' \donttest{
#' ## safer procedure with more iterations (it may take some time)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## Fit the model for continuous y 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit1 <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#' 
#' post.pred.check(y = gestage, x = dde, fit = fit1, mcmc = mcmc, J = 10, H = 10, a = 37, 
#'                 max.x = max(dde), xlim = c(0,150))
#' 
#' }
#' }

post.pred.check <- function(y, x, z, fit, mcmc, J=10, H=10, a, max.x=max(x), 
                            xlim=c(0, max(x)), bandwidth=20, oneevery = 20)
{
  if(oneevery>(mcmc$nrep/mcmc$thin)) stop(paste("Too few MCMC iterations to print one every", oneevery, "iterations"))
  ggplot2::theme_set(ggplot2::theme_bw())
  if(!fit$bin){
    if(is.null(fit$z)){
      index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
      knots <- seq(min(x)+1, max.x, length=J-3)
      phi <- function(x) splines2::iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
      res <- rep(NA, length(index))
      tailp.mcmc1 <- function(i, pi0, mu0, tau0, mu1, tau1, theta, phi_x)
        {
        beta_x <-  as.double(phi_x %*% theta[i,])
        res <-  beta_x * stats::pnorm(a, mu1[i], 1/sqrt(tau1[i]))
        for(h in 1:ncol(pi0))
          {
          res <- res + (1-beta_x)*pi0[i,h] * stats::pnorm(a, mu0[i,h], 1/sqrt(tau0[i,h]))     
          }
        below <- stats::rbinom(length(res), 1, res)
        smoothed <- KernSmooth::locpoly(x=x, y=below, degree=0, bandwidth = bandwidth, 
                        gridsize = 100, range.x = c(0, max.x))$y
        }
      
      belowa.comire <- sapply(index, tailp.mcmc1, fit$mcmc$nu0, fit$mcmc$th0, fit$mcmc$tau0,
                fit$mcmc$th1, fit$mcmc$tau1, fit$mcmc$w, phi(x))
      
      belowa.true <- KernSmooth::locpoly(x=x, y=y<a, degree=0, bandwidth = bandwidth, gridsize = 100, range.x = c(0, max.x))$y
      
      ppc.data <- data.frame(cbind(x=seq(0,max.x, length=100), 
                                   Fx=c(c(belowa.comire[,1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery])), 
                                   repl=c(rep(1:length(1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery),each=100))))
      
      ppc <-  ggplot2::ggplot(ppc.data, ggplot2::aes(x=.data$x, y=.data$Fx)) + 
        ggplot2::geom_line(alpha=0.25, ggplot2::aes(group=factor(.data$repl)), col="grey") +
        ggplot2::geom_line(data=data.frame(x=seq(0,max.x, length=100), y=belowa.true), ggplot2::aes(x=.data$x, y=.data$y), col=1)+
        ggplot2::labs(x = "", y = expression(F(a *"| x")*" | data")) +
        ggplot2::geom_point(data=data.frame(x, zero=rep(0,length(x))), ggplot2::aes(.data$x, .data$zero), alpha=1, cex=.5, pch="|", na.rm=TRUE) 
      
      ppc + ggplot2::coord_cartesian(ylim=c(0,1), xlim=xlim) 
      
    }
    else {
      if(fit$univariate){
        index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
        knots <- seq(min(x)+1, max.x, length=J-3)
        phi <- function(x) splines2::iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
        res <- rep(NA, length(index))
        
        tailp.mcmc2 <- function(i, pi0, th0, tau0, th1, tau1, w, ga, phi_x)
        {
          beta_x <-  as.double(phi_x %*% w[i,])
          res <-  beta_x * stats::pnorm(a, th1[i]+ga[i]*fit$z, 1/sqrt(tau1[i]))
          
          for(h in 1:ncol(pi0))
          {
            res <- res + (1-beta_x) * pi0[i,h] * stats::pnorm(a, th0[i,h]+ga[i]*fit$z , 1/sqrt(tau0[i,h]))     
          }
          
          below <- stats::rbinom(length(res), 1, res) 
          smoothed <- NonpModelCheck::localpoly.reg( X=as.matrix(cbind(x,fit$z)), Y=below, 
                                     points=as.matrix(cbind(seq(0, max.x, length=100), 
                                                            seq(min(fit$z), max(fit$z), length=100))),
                                     bandwidth = rep(bandwidth, 2), degree.pol = 0)$predicted
        }
        
        belowa.comire <- sapply(index, tailp.mcmc2, fit$mcmc$nu0, fit$mcmc$th0, fit$mcmc$tau0,
                                fit$mcmc$th1, fit$mcmc$tau1, fit$mcmc$w, fit$mcmc$ga, phi(x))
        
        belowa.true <- NonpModelCheck::localpoly.reg( X=as.matrix(cbind(x,fit$z)), Y=y<a, 
                                      points=as.matrix(cbind(seq(0, max.x, length=100), 
                                                             seq(min(fit$z), max(fit$z), length=100))), 
                                      bandwidth=rep(bandwidth,2), degree.pol = 0)$predicted
       
        ppc.data <- data.frame(cbind(x=seq(0,max.x, length=100),
                                     Fx=c(c(belowa.comire[,1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery])), 
                                     repl=c(rep(1:length(1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery),each=100))) )
        
        ppc <-  ggplot2::ggplot(ppc.data, ggplot2::aes(x=.data$x, y=.data$Fx)) + 
          ggplot2::geom_line(alpha=0.25, ggplot2::aes(group=factor(.data$repl)), col="grey") +
          ggplot2::geom_line(data=data.frame(x=seq(0,max.x, length=100), y=belowa.true), ggplot2::aes(x=.data$x, y=.data$y), col=1)+
          ggplot2::labs(x = "", y = expression(F(a *"| x,z")*" | data")) +
          ggplot2::geom_point(data=data.frame(x, zero=rep(0,length(x))), ggplot2::aes(.data$x, .data$zero), alpha=1, cex=.5, pch="|", na.rm=TRUE) 
        
        ppc + ggplot2::coord_cartesian(ylim=c(0,1), xlim=xlim)

      }
      else{
        index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
        knots <- seq(min(x)+1, max.x, length=J-3)
        phi <- function(x) splines2::iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
        res <- rep(NA, length(index))
        
        tailp.mcmc3 <- function(i, pi0, th0, tau0, th1, tau1, w, ga, phi_x)
        {
          
          beta_x <-  as.double(phi_x %*% w[i,])
          res <-  beta_x * stats::pnorm(a, th1[i]+crossprod(t(fit$z), ga[i,]), 1/sqrt(tau1[i]))
          
          for(h in 1:ncol(pi0))
          {
            res <- res + (1-beta_x) * pi0[i,h] * stats::pnorm(a, th0[i,h]+crossprod(t(fit$z), ga[i,]) , 1/sqrt(tau0[i,h]))     
          }
          
          below <- stats::rbinom(length(res), 1, res) 
          
          p <- ncol(fit$z)
          z.points<- matrix(NA, nrow=100, ncol=p)
          for(j in 1:p){
            z.points[,j]<- seq(min(fit$z[,j]), max(fit$z[,j]), length=100)
          }
          smoothed <- NonpModelCheck::localpoly.reg( X=as.matrix(cbind(x,fit$z)), Y=below, 
                                     points=as.matrix(cbind(seq(0, max.x, length=100), z.points)),
                                     bandwidth = rep(bandwidth, p+1), degree.pol = 0)$predicted
        }
        
        belowa.comire <- sapply(index, tailp.mcmc3, fit$mcmc$nu0, fit$mcmc$th0, fit$mcmc$tau0,
                      fit$mcmc$th1, fit$mcmc$tau1, fit$mcmc$w, fit$mcmc$ga, phi(x))
        points <- matrix(NA, ncol=ncol(fit$z), nrow=100)
        for(j in 1:ncol(fit$z)) points[,j]<- seq(min(fit$z[,j]), max(fit$z[,j]), length=100)
        belowa.true <- NonpModelCheck::localpoly.reg(as.matrix(cbind(x,fit$z)), Y=y<a, 
                                       points=as.matrix(cbind(seq(0, max.x, length=100), points)), 
                                       bandwidth=rep(bandwidth,(ncol(fit$z))+1), degree.pol = 0)$predicted
        
        ppc.data <- data.frame(cbind(x=seq(0,max.x, length=100),
                                     Fx=c(c(belowa.comire[,1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery])), 
                                     repl=c(rep(1:length(1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery),each=100))))
        ppc <-  ggplot2::ggplot(ppc.data, ggplot2::aes(x=.data$x, y=.data$Fx)) + 
          ggplot2::geom_line(alpha=0.25, ggplot2::aes(group=factor(.data$repl)), col="grey") +
          ggplot2::labs(x = "", y = expression(F(a*";x,z")*" | data")) +
          ggplot2::geom_point(data=data.frame(x, zero=rep(0,length(x))), 
                              ggplot2::aes(.data$x, .data$zero), alpha=1, cex=.5, pch="|", na.rm=TRUE) +
          ggplot2::geom_line(data=data.frame(x=seq(0,max.x, length=100), y=belowa.true), ggplot2::aes(x=.data$x, y=.data$y), col=1)
        ppc + ggplot2::coord_cartesian(ylim=c(0,1), xlim=xlim) 
        
      }
    }
  }
  else{
    index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
    x.grid <- seq(0, max(x), length=100)
    pi_x_prova <- matrix(NA, nrow=100, ncol = length(index))
    ite=1
    for(i in index){
      pi_x_prova[,ite] <- (1-fit$mcmc$beta[i,])*fit$mcmc$pi0[i]+fit$mcmc$beta[i,]
      ite=ite+1
    }
    ppc.data <- data.frame(cbind(x=seq(0,max(x), length=100), 
                                 pi_x=c(c(pi_x_prova[,1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery])), 
                                 repl=c(rep(1:length(1:(((mcmc$nrep)/mcmc$thin)/oneevery)*oneevery),each=100))))
    
    true <- KernSmooth::locpoly(x=x, y=y, degree=1, bandwidth = bandwidth, gridsize = 500, range.x = c(0, max(x)))
  
    ppc <- ggplot2::ggplot(ppc.data, ggplot2::aes(x=.data$x, y=.data$pi_x)) + 
      ggplot2::geom_line(alpha=0.25, ggplot2::aes(group=factor(.data$repl)), col="grey")+
      ggplot2::geom_line(data=data.frame(x=seq(0,max(x), length=500), y=true$y), ggplot2::aes(x=.data$x, y=.data$y), col=1)+
      ggplot2::labs(x = "", y = expression(pi[x])) +
      ggplot2::geom_point(data=data.frame(x, zero=rep(0,length(x))), ggplot2::aes(.data$x, .data$zero), alpha=1, cex=.5, pch="|", na.rm=TRUE) 
    
    ppc + ggplot2::coord_cartesian(ylim=c(0,1), xlim=xlim)

  }
}

# -----------------------------------------------------------------------
# PLOT 2: posterior mean density
# -----------------------------------------------------------------------

#' @rdname fit.pdf.mcmc
#' 
#' @importFrom ggplot2 theme_set theme_bw ggplot aes geom_line geom_ribbon ylim xlim labs theme element_text
#' @importFrom splines2 iSpline
#' @importFrom stats dnorm quantile var
#' @importFrom grid unit
#' @export 
#' 
#' @usage fit.pdf.mcmc(y, x, fit, mcmc, J=10, H = 10, alpha = 0.05, 
#' max.x = max(x), x.val, y.grid = NULL, xlim = c(0, max(x)), 
#' ylim = c(0, 1), xlab = NULL)
#' 
#' @title Posterior mean density plot for dose intervals
#' 
#' @description Pointwise posterior mean (continuous blue lines), and credible bands (shaded blue areas) for f (y | x, z) 
#' calculated in \code{x.val} under the the model fitted in \code{fit}.
#' 
#' @param y optional numeric vector for the response used in \code{comire.gibbs}. If \code{y} is missing, \code{y.grid} must be provided.
#' @param x numeric vector for the covariate relative to the dose of exposure used in \code{comire.gibbs}.
#' @param fit the output of \code{comire.gibbs} opportunely trasformed in \code{classCoMiRe} class.
#' @param mcmc a list giving the MCMC parameters.
#' @param J parameter controlling the number of elements of the I-spline basis
#' @param H total number of components in the mixture at \eqn{x_0}.
#' @param alpha level of the credible bands.
#' @param max.x maximum value allowed for x.
#' @param x.val central points of each dose interval to be used in the posterior estimation of the probability density function.
#' @param y.grid optional numerical vector giving the actual values of the grid for y for plotting the posterior mean density. If \code{y.grid} is not provided, standard grids are automatically used.
#' @param xlim,ylim numeric vectors of length 2, giving the x and y coordinates ranges for the plot.
#' @param xlab the title of the x axis.
#' 
#' @author Antonio Canale, Arianna Falcioni
#' 
#' @examples{
#' data(CPP)
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' premature <- as.numeric(gestage<=37)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## too few iterations to be meaningful. see below for safer and more comprehensive results
#' 
#' mcmc <- list(nrep=10, nb=2, thin=1, ndisplay=4) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit.dummy <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=1, max.x=180)
#'                      
#' fit.pdf.mcmc(y = gestage, x = dde, fit = fit.dummy, mcmc = mcmc, J = 10, H = 10, 
#'                          alpha = 0.05, max.x = max(dde), x.val = 125, 
#'                          xlim = c(25,48), ylim = c(0,0.25),
#'                          xlab = "Gest. age. for DDE = 125")
#'                          
#' \donttest{
#' ## safer procedure with more iterations (it may take some time)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## Fit the model for continuous y 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit1 <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#'           
#' fit.pdf.mcmc(y = gestage, x = dde, fit = fit1, mcmc = mcmc, J = 10, H = 10, 
#'                          alpha = 0.05, max.x = max(dde), x.val = 125, 
#'                          xlim = c(25,48), ylim = c(0,0.25),
#'                          xlab = "Gest. age. for DDE = 125")
#'                          
#' }
#' }

fit.pdf.mcmc <- function(y, x, fit, mcmc, J=10, H=10, alpha=0.05, max.x=max(x), x.val, y.grid=NULL, xlim=c(0, max(x)), ylim=c(0,1), xlab=NULL)
{
  
  if(fit$bin) stop("No function fit.pdf.mcmc defined for family = 'binary'
         (no density function for a binary variable!)\n")
  y.grid <- if(is.null(y.grid)) seq(min(y)-sqrt(stats::var(y)), max(y) + sqrt(stats::var(y)), length = 100) else y.grid 
  ggplot2::theme_set(ggplot2::theme_bw())
  index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
  knots <- seq(min(x)+1, max.x, length=J-3)
  phi <- function(x) splines2::iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
  
  if(is.null(fit$z)){
      # comire senza confunders
      dens.mcmc1 <- function(i, nu0, mu0, tau0, mu1, tau1, w, phi_x, y.grid)
      {
        beta_x <-  as.double(phi_x %*% w[i,])
        res <-  beta_x * stats::dnorm(y.grid, mu1[i], 1/sqrt(tau1[i]))
        for(h in 1:ncol(nu0))
        {
          res <- res + (1-beta_x)*nu0[i,h] * stats::dnorm(y.grid, mu0[i,h], 1/sqrt(tau0[i,h]))     
        }
        res
      }
      densapply <- function(x){sapply(index, dens.mcmc1, nu0=fit$mcmc$nu0, mu0=fit$mcmc$th0, tau0=fit$mcmc$tau0,
                                mu1=fit$mcmc$th1, tau1=fit$mcmc$tau1, w=fit$mcmc$w, phi_x=phi(x), y.grid=y.grid)}
  }
  #con confunders
    else {
      #comire con un solo confunder
      if(fit$univariate){
        dens.mcmc2 <- function(i, nu0, th0, tau0, th1, tau1, w, phi_x, ga, y.grid)
        {
          beta_x <-  as.double(phi_x %*% w[i,])
          res <-  beta_x * stats::dnorm(y.grid, th1[i]+ga[i]*fit$z.val, 1/sqrt(tau1[i]) )
          for(h in 1:ncol(nu0))
          {
            res <- res + (1-beta_x)*nu0[i,h] * stats::dnorm(y.grid, th0[i,h]+ga[i]*fit$z.val , 1/sqrt(tau0[i,h]))     
          }
          res
        }
        densapply <- function(x){sapply(index, dens.mcmc2, nu0=fit$mcmc$nu0, th0=fit$mcmc$th0, tau0=fit$mcmc$tau0,
                                  th1=fit$mcmc$th1, tau1=fit$mcmc$tau1, w=fit$mcmc$w, phi_x=phi(x), 
                                  ga = fit$mcmc$ga, y.grid=y.grid)}
        }
      else{
        #comire con piu' cunfunders
        dens.mcmc3 <- function(i, nu0, th0, tau0, th1, tau1, w, phi_x, ga, y.grid)
        {
          beta_x <-  as.double(phi_x %*% w[i,])
          res <-  beta_x * stats::dnorm(y.grid, th1[i]+crossprod(t(fit$z.val), ga[i,]), 1/sqrt(tau1[i]) )
          for(h in 1:ncol(nu0))
          {
            res <- res + (1-beta_x)*nu0[i,h] * stats::dnorm(y.grid, th0[i,h]+crossprod(t(fit$z.val), ga[i,]) , 1/sqrt(tau0[i,h]))     
          }
          res
        }
        densapply <- function(x){sapply(index, dens.mcmc3, nu0=fit$mcmc$nu0, th0=fit$mcmc$th0, tau0=fit$mcmc$tau0,
                                  th1=fit$mcmc$th1, tau1=fit$mcmc$tau1, w=fit$mcmc$w, phi_x=phi(x), 
                                  ga = fit$mcmc$ga, y.grid=y.grid)}
      }
    }
  # plot the length(x.val) conditional densities
  all.pdf <- list()
  for(j in 1:length(x.val))
  {
    resj <- densapply(x.val[j])
    pdf_fit <- cbind( rowMeans(resj), t(apply(resj, 1, stats::quantile, probs=c(alpha/2,1-alpha/2))))
    data <- data.frame(pdf_fit, y.grid)
    names(data)[1:3] <- c("mean","low","upp")
    
    pdf.j <- ggplot2::ggplot(data) +  
      ggplot2::geom_line(ggplot2::aes(x=y.grid, y=mean), col="blue", na.rm=TRUE) +
      ggplot2::geom_ribbon(ggplot2::aes(ymax=.data$upp, ymin=.data$low, x=y.grid), fill=4, alpha=.1) + 
      ggplot2::ylim(ylim) +  ggplot2::xlim(xlim)+  ggplot2::labs(x=xlab[j], y="") + 
      ggplot2::theme(plot.margin=grid::unit(c(1,1,1,1),"lines"), axis.title=ggplot2::element_text(size=6))
    all.pdf[[j]] <- pdf.j
  }
  all.pdf
  }



# -----------------------------------------------------------------------
# PLOT 3: beta plot
# -----------------------------------------------------------------------

#' @rdname betaplot
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_bw ylim xlim geom_point
#' @importFrom rlang .data
#' @export 
#' 
#' @usage betaplot(x, fit, x.grid = NULL, xlim = c(0, max(x)), xlab = "x")
#' 
#' @title \eqn{\beta(x)} plot
#' 
#' @description Posterior mean (continuous lines) and pointwise credible bands (shaded areas) for \eqn{\beta(x)}.
#' 
#' @param x numeric vector for the covariate relative to the dose of exposure used in \code{comire.gibbs}.
#' @param fit the output of \code{comire.gibbs} opportunely trasformed in \code{classCoMiRe} class.
#' @param x.grid optional numerical vector giving the actual values of the grid for x for plotting \eqn{\beta(x)}. 
#' If \code{x.gird} is not provided, standard grids are automatically used.
#' @param xlim numeric vectors of length 2, giving the x coordinates ranges for the plot.
#' @param xlab the title of the x axis.
#' 
#' @author Antonio Canale
#'  
#' @examples{
#' data(CPP)
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' premature <- as.numeric(gestage<=37)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## too few iterations to be meaningful. see below for safer and more comprehensive results
#' 
#' mcmc <- list(nrep=10, nb=2, thin=1, ndisplay=4) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit.dummy <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=1, max.x=180)
#'                      
#' betaplot(x=dde, fit=fit.dummy, x.grid=seq(0,180, length=100), xlim=c(0,150))
#' 
#' \donttest{
#' ## safer procedure with more iterations (it may take some time)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## Fit the model for continuous y 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit1 <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#' 
#'                          
#' betaplot(x=dde, fit=fit1, x.grid=seq(0,180, length=100), xlim=c(0,150))
#' 
#' }
#' }
betaplot <- function(x, fit, x.grid=NULL, xlim=c(0,max(x)), xlab="x")
{
  x.grid <- if(is.null(x.grid)) seq(0, max(x), length = 100) else x.grid 
  beta.data <- data.frame(beta=fit$post.means$beta,
                          low=fit$ci$beta[1,], upp=fit$ci$beta[2,],
                          x.grid=x.grid)
  betaplot <- ggplot2::ggplot(beta.data, ggplot2::aes(x.grid,beta)) +
    ggplot2::geom_line(lty=1, col=4, na.rm=TRUE) + 
    ggplot2::geom_ribbon(ggplot2::aes(ymax=.data$upp, ymin=.data$low), fill=4, alpha=.1) +
    ggplot2::labs(y=expression(beta(x)), x=xlab) + 
    ggplot2::theme_bw() + ggplot2::ylim(c(0,1)) +  ggplot2::xlim(xlim) +
    ggplot2::geom_point(data=data.frame(x, zero=rep(0,length(x))), ggplot2::aes(.data$x, .data$zero), alpha=1, cex=.5, pch="|", na.rm=TRUE) 
  betaplot
}


# -----------------------------------------------------------------------
# PLOT 4: additional risk
# -----------------------------------------------------------------------

#' @rdname add.risk
#' 
#' @importFrom stats approxfun integrate quantile var
#' @export  
#' 
#' @usage add.risk(y, x, fit, mcmc, a, alpha=0.05, 
#' x.grid=NULL, y.grid=NULL)
#' 
#' @title Additional risk function 
#' 
#' @description Additional risk function estimated from the object \code{fit}
#' 
#' @param y optional numeric vector for the response used in \code{comire.gibbs}. If \code{y} is missing, \code{y.grid} must be provided.
#' @param x numeric vector for the covariate relative to the dose of exposure used in \code{comire.gibbs}.
#' @param fit the output of \code{comire.gibbs}. an object of the class \code{classCoMiRe}.
#' @param mcmc a list giving the MCMC parameters.
#' @param a threshold of clinical interest for the response variable
#' @param alpha level of the credible bands.
#' @param x.grid optional numerical vector giving the actual values of the grid for x for plotting the additional risk function. If \eqn{x.gird} is not provided, standard grids are automatically used.
#' @param y.grid optional numerical vector giving the actual values of the grid for y for plotting the additional risk function. If \eqn{y.gird} is not provided, standard grids are automatically used.
#' 
#' @return A list of arguments for generating posterior output. It contains:
#' \itemize{
#' \item{\code{mcmc.risk}}{ a matrix containing in the lines the MCMC chains, after thinning, of the additional risk function over \code{x.grid}, in the columns. }
#' \item{\code{summary.risk}}{ a data frame with four variables: the posterior means of the additional risk function over \code{x.grid}, the respective \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles, and \code{x.grid}.}
#' }
#' 
#' @author Antonio Canale, Arianna Falcioni
#' 
#' @examples{
#' data(CPP)
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' premature <- as.numeric(gestage<=37)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## too few iterations to be meaningful. see below for safer and more comprehensive results
#' 
#' mcmc <- list(nrep=10, nb=2, thin=1, ndisplay=4) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit.dummy <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=1, max.x=180)
#'                      
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit.dummy, mcmc = mcmc, 
#'     a = 37, x.grid = seq(0, max(dde), length = 100))
#' riskplot(risk.data$summary.risk, xlab="DDE", x = dde, xlim = c(0,150))
#'                     
#' \donttest{
#' ## safer procedure with more iterations (it may take some time)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## Fit the model for continuous y 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit1 <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#'  
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit1, mcmc = mcmc, 
#' a = 37, x.grid = seq(0, max(dde), length = 100))
#' riskplot(risk.data$summary.risk, xlab="DDE", x = dde, xlim = c(0,150))
#' 
#' }
#' }
add.risk <- function(y, x, fit, mcmc, a, alpha=0.05, x.grid=NULL, y.grid=NULL)
{
  y.grid <- if(is.null(y.grid)) seq(min(y)-sqrt(stats::var(y)), max(y) + sqrt(stats::var(y)), length = 100) else y.grid 
  x.grid <- if(is.null(x.grid)) seq(0, max(x), length = 100) else x.grid 
  if(!fit$bin){
    if(is.null(fit$z)){
      index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
      r_a <- matrix(NA, length(index), length(x.grid))
      low.int <- floor(min(y))
      ite <- 1
      for(i in index)
      {
        f0 <- stats::approxfun(x=y.grid, y=fit$mcmc$f0[i,])
        f1 <- stats::approxfun(x=y.grid, y=fit$mcmc$f1[i,])
        r_a[ite,] <- fit$mcmc$beta[i,]*(stats::integrate(function(x)f1(x), low.int, a)$value-
                                          stats::integrate(function(x)f0(x), low.int, a)$value )
        
        ite <- ite + 1
      }
      risk.data <- data.frame(apply(r_a,2,mean), 
                              apply(r_a,2, stats::quantile, probs=alpha/2), 
                              apply(r_a,2, stats::quantile, probs=1-alpha/2), 
                              x.grid)
      colnames(risk.data) <- c("risk","low","upp", "x")
      list(mcmc.risk=r_a, summary.risk= risk.data)
    }
    else {
        index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
        r_a <- matrix(NA, length(index), length(x.grid))
        low.int <- floor(min(y))
        ite <- 1
        for(i in index)
        {
          f0_a <- stats::approxfun(x=y.grid, y=fit$f0[i,])
          f1_a <- stats::approxfun(x=y.grid, y=fit$f1[i,])
          r_a[ite,] <- fit$mcmc$beta[i,]*(stats::integrate(function(x) f1_a(x), low.int, a)$value-
                                            stats::integrate(function(x) f0_a(x), low.int, a)$value )
          
          ite <- ite + 1
        }
        risk.data <- data.frame(apply(r_a,2,mean), 
                                apply(r_a,2, stats::quantile, probs=alpha/2), 
                                apply(r_a,2, stats::quantile, probs=1-alpha/2), 
                                x.grid)
        colnames(risk.data) <- c("risk","low","upp", "x")
        list(mcmc.risk=r_a, summary.risk= risk.data)
    }
  }
}

# -----------------------------------------------------------------------
# additional risk plot
# -----------------------------------------------------------------------
#' @rdname riskplot
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_bw theme ylim xlim geom_point
#' @importFrom grid unit
#' @export 
#' 
#' @usage riskplot(risk.data, xlab = NULL, x = NULL, ylim=c(0,1), xlim=c(0, max(x)))
#' 
#' @title Additional risk function plot
#' 
#' @description Posterior mean (continuous lines) and pointwise credible bands (shaded areas) for \eqn{Ra(x, a)}.
#' 
#' @param risk.data output of \code{add.risk} function.
#' @param xlab the title of the x axis.
#' @param x numeric vector for the covariate relative to the dose of exposure used in \code{comire.gibbs}.
#' @param xlim,ylim numeric vectors of length 2, giving the x and y coordinates ranges for the plot.
#' 
#' @author Antonio Canale
#' 
#' @examples{
#' data(CPP)
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' premature <- as.numeric(gestage<=37)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## too few iterations to be meaningful. see below for safer and more comprehensive results
#' 
#' mcmc <- list(nrep=10, nb=2, thin=1, ndisplay=4) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit.dummy <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=1, max.x=180)
#'                      
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit.dummy, mcmc = mcmc, 
#'     a = 37, x.grid = seq(0, max(dde), length = 100))
#' riskplot(risk.data$summary.risk, xlab="DDE", x = dde, xlim = c(0,150))
#' 
#' \donttest{
#' ## safer procedure with more iterations (it may take some time)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## Fit the model for continuous y 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#'  
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit, mcmc = mcmc, 
#' a = 37, x.grid = seq(0, max(dde), length = 100))
#' riskplot(risk.data$summary.risk, xlab="DDE", 
#' x = dde, xlim = c(0,150))
#' 
#' }
#' }
riskplot <- function(risk.data, xlab=NULL, x=NULL, ylim=c(0,1), xlim=c(0, max(x)))
{
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  out <- ggplot2::ggplot(risk.data, ggplot2::aes(.data$x,.data$risk)) + ggplot2::geom_line(lty=1, col=4, na.rm=TRUE) + 
    ggplot2::geom_ribbon(ggplot2::aes(ymax=.data$upp, ymin=.data$low), fill=4, alpha=.1) +
    ggplot2::labs(y=expression(R[A](x, a)), x=xlab)+ ggplot2::theme_bw() + 
    ggplot2::theme(plot.margin=grid::unit(c(1,1,1,1),"lines")) + ggplot2::ylim(ylim) + ggplot2::xlim(xlim)
  if(is.null(x)) return(out)
  else{
    onlyx <- data.frame(x, zero=rep(0,length(x)))
    out <- out + ggplot2::geom_point(data=onlyx, ggplot2::aes(x, .data$zero), alpha=1, cex=.5, pch="|", na.rm=TRUE)
  }
  return(out)
}

# -----------------------------------------------------------------------
# benchmark dose
# -----------------------------------------------------------------------
#' @rdname BMD
#' 
#' @importFrom stats splinefun uniroot quantile
#' @export 
#' 
#' @usage BMD(level, risk, x, alpha=0.05)
#' 
#' @title Benchmark dose 
#' 
#' @description Benchmark dose associated to a particular risk
#' 
#' @param level dose level of interest.
#' @param risk \code{summary.risk$mcmc.risk} from the output of \code{add.risk} function.
#' @param x numeric vector for the covariate relative to the dose of exposure used in \code{comire.gibbs}.
#' @param alpha level of the credible bands.
#' 
#' @return A dataframe containing as variables:
#' \itemize{
#' \item{\code{q}}{ the dose level of interest.}
#' \item{\code{BMD}}{ the benchmark dose.}
#' \item{\code{low}}{ lower credible limit.}
#' \item{\code{upp}}{ upper credible limit.}
#' \item{\code{BMDL}}{ a more conservative benchmark dose.}
#' }
#' 
#' @author Antonio Canale
#' 
#' @examples{
#' data(CPP)
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' premature <- as.numeric(gestage<=37)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## too few iterations to be meaningful. see below for safer and more comprehensive results
#' 
#' mcmc <- list(nrep=10, nb=2, thin=1, ndisplay=4) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit.dummy <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=1, max.x=180)
#'                      
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit.dummy, mcmc = mcmc, 
#'     a = 37, x.grid = seq(0, max(dde), length = 100))
#' bmd.data <- BMD(seq(0,.20, length=50), risk.data$mcmc.risk, 
#' x=seq(0,max(dde), length=100), alpha=0.05)
#' bmd.plot(bmd.data)       
#'                      
#' \donttest{
#' ## safer procedure with more iterations (it may take some time)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## Fit the model for continuous y 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#' 
#' 
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit, mcmc = mcmc,
#' a = 37, x.grid = seq(0, max(dde), length = 100))
#' bmd.data <- BMD(seq(0,.20, length=50), risk.data$mcmc.risk, 
#' x=seq(0,max(dde), length=100), alpha=0.05)
#' bmd.plot(bmd.data)       
#' 
#' }
#' }

BMD <- function(level, risk, x, alpha=0.05)
{
  bmd <- function(r,q,range=c(0,max(x))) 
  {
    ris <- stats::splinefun(y=r,x=x)
    if(ris(range[2])-q<0) return(NA)
    else
      stats::uniroot(function(x) ris(x)-q,range)$root
  }
  bmd.apply <- function(q) apply(risk, 1, bmd, q=q, range=c(0,180))
  if(length(level)==1){
    bmd.data <- bmd.apply(level)
    return(c(mean(bmd.data), stats::quantile(bmd.data, probs=c(alpha/2,1-alpha/2,alpha))))
  }
  else
  {
    bmd.mcmc.comire <-  sapply(level, bmd.apply)
    bmd.data <- data.frame(level,
                           colMeans(bmd.mcmc.comire, na.rm=TRUE), 
                           t(apply(bmd.mcmc.comire,2, stats::quantile, prob=c(alpha/2,1-alpha/2), na.rm=TRUE)),
                           apply(bmd.mcmc.comire,2, stats::quantile, prob=alpha, na.rm=TRUE))
    colnames(bmd.data) <- c("q", "BMD","low","upp", "BMDL")
    return(bmd.data)
  }
}

# -----------------------------------------------------------------------
# bmd plot
# -----------------------------------------------------------------------
#' @rdname bmd.plot
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_bw theme
#' @importFrom rlang .data
#' @importFrom grid unit
#' @export 
#' 
#' @title Benchmark dose plot
#' 
#' 
#' @description Posterior mean (continuous lines) and pointwise credible bands (shaded areas) for the benchmark dose in function of the increase in risk.
#' 
#' @param bmd.data output of \code{BMD} function.
#'
#' @author Antonio Canale
#' 
#' @examples{
#' data(CPP)
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' premature <- as.numeric(gestage<=37)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## too few iterations to be meaningful. see below for safer and more comprehensive results
#' 
#' mcmc <- list(nrep=10, nb=2, thin=1, ndisplay=4) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit.dummy <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=1, max.x=180)
#'                      
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit.dummy, mcmc = mcmc, 
#'     a = 37, x.grid = seq(0, max(dde), length = 100))
#' bmd.data <- BMD(seq(0,.20, length=50), risk.data$mcmc.risk, 
#' x=seq(0,max(dde), length=100), alpha=0.05)
#' bmd.plot(bmd.data)       
#' 
#' \donttest{
#' ## safer procedure with more iterations (it may take some time)
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' ## Fit the model for continuous y 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#'               
#' fit <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#'                      
#' risk.data <- add.risk(y = gestage, x = dde, fit = fit, mcmc = mcmc, 
#' a = 37, x.grid = seq(0, max(dde), length = 100))
#' bmd.data <- BMD(seq(0,.20, length=50), risk.data$mcmc.risk, 
#' x=seq(0,max(dde), length=100), alpha=0.05)
#' bmd.plot(bmd.data)       
#' 
#' }
#' }

bmd.plot <- function(bmd.data)
{
  ggplot2::ggplot(bmd.data, ggplot2::aes(q,BMD)) + 
    ggplot2::geom_line(ggplot2::aes(q, BMD), lty=1, col=4) + ggplot2::geom_ribbon(ggplot2::aes(ymax=.data$upp, ymin=.data$low), fill=4, alpha=.1) +
    ggplot2::labs(y=expression(BMD[q]), x="q")+ ggplot2::theme_bw() + 
    ggplot2::theme(plot.margin=grid::unit(c(1,1,1,1),"lines"))  
}


# -----------------------------------------------------------------------
# classCoMiRe plot method
# -----------------------------------------------------------------------
#' CoMiRe plot 
#' 
#' 
#' @description An S3 plot method for an object of \code{classCoMiRe} class.
#' 
#' @details  The output is a list of 
#' \code{ggplot2} plots containing the result of the \code{betaplot} function and, if the threshold \code{a} is provided, of 
#' \code{post.pred.check}, \code{riskplot}, \code{bmd.plot}.
#'  
#' @param x the output of \code{comire.gibbs}, an object of the \code{classCoMiRe} class.
#' @param y numeric vector for the response used in \code{comire.gibbs}.
#' @param xobs numeric vector for the covariate relative to the dose of exposure used in \code{comire.gibbs}.
#' @param mcmc a list giving the MCMC parameters.
#' @param J parameter controlling the number of elements of the I-spline basis
#' @param H total number of components in the mixture at \eqn{x_0}.
#' @param a optional threshold of clinical interest for the response variable.
#' @param max.x maximum value allowed for x.
#' @param bandwidth the kernel bandwidth smoothing parameter for the \code{post.pred.check} plot.
#' @param x.grid optional numerical vector giving the actual values of the grid for x for plotting the additional risk function. If \eqn{x.gird} is not provided, standard grids are automatically used.
#' @param xlim,ylim numeric vectors of length 2, giving the x and y coordinates ranges for the plot.
#' @param xlab the title of the x axis.
#' @param alpha level of the credible bands, default 0.05
#' @param risk if \code{TRUE} the additional risk plot via \code{riskplot} is computed.
#' @param bmd if \code{TRUE} the benchmark dose plot via \code{bmd.plot} is computed.
#' @param level if \code{bmd=TRUE}, dose levels of interest for BMD plot.
#' @param oneevery integer number representing how many MCMC draws to plot in the posterior predictive check. It draws one sample every \code{oneevery}.
#' @param ... additional arguments to be passed.
#' 
#' @rdname plot.classCoMiRe
#' 
#' @return If \code{a=NULL} returns only \code{betaplot} otherwise, if \code{risk=FALSE} and 
#' \code{bmd=FALSE} returns a list containing \code{betaplot} (which is automatically plotted) 
#' and \code{post.pred.check} plot. Finally, if \code{a} is provided, \code{risk=TRUE} and \code{bmd=TRUE}
#' returns a list with \code{betaplot}, \code{post.pred.check}, \code{riskplot} and \code{bmd.plot}.
#' 
#' @author Antonio Canale, Arianna Falcioni
#' @export 


plot.classCoMiRe<- function(x, y, xobs, mcmc, J = 10 , H=10, 
                          a=NULL, max.x=max(xobs), bandwidth=20, x.grid=NULL,
                          xlim=c(0, max(xobs)), ylim=c(0,1), xlab="x", alpha=0.05, 
                          risk=TRUE, bmd=TRUE, level, oneevery=20, ...){
  fit <- x
  oneevery <- ifelse(is.null(oneevery), 20, oneevery)
  if(fit$bin){
    index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
    plot1 <- post.pred.check(y=y, x=xobs, fit=fit, bandwidth = bandwidth, 
                             mcmc=mcmc, J=J,  H=H, max.x=max.x, xlim=xlim, oneevery=oneevery)
    
    plot2 <- betaplot(x=xobs, fit=fit, x.grid=x.grid, xlim=xlim)
    
    bmd.data <- BMD(level, fit$mcmc$beta[index,], x=seq(0,max(xobs), length=100), 
                    alpha=alpha)
    plot3 <- bmd.plot(bmd.data) 
    
    out<- list(betaplot=plot2, post.pred.check=plot1, bmdplot=plot3)
    
    plot2
    return(out)
  }
  else{
  if(is.null(a)){
    plot1<- betaplot(x=xobs, fit=fit, x.grid=x.grid, xlim=xlim, xlab=xlab)
  }
  else if(!is.null(a) & !risk & !bmd){
    plot1<- betaplot(x=xobs, fit=fit, x.grid=x.grid, xlim=xlim, xlab=xlab)
    plot2<- post.pred.check(y=y, x=xobs, fit=fit, mcmc=mcmc, J=J, H=H, a=a, max.x=max.x, 
                            xlim=xlim, bandwidth=bandwidth, oneevery = oneevery)
    out<- list(betaplot=plot1, post.pred.check=plot2)
    plot1
    return(out)
  }
  else if(!is.null(a) & risk & bmd){
    plot1<- betaplot(x=xobs, fit=fit, x.grid=x.grid, xlim=xlim, xlab=xlab)
    plot2<- post.pred.check(y=y, x=xobs, fit=fit, mcmc=mcmc, J=J, H=H, a=a, max.x=max.x, 
                          xlim=xlim, bandwidth=bandwidth,  oneevery = oneevery)
    risk.data <- add.risk(y=y, x=xobs, fit=fit, mcmc=mcmc, a=a, x.grid=x.grid, alpha=alpha)
    plot3<- riskplot(risk.data$summary.risk, xlab=xlab, x=xobs, ylim=ylim, xlim=xlim)
    
    bmd.data <- BMD(level=level, risk.data$mcmc.risk, 
                    x=seq(0,max(xobs), length=100), 
                    alpha=alpha)
    plot4<- bmd.plot(bmd.data=bmd.data)

    out<- list(betaplot=plot1, post.pred.check=plot2, riskplot=plot3, bmdplot=plot4)
    
    plot1
    return(out)
  }
  }
}

# -----------------------------------------------------------------------
# classCoMiRe print method
# -----------------------------------------------------------------------

#' @description The \code{print.classCoMiRe} method prints the type of a \code{classCoMiRe} object.
#' @title CoMiRe print 
#' @param x an object of class \code{classCoMiRe};
#' @param ... additional arguments.
#' @rdname print.classCoMiRe
#' @export
#' 
#' @author Antonio Canale, Arianna Falcioni
#' 


print.classCoMiRe <- function(x, ...) {
  fit <- x
  if(!fit$bin){
    if(is.null(fit$z)){
      cat("CoMiRe model fit via Gibbs Sampler\n")
      cat("Family: Continuous")
      }
    else {
      if(fit$univariate){
        cat("CoMiRe model fit via Gibbs Sampler\n")
        cat("Family: Continuous")
      }
      else{
        cat("CoMiRe model fit via Gibbs Sampler\n")
        cat("Family: Continuous")
      }
    }
  }
  else{
    cat("CoMiRe model fit via Gibbs Sampler\n")
    cat("Family: Binary")
  }

}
  

# -----------------------------------------------------------------------
# classCoMiRe summary method
# -----------------------------------------------------------------------

#'
#' @description The \code{summary.classCoMiRe} method provides summary information on \code{classCoMiRe} objects.
#' @title CoMiRe summary 
#' 
#' @param object an object of class \code{classCoMiRe};
#' @param ... additional arguments
#'
#' @rdname summary.classCoMiRe
#' @export
#' 
#' @author Antonio Canale Arianna Falcioni


summary.classCoMiRe <- function(object, ...) {
  fit <- object
    if(!fit$bin){
      if(is.null(fit$z)){
        cat("CoMiRe model fit via Gibbs Sampler\n")
        cat("Family: Continuous\n")
        cat(paste(c("Formula: ", fit$call), collapse =""))
        cat("\nPosterior approximation based on", fit$nrep ,"iterations")
            }
      else {
        if(fit$univariate){
          cat("CoMiRe model fit via Gibbs Sampler\n")
          cat("Family: Continuous\n")
          cat(paste(c("Formula: ", fit$call), collapse =""))
          cat("\nPosterior approximation based on", fit$nrep ,"iterations")
        }
        else{
          cat("CoMiRe model fit via Gibbs Sampler\n")
          cat("Family: Continuous\n")
          cat(paste(c("Formula: ", fit$call), collapse =""))
          cat("\nPosterior approximation based on", fit$nrep ,"iterations")
        }
      }
    }
    else{
      cat("CoMiRe model fit via Gibbs Sampler\n")
      cat("Family: Binary\n")
      cat(paste(c("Formula: ", fit$call), collapse =""))
      cat("\nPosterior approximation based on", fit$nrep ,"iterations")
    }
}
    



  
