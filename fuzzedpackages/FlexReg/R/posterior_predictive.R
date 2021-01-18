#' Posterior Predictive
#'
#' The function takes an object of class \code{`flexreg`} and generates values from the posterior predictive distribution.
#'
#' @param model an object of class \code{`flexreg`}
#' @param newdata an optional data frame containing variables with which to predict. If omitted, the fitted values are used.
#'
#' @details The function generates values from the posterior predictive distribution, which is the distribution of a  future outcome given the observed data.
#' @return An object of class \code{`flexreg_postpred`} containing a matrix with the simulated posterior predictions. Each column refers to a statistical unit to predict.
#'
#' @examples
#' data("Reading")
#' dataset <- Reading
#' FB <- flexreg(accuracy ~ iq, dataset, n.iter=1000)
#' pp <- posterior_predict(FB)
#' plot(pp)
#'
#' @import rstan stats Formula
#'
#' @references{Gelman, A.; Carlin, J. B.; Stern, H. S. and Rubin, D. B. (2014), Bayesian Data Analysis, 3th edition. Chapman and Hall/CRC. doi:10.1201/b16018}
#' @export
#'
#'

posterior_predict.flexreg <- function(model, newdata=NULL)
{
  model.name <- model$model@model_name
  N <- length(model$response)
  nsim <- dim(model$model)[1]*dim(model$model)[2]

  link.mu <- model$link.mu
  link.phi <- model$link.phi

  covariate.names.mu <- colnames(model$design.X)
  covariate.names.phi <- colnames(model$design.Z)

  ##############################
  # Estraggo e salvo le catene:
  ##############################


  # Se non ho dei nuovi dati da prevedere, predico i dati osservati. Le catene sono giÃ  salvate nel modello
  # quindi estraggo la catena ed applico la funzione di stima scelta (di default la media a posteriori)
  if(is.null(newdata)){
    #Ncol <- N
    mu.chain <- rstan::extract(model$model, pars="mu", permuted=T)[[1]]
    phi.chain <- rstan::extract(model$model, pars="phi", permuted=T)[[1]]
  } else { # --> Se ho dati da prevedere
    if(!all(colnames(newdata) %in% unique(c(covariate.names.mu, covariate.names.phi))))
      stop("`newdata` must contain the same predictors as the one ones in formula")
    N <- nrow(newdata)
    if((c("(Intercept)") %in% colnames(newdata)) ==F & # If newdata does not have the intercept..
       ((c("(Intercept)") %in% colnames(model$design.X)) ==T) |
       (c("(Intercept)") %in% colnames(model$design.Z) ==T)) # .. and the intercept is required..
      newdata$`(Intercept)` <- rep(1, nrow(newdata)) # .. add the intercept

    newdata.mu  <- newdata[,which(colnames(newdata) %in% colnames(model$design.X) )]
    newdata.phi <- newdata[,which(colnames(newdata) %in% colnames(model$design.Z) )]
    mu.chain <-   mu.chain.nd(model$model, newdata.mu , link.mu)
    phi.chain <-   phi.chain.nd(model$model, newdata.phi , link.phi)
  }

  #additional parameters
  if(model.name %in% c("FB", "FB_phi")){
    p.chain <- rstan::extract(model$model, pars="p", permuted=T)[[1]]
    w.chain <- rstan::extract(model$model, pars="w", permuted=T)[[1]]
  }

  if(model.name %in% c("VIB", "VIB_phi")){
    p.chain <- rstan::extract(model$model, pars="p", permuted=T)[[1]]
    k.chain <- rstan::extract(model$model, pars="k", permuted=T)[[1]]
  }

  post.pred <- matrix(NA, ncol=N, nrow=nsim)
  for(l in 1:N){
    mu.post <- mu.chain[,l]
    phi.post <- ifelse(model.name %in% c("Beta_phi", "FB_phi", "VIB_phi"), phi.chain[,l], phi.chain)
    if(model.name %in% c("Beta", "Beta_phi")) {
      param <- cbind(mu.post, phi.post)
      type <- "Beta"
    } else if(model.name %in% c("FB", "FB_phi")) {
      param <- cbind(mu.post, phi.post,p.chain, w.chain)
      type <- "FB"
    } else if(model.name %in% c("VIB", "VIB_phi")) {
      param <- cbind(mu.post, phi.post,p.chain, k.chain)
      type <- "VIB"
    } else stop("Model not defined")

    post.pred[,l] <- if(type == "Beta") rBeta_mu(n=nsim, mu=param[,1], phi=param[,2]) else
      if(type == "FB") rFB(n=nsim, mu=param[,1], phi=param[,2], p=param[,3], w=param[,4]) else
        if(type=="VIB") rVIB(n=nsim, mu=param[,1], phi=param[,2], p=param[,3], k=param[,4])

  }
  class(post.pred) <- "flexreg_postpred"
  return(post.pred)
}

#' @title posterior_predict
#' @export
#' @keywords internal
#'
# sistemare qui...gli argomenti della funzione devono matchare con quelli della precedente
posterior_predict <- function(model, newdata=NULL)
{
  UseMethod("posterior_predict")
}


#' Posterior Predictives Plot
#'
#' A method for plotting the simulated posterior predictive distribution. It takes an object of class \code{`flexreg_postpred`} (i.e., the output of \code{posterior_predictive.flexreg}) and plots the posterior predictive interval for each statistical unit.
#' Additionally, the mean of the posterior predictives and the value of the response variable can be added.
#'
#'
#' @param x an object of class \code{`flexreg_postpred`} containing the simulated posterior predictives.
#' @param prob the interval probability for the posterior predictives (default is 0.9).
#' @param p_mean a logical value indicating whether the posterior predictives' mean should be plotted.
#' @param response a numerical vector containing the response variable to be added to the plot. if \code{NULL}, response values are omitted.
#' @param ... additional arguments. Currently not used.
#'
#' @examples
#' data("Reading")
#' dataset <- Reading
#' FB <- flexreg(accuracy ~ iq, dataset, n.iter=1000)
#' pp <- posterior_predict(FB)
#' plot(pp)
#'
#' @import stats graphics
#'
#' @method plot flexreg_postpred
#'
#' @export
#'
#'
#'
plot.flexreg_postpred <- function(x, prob=0.9, p_mean=F, response=NULL, ...){
  id <- pp_inf <- pp_sup <- y <- NULL
  post.pred <- x
  N <- ncol(post.pred)
  # limits computation (quantile):
  pred.int <- t(apply(post.pred, 2, function(x) quantile(x, p=c((1-prob)/2, prob+(1-prob)/2))))

  #Definition of the data.frame involved by ggplot:
  dd <- data.frame(pp_inf=pred.int[,1],
                   pp_sup=pred.int[,2],
                   id=c(1:N))

  ll <-
    ggplot(dd) + theme_minimal() +
    scale_y_continuous(name="Posterior Predictive", breaks=seq(0,1,by=.25), limits=c(0,1.1))+
    geom_point(aes(x=id, y=pp_inf), shape="-",size=5)+
    geom_point(aes(x=id, y=pp_sup), shape="-",size=5)+
    labs(x="Unit ID")
  # Add segments:
  for(e in dd$id) {
    ll <- ll + geom_segment(x=e, y=dd$pp_inf[e], xend=e, yend=dd$pp_sup[e], linetype="dashed",size=.1)
  }

  # If p_mean = T, add the posterior predictive mean
  if(p_mean){
    pred.mean <- as.numeric(t(apply(post.pred, 2, function(x) mean(x))))
    ll <- ll + geom_point(data=data.frame(x=c(1:N),y=pred.mean),aes(x=x,y=y), color="#0072B2", size=1.5)
  }
  # Add the response:
  if(!is.null(response)) {
    ll <- ll + geom_point(data=data.frame(x=c(1:N),y=response),aes(x=x,y=y), color="#D55E00", size=1.5)
  }

  return(ll)
}
# plot.flexreg_postpred <- function(object, prob=0.9, p_mean=F, response=NULL, ...){
#   post.pred <- object
#   N <- ncol(post.pred)
#   pred.int <- t(apply(post.pred, 2, function(x) quantile(x, p=c((1-prob)/2, prob+(1-prob)/2))))
#
#
#   plot(pred.int[,1] ~ c(1:N), pch=-95, ylim=c(0,1.1), ylab="Posterior Predictive", xlab="Unit ID", main="Posterior Predictive")
#   points(c(1:N),pred.int[,2], pch=-95)
#   if(p_mean){
#     pred.mean <- t(apply(post.pred, 2, function(x) mean(x)))
#     points(c(1:N),pred.mean, pch=20, col=4)
#   }
#   if(!is.null(response)) points(c(1:N), response, pch=20, col=2)
#   for(i in 1:N) lines(c(i,i), c(pred.int[i,1],pred.int[i,2]), lty=3)
# }
