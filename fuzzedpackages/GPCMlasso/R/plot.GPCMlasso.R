#' Plot function for GPCMlasso
#' 
#' Plot function for a \code{GPCMlasso} object. Plots show coefficient paths
#' of DIF (or DSF) parameters along (a transformation of) the tuning parameter lambda.
#' One plot per item is created, every single parameter corresponding to this item
#' is depicted by a single path. 
#' The optimal model is highlighted with a red dashed line. 
#' 
#' @usage \method{plot}{GPCMlasso}(x, select = c("BIC", "AIC", "cAIC", "cv"),
#' log.lambda = TRUE, items_per_page = 1, items = "all", 
#' columns = NULL, ask_new = TRUE, lambda.lines = TRUE,
#' equal_range = TRUE, \dots)
#' @param x \code{GPCMlasso} object
#' @param select Specifies which criterion to use for the optimal model, we recommend the 
#' default value "BIC". If cross-validation was performed, automatically the optimal
#' model according to cross-validation is used. The chosen optimal model is 
#' highlighted with a red dashed line. 
#' @param log.lambda A logical value indicating whether lambda or a log-transformation of 
#' lambda should be used as x-axis in the plots.
#' @param items_per_page By default, each plot/item is put on a separate page. For example,
#' \code{items_per_page=4} would put four plots/items on one page.
#' @param items By default, all items are plotted. If \code{items=c(1,3)}, 
#' only the first and the third item are plotted.
#' @param columns Specifies the number of columns to use when several
#' plots are on one page. Only relevant if \code{items_per_page}>1. 
#' @param ask_new If TRUE, the user is asked to confirm before the next item is plotted.
#' @param lambda.lines A logical value indicating whether a thin gray line plotted
#' for each value from the vector of tuning parameters from \code{object}
#' @param equal_range A logical value indicating whether for each plot equal limits 
#' on the y-axis shall be used.
#' @param ... Further plot arguments.
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @references Schauberger, Gunther and Mair, Patrick (2019): A Regularization Approach for the Detection of Differential 
#' Item Functioning in Generalized Partial Credit Models, \emph{Behavior Research Methods}, \url{https://link.springer.com/article/10.3758/s13428-019-01224-2}
#' @seealso \code{\link{GPCMlasso}}
#' @examples
#' data(tenseness_small)
#' 
#' ## formula for simple model without covariates
#' form.0 <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~0"))
#' 
#' ######
#' ## fit simple RSM where loglikelihood and score function are evaluated parallel on 2 cores
#' rsm.0 <- GPCMlasso(form.0, tenseness_small, model = "RSM", 
#' control= ctrl_GPCMlasso(cores=2))
#' rsm.0
#' 
#' \dontrun{
#' ## formula for model with covariates (and DIF detection)
#' form <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~."))
#' 
#' ######
#' ## fit GPCM model with 10 different tuning parameters
#' gpcm <- GPCMlasso(form, tenseness_small, model = "GPCM", 
#'                   control = ctrl_GPCMlasso(l.lambda = 10))
#' gpcm
#' plot(gpcm)
#' pred.gpcm <- predict(gpcm)
#' trait.gpcm <- trait.posterior(gpcm)
#' 
#' ######
#' ## fit RSM, detect differential step functioning (DSF)
#' rsm.DSF <- GPCMlasso(form, tenseness_small, model = "RSM", DSF = TRUE, 
#'                      control = ctrl_GPCMlasso(l.lambda = 10))
#' rsm.DSF
#' plot(rsm.DSF)
#' 
#' ## create binary data set
#' tenseness_small_binary <- tenseness_small
#' tenseness_small_binary[,1:5][tenseness_small[,1:5]>1] <- 2
#' 
#' ######
#' ## fit and cross-validate Rasch model
#' set.seed(1860)
#' rm.cv <- GPCMlasso(form, tenseness_small_binary, model = "RM", cv = TRUE, 
#'                    control = ctrl_GPCMlasso(l.lambda = 10))
#' rm.cv
#' plot(rm.cv)
#' }
plot.GPCMlasso <- function(x, select = c("BIC", "AIC", "cAIC", "cv"),
                           log.lambda = TRUE, items_per_page = 1, items = "all", 
                           columns = NULL, ask_new = TRUE, lambda.lines = TRUE,
                           equal_range = TRUE, ...){
  op <- par(no.readonly = TRUE)

  select <- match.arg(select, c("BIC", "AIC", "cAIC", "cv"))
  
  if(select=="BIC"){
    criterion <- x$BIC
  }
  if(select=="AIC"){
    criterion <- x$AIC
  }
  if(select=="cAIC"){
    criterion <- x$cAIC
  }
  if(select=="cv" | !is.null(x$cv_error)){
    criterion <- x$cv_error
    if(is.null(x$cv_error)){
      warning("No cross-validation was performed for the model. Instead of the cross-validation error the optimal model is selected by BIC")
      criterion <- x$BIC
    }
  }
  
    with(x$design_list,{
    if(m==0){
      stop("No covariates, nothing to plot!")
    }
    if(length(x$control$lambda)==1){
      stop("Only one tuning parameter, nothing to plot!")
    }
    
      n.dif.par <- I*m
      if(x$DSF){
        n.dif.par <- sum(q)*m
      }

  if(identical(items, "all")){items <- 1:I}
    
  n.plots <- length(items)
  pages <- ceiling(n.plots/items_per_page)
    
  coefs <- x$coefficients
  
  gamma.start <- sum(q)+1
  if(RSM){
    gamma.start <- q[1]+I
  }
  
  if(x$main.effects){
    gamma.start <- gamma.start + m
  }
  
  gamma <- coefs[,gamma.start:(gamma.start+n.dif.par-1)]

  
  if (is.null(columns)) {
    cols <- floor(sqrt(items_per_page))
  } else {
    cols <- columns
  }
  rows <- ceiling(items_per_page/cols)
  
  g.range = range(gamma, na.rm = TRUE)

  start.gamma <-1
  plots_on_page <- 0
  pages_done <- 0
  par(mfrow=c(rows,cols),xpd=TRUE)
  
  for(u in 1:I){


      if(x$DSF){
        par.item <- q[u]*m
      }else{
        par.item <- m
      }
    if(u %in% items){
    plot.item(gamma[,start.gamma:(start.gamma+par.item-1), drop = FALSE], 
           x$item.names[u], par.item, x$control$lambda, g.range, 
           equal_range, criterion, x.names, log.lambda,
           lambda.lines)
    plots_on_page <- plots_on_page+1
    if(plots_on_page==items_per_page & pages_done<(pages-1)){
      plots_on_page <- 0
      pages_done <- pages_done+1
      if(interactive() & ask_new)
      {readline("Press enter for next plot!")}
      par(mfrow=c(rows,cols),xpd=TRUE)
    }
    }
    start.gamma <- start.gamma+par.item    }
  
  })
  par(op)
  invisible(x)
}

plot.item <- function(item, name, par.item, lambda, g.range, equal_range, 
                      criterion,  x.names, log.lambda, lambda.lines){

    if(equal_range){
    i.lim <- g.range
  }else{
    i.lim <- range(item, na.rm = TRUE)
  }

  if(log.lambda){
    xlab <- expression(log(lambda+1))
    lambda <- log(lambda+1)
  }else{
    xlab <- expression(lambda)
  }
  
  x.lim <- rev(range(c(lambda,min(lambda)-abs(diff(range(lambda)))*0.1)))
  
  plot(lambda, item[,1],ylim=i.lim,type="l", main=name,xlim=x.lim,
       xlab=xlab,ylab="parameters",frame.plot=FALSE,
       lwd=par()$lwd)
  if(par.item>1){
    for(uu in 2:par.item){
      lines(lambda, item[,uu],lwd=par()$lwd)
    }
  }
  if(lambda.lines){
    segments(x0=lambda,y1=i.lim[1],y0=i.lim[2],lwd=.5, lty=2, col="lightgray")
    # abline(v=lambda, lwd=.5, lty=2, col="lightgray")
  }
  
  
  if(!is.null(criterion)){
  segments( lambda[which.min(criterion)], max(i.lim),
            lambda[which.min(criterion)], min(i.lim) ,
            col=2,lty=2,lwd=par()$lwd)
  }

  x.lab1 <- min(lambda)-abs(diff(range(lambda)))*0.02
  x.lab2 <- min(lambda)-abs(diff(range(lambda)))*0.005
  y.lab1 <- item[length(lambda),]
  y.lab2 <- spread.labs(y.lab1, 1.2*strheight("A"))
  
  ncat <- par.item/length(x.names)
  if(ncat>1){
    x.names <- paste(rep(x.names,ncat),rep(1:ncat,each=length(x.names)),sep=".")
  }
  
  text( x.lab1, y.lab2, x.names,pos=4)
  segments( x.lab2, y.lab1,
            x.lab1, y.lab2 ,col="gray")
  
  
  ## used for illustrating plot in paper
  # segments( 3, i.lim[1],
  #           3, i.lim[2] ,col="gray",lty=2, lwd = 2)
  
}
