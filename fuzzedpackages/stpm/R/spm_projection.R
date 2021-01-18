#'A data projection with previously estimated or user-defined parameters. 
#'Projections are constructed for a cohort with fixed or 
#'normally distributed initial covariates. 
#'@references Yashin, A. et al (2007), 
#'Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.
#'@references Akushevich I., Kulminski A. and Manton K. (2005). 
#'Life tables with covariates: Dynamic model 
#'for Nonlinear Analysis of Longitudinal Data. 
#'Mathematical Popu-lation Studies, 12(2), pp.: 51-80.
#'<DOI: 10.1080/08898480590932296>.
#'@references Yashin, A. et al (2007), Health decline, 
#'aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#'@param x A list of parameters from output of the \code{spm(...)} function.
#'@param N A number of individuals to simulate, N=100 by default.
#'@param ystart A vector of starting values of covariates (variables), ystart=80 by default.
#'@param model A model type. Choices are: "discrete", "continuous" or "time-dependent".
#'@param tstart Start time (age), default=30. Can be an interval: c(a, b) - in this case,
#'the starting time is sumulated via \code{runif(1, a, b)}.
#'@param tend End time (age), default=105.
#'@param dt A time interval between observations, dt=1 by default.
#'@param sd0 A standard deviation value for simulation of the next value of variable.
#'sd0=1 by default.
#'@param nobs A number of observations (lines) for i-th individual.
#'@param gomp A flag (FALSE by default). 
#'When it is set, then time-dependent exponential form of mu0 and Q are used:
#' mu0 = mu0*exp(theta*t), Q = Q*exp(theta*t). 
#' Only for continous-time SPM.
#' @param format Data format: "short" (default), "long".
#'@return An object of 'spm.projection' class with two elements. 
#'(1) A simulated data set.
#'(2) A summary statistics which includes (i) age-specific means of state variables and
#'(ii) Survival probabilities.
#'@export
#'@examples \dontrun{ 
#'library(stpm)
#'# Setting up the model
#'model.par <- list()
#'model.par$a <- matrix(c(-0.05, 1e-3, 2e-3, -0.05), nrow=2, ncol=2, byrow=TRUE)
#'model.par$f1 <- matrix(c(90, 35), nrow=1, ncol=2)
#'model.par$Q <- matrix(c(1e-8, 1e-9, 1e-9, 1e-8), nrow=2, ncol=2, byrow=TRUE)
#'model.par$f <- matrix(c(80, 27), nrow=1, ncol=2)
#'model.par$b <- matrix(c(6, 2), nrow=2, ncol=2)
#'model.par$mu0 <- 1e-6
#'model.par$theta <- 0.09
#'# Projection
#'# Discrete-time model
#'data.proj.discrete <- spm_projection(model.par, N=5000, ystart=c(80, 27))
#'plot(data.proj.discrete$stat$srv.prob)
#'# Continuous-time model
#'data.proj.continuous <- spm_projection(model.par, N=5000, 
#'ystart=c(80, 27), model="continuous")
#'plot(data.proj.continuous$stat$srv.prob)
#'# Time-dependent model
#'model.par <- list(at = "-0.05", f1t = "80", Qt = "2e-8", 
#'ft= "80", bt = "5", mu0t = "1e-5*exp(0.11*t)")
#'data.proj.time_dependent <- spm_projection(model.par, N=500, 
#'ystart=80, model="time-dependent")
#'plot(data.proj.time_dependent$stat$srv.prob, xlim = c(30,105))
#'}
spm_projection <- function(x, 
                           N=100, 
                           ystart=80, 
                           model="discrete", 
                           tstart=30, tend=105, 
                           dt=1, 
                           sd0=1, 
                           nobs=NULL, 
                           gomp=TRUE, 
                           format="short") {
  
  avail.models <- c("discrete", "continuous", "time-dependent")
  if(!(model %in% avail.models)) {
    stop(paste("Provided model", model, "not found in the list of available models."))
  }
  
  if(length(tstart) > 2) {
    stop(paste("Incorrect tstart:", tstart))
  }
  
  res <- list()
  res.sim <- list()
  # Statistics
  stat <- list()
  
  if(model == "time-dependent") {
    # Data simulation for time-dependent model
    
    formulas.work <- list(at = "-0.05", f1t = "80", Qt = "2e-8", 
                          ft= "80", bt = "5", mu0t = "1e-5*exp(0.11*t)")
    
    if (!is.null(x)) {
      formulas.work <- x
    }
    
    #Simulate (project) data:
    res.sim <- simdata_time_dep(N=N,f=formulas.work,
                                    step=dt, 
                                    tstart=tstart, 
                                    tend=tend, 
                                    ystart=ystart, 
                                    sd0=sd0, nobs=nobs,
                                format=format)
    #Computing summary statistics
    ## Age-specific means:
    bins<-10
    cutpoints<-quantile(res.sim[,3],(0:bins)/bins)
    binned <-cut(res.sim[,3],cutpoints,include.lowest=TRUE)
    for(i in seq(5,length(colnames(res.sim)),by=2)) {
      mean.cov <- tapply(res.sim[,5], binned, mean)
      stat[["mean.by.age"]][[colnames(res.sim)[i]]] <- mean.cov
    }
    
    ##Survival probabilities:
    vt <- data.frame(matrix(nrow=N, ncol=5))
    df <- data.frame(res.sim)
    ddg <- split(df, f=df$id)
    
    invisible(sapply(1:N, 
                     function(i) {  dg <- ddg[[i]]
                     vt[i,] <<- c(dg[1,1], # id
                                  dg[dim(dg)[1], 4] - dg[1, 3], # time
                                  dg[dim(dg)[1], 4], # age
                                  dg[dim(dg)[1], 2], # case
                                  dg[1, 3])} # start
    ))
    
    colnames(vt) = c("id", "time", "age","case", "start")
    
    ##Survival probabilities:
    srv.prob <- survfit( Surv(start, age, case) ~ 1, data = vt, conf.type = "log-log")
    stat[["srv.prob"]] <- srv.prob
    
  } else if(model == "discrete") {
    
    if(length(x$a) != length(ystart)^2) {
      stop("Number of dimensions does not match with the number of values provided in ystart.")
    }
    
    res.sim <- simdata_discr(N=N, 
                               a=x$a, 
                               f1=x$f1, 
                               Q=x$Q, 
                               f=x$f, 
                               b=x$b, 
                               mu0=x$mu0, 
                               theta=x$theta, 
                               ystart=ystart, 
                               tstart=tstart, tend=tend, 
                               dt=dt,
                             format=format)
    
   
    # Age-specific means:
    bins<-10
    cutpoints<-quantile(res.sim[,3],(0:bins)/bins)
    binned <-cut(res.sim[,3],cutpoints,include.lowest=TRUE)
    for(i in seq(5,length(colnames(res.sim)),by=2)) {
      mean.cov <- tapply(res.sim[,5], binned, mean)
      stat[["mean.by.age"]][[colnames(res.sim)[i]]] <- mean.cov
    }
    
    vt <- data.frame(matrix(nrow=N, ncol=5))
    df <- data.frame(res.sim)
    ddg <- split(df, f=df$id)
    
    invisible(sapply(1:N, 
                     function(i) {  dg <- ddg[[i]]
                     vt[i,] <<- c(dg[1,1], # id
                                  dg[dim(dg)[1], 4] - dg[1, 3], # time
                                  dg[dim(dg)[1], 4], # age
                                  dg[dim(dg)[1], 2], # case
                                  dg[1, 3])} # start
    ))
    
    colnames(vt) = c("id", "time", "age","case", "start")
    
    #Survival probabilities:
    srv.prob <- survfit( Surv(start, age, case) ~ 1, data = vt, conf.type = "log-log")
    stat[["srv.prob"]] <- srv.prob
    
    
  } else if(model == "continuous") {
    
    if(length(x$a) != length(ystart)^2) {
      stop("Number of dimensions does not match with the number of values provided in ystart.")
    }
    
    # Data simulation for discrete and continuous models
    res.sim <- simdata_cont(N=N, 
                             a=x$a, 
                             f1=x$f1, 
                             Q=x$Q, 
                             f=x$f, 
                             b=x$b, 
                             mu0=x$mu0, 
                             theta=x$theta,
                             dt=dt, 
                             ystart=ystart,
                             tstart=tstart, tend=tend, 
                             sd0=sd0, nobs=nobs, gomp=gomp,
                            format=format)
    
    
    # Age-specific means:
    bins<-10
    cutpoints<-quantile(res.sim[,3],(0:bins)/bins)
    binned <-cut(res.sim[,3],cutpoints,include.lowest=TRUE)
    for(i in seq(5,length(colnames(res.sim)),by=2)) {
      mean.cov <- tapply(res.sim[,5], binned, mean)
      stat[["mean.by.age"]][[colnames(res.sim)[i]]] <- mean.cov
    }
    
    ####### Survival probabilities ########
    vt <- data.frame(matrix(nrow=N, ncol=5))
    df <- data.frame(res.sim)
    ddg <- split(df, f=df$id)
    
    invisible(sapply(1:N, 
                     function(i) {  dg <- ddg[[i]]
                                    vt[i,] <<- c(dg[1,1], # id
                                                dg[dim(dg)[1], 4] - dg[1, 3], # time
                                                dg[dim(dg)[1], 4], # age
                                                dg[dim(dg)[1], 2], # case
                                                dg[1, 3])} # start
    ))
    
    colnames(vt) = c("id", "time", "age","case", "start")
    
    srv.prob <- survfit( Surv(start, age, case) ~ 1, data = vt, conf.type = "log-log")
    stat[["srv.prob"]] <- srv.prob
    
  }
  
  res <- list(data=res.sim, stat=stat)
  
  class(res) <- "spm.projection"
  invisible(res)
}