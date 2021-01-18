#' @export
print.mumm <- function(x, ...) {

  cat("Multiplicative mixel model fit by ML \n Formula:")

  for(i in 1:length(deparse(x$call$formula,width.cutoff = 50))) {
    if(i>1){cat("      ")}
    cat(deparse(x$call$formula,width.cutoff = 50)[i],"\n")
  }


  cat("Data:", x$call$data, "\n")

  cat("Log-likelihood at convergence:", -x$objective, "\n")

  cat("Random effects: \n")
  #table with variance components
  table_sd = data.frame(Groups = names(x$sigmas), Std.Dev. = x$sigmas)
  print(table_sd, right = FALSE, row.names = FALSE)

  cat("Correlation:", x$est_cor, "\n")

  cat("Number of obs:", x$nobs, "\n")

  cat("Fixed Effects: \n")
  print.default(x$par_fix)

}

#' Extract Random Effects
#'
#' A function to extract the estimated random effects from a model object of class \code{mumm}.
#'
#'
#' @param object an object of class "mumm"
#'
#' @param ... Currently not used
#'
#' @return A named list with the estimated random effects, where each element in the list is
#' a numeric vector consisting of the estimated random effect coefficients for a random factor in the model.
#'
#' @examples
#' set.seed(100)
#' sigma_e <- 1.5
#' sigma_a <- 0.8
#' sigma_b <- 0.5
#' sigma_d <- 0.7
#' nu <- c(8.2, 6.2, 2.3, 10.4, 7.5, 1.9)
#'
#' nA <- 15
#' nP <- 6
#' nR <- 5
#'
#' a <- rnorm(nA, mean = 0, sd = sigma_a)
#' b <- rnorm(nA, mean = 0, sd = sigma_b)
#' d <- rnorm(nA*nP, mean = 0, sd = sigma_d)
#' e <- rnorm(nA*nP*nR, mean = 0, sd = sigma_e)
#'
#' Assessor <- factor(rep(seq(1,nA),each = (nP*nR)))
#' Product <- factor(rep(rep(seq(1,nP),each = nR), nA))
#' AssessorProduct <- (Assessor:Product)
#'
#' y <- nu[Product] + a[Assessor] + b[Assessor]*(nu[Product]-mean(nu)) + d[AssessorProduct] + e
#'
#' sim_data <- data.frame(y, Assessor, Product)
#'
#' fit <- mumm(y ~ 1 + Product + (1|Assessor) + (1|Assessor:Product) +
#'              mp(Assessor,Product) ,data = sim_data)
#'
#' ranef(fit)
#'
#' @importFrom lme4 ranef
#' @export
ranef.mumm <- function(object, ...) {

  rand_ef = names(object$sigmas)
  names_random_par = names(object$par_rand)
  nlevels_par = object$nlevels_par_rand
  par_rand_list = list()

  index = 1
  for(i in 1:length(object$nlevels_par_rand)) {
    par_rand_list[[as.character(rand_ef[i])]] = object$par_rand[index:(index+nlevels_par[[i]]-1)]
    index = index + nlevels_par[[i]]
  }

  print.default(par_rand_list)

}

#' Confidence Intervals for Model Parameters
#'
#' Computes confidence intervals for the fixed effect parameters and the variance components
#' for an object of class \code{mumm}.
#'
#' @param object an object of class mumm.
#'
#' @param parm a vector of parameter names or a matrix, where the rows specify linear combinations of the model parameters. If missing,
#' confidence intervals will be computed for all of the fixed effect parameters and all of the variance components.
#'
#' @param level the confidence level.
#'
#' @param ... Currently not used.
#'
#' @details The confidence intervals are computed by the profile likelihood method.
#'
#' @return A matrix with the first column showing the lower confidence limit and the second column showing the
#' upper limit for each parameter or linear combination of parameters.
#'
#' @examples
#' set.seed(100)
#' sigma_e <- 1.5
#' sigma_a <- 0.8
#' sigma_b <- 0.5
#' sigma_d <- 0.7
#' nu <- c(8.2, 6.2, 2.3, 10.4, 7.5, 1.9)
#'
#' nA <- 15
#' nP <- 6
#' nR <- 5
#'
#' a <- rnorm(nA, mean = 0, sd = sigma_a)
#' b <- rnorm(nA, mean = 0, sd = sigma_b)
#' d <- rnorm(nA*nP, mean = 0, sd = sigma_d)
#' e <- rnorm(nA*nP*nR, mean = 0, sd = sigma_e)
#'
#' Assessor <- factor(rep(seq(1,nA),each = (nP*nR)))
#' Product <- factor(rep(rep(seq(1,nP),each = nR), nA))
#' AssessorProduct <- (Assessor:Product)
#'
#' y <- nu[Product] + a[Assessor] + b[Assessor]*(nu[Product]-mean(nu)) + d[AssessorProduct] + e
#'
#' sim_data <- data.frame(y, Assessor, Product)
#'
#' fit <- mumm(y ~ 1 + Product + (1|Assessor) + (1|Assessor:Product) +
#'              mp(Assessor,Product) ,data = sim_data)
#' \donttest{
#' confint(fit, parm = c('Product3', 'mp Assessor:Product'), level = 0.90)
#' }
#'
#' @export
confint.mumm <- function(object, parm = "all", level = 0.95, ...){

  confints = c()
  n_parfix = length(object$par_fix)

  if(is.na(object$est_cor)){
    name_vector = c(names(object$par_fix),names(object$sigmas))
  }else{
    name_vector = c(names(object$par_fix),names(object$sigmas),names(object$est_cor))
  }
  names(object$obj$par) = name_vector

  if(typeof(parm) == "character"){

    if(parm[1] == "all"){
      index = 1:length(name_vector)
    }else{
      index = match(parm,name_vector)
      name_vector <- name_vector[index]
    }


    for(i in 1:length(index)){
      profile = tmbprofile(object$obj, name = index[i] , trace = FALSE)
      c = stats::confint(profile, level = level)

      #If parameter is a variance component
      if(index[i]>n_parfix){
        if(name_vector[i]=="Correlation"){
          confints = rbind(confints,c/sqrt(1. + c))
        }else{
          confints = rbind(confints,exp(c))
        }
      }else{
        confints = rbind(confints,c)
      }


    }
    dimnames(confints)[[1]] <- name_vector

  }else{ #if parm is a matrix with rows containing linear comb.

    if(is.vector(parm)){
    parm = as.matrix(t(parm))
    }

    for(i in 1:nrow(parm)){

      index_vec = 1:length(object$par)
      index = index_vec[as.logical(parm[i,])]

      profile = tmbprofile(object$obj, lincomb = parm[i,] , trace = FALSE)
      c = stats::confint(profile, level = level)

      #If parameter is a variance component
      if(index[1]>n_parfix){
        if(name_vector[index[1]]=="Correlation"){
          confints = rbind(confints,c/sqrt(1. + c))
        }else{
          confints = rbind(confints,exp(c))
        }
      }else{
        confints = rbind(confints,c)
      }

    }
  }


  print.default(confints)
}
