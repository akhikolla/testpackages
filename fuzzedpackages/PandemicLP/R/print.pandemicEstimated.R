
#### auxiliar function:
## display pameters selection: especially excluding display estimates seasonal parameters without interpretation.
## applicable to any model: "glogistic", "seasonal", "multi_waves"

## Obsolete ##
# excluding_SeasonalParameters=function(out){
#   #out = object of PandemicEstimated class ( output of pandemic_model)
#
#   aux=paste0("mu[",1:nrow(out$Y$data),"]")
#   if(length(out$seasonal_effect)==3 | is.null(out$seasonal_effect)){ #if gen logistic or multiwaves or gen logistic with three days-effect
#     aux=c(aux,"lp__")
#   } else if( length(out$seasonal_effect)==1 & !is.null(out$seasonal_effect)){   # one day-effect
#     aux=c(aux,"lp__","d_2","d_3")
#   } else if( length(out$seasonal_effect)==2 & !is.null(out$seasonal_effect)){   #two days-effect
#     aux=c(aux,"lp__","d_3")
#   }
#   i=which(!(names(out$fit) %in% aux) )
#   pars=names(out$fit)[i]
#
#   sort(pars)  #alphabetical order
# }



#' Print method for \code{pandemicEstimated} objects
#'
#' The print method for \code{pandemicEstimated} object of class S3 displays a compact summary of the
#' fitted model. See the \strong{Details} section below for descriptions of the different components of the printed
#' output.  For additional summary statistics and diagnostics use \code{\link{summary.pandemicEstimated}}.
#'
#'
#' @method print pandemicEstimated
#' @templateVar pandemicEstimatedArg x
#' @param x An object of S3 class \code{\link{pandemicEstimated-objects}}.
#' @param probs A numeric vector of quantiles of interest. The default is
#' \code{c(0.025,0.5,0.975)}.
#' @param digits Number of digits to use for formatting numbers.
#' @param info TRUE or FALSE: more details for output interpretation. The Default is TRUE.
#' @param ... Currently unused.
#' @return Returns \code{x}, invisibly.
#' @details
#'
#' \subsection{Point estimates}{
#' Regardless of the estimation algorithm, point estimates are mean and (or) quantiles computed from simulations.
#' For models fit using MCMC ("sampling", this is default algorithim of \code{pandemic_model} function),  the posterior sample
#' is used. For others estimation algorithm see \code{\link[rstan]{sampling}}  (\pkg{rstan} package).
#'
#' }
#'
#' \subsection{Convergence and efficiency diagnostics for Markov Chains}{
#'
#' Included in the print are: split effective sample sizes (n_eff) and split Rhats.
#'
#'
#' The R-hat convergence diagnostic compares the
#' between- and within-chain estimates for model parameters and other univariate
#' quantities of interest. If chains have not mixed well (ie, the between- and
#' within-chain estimates don't agree), R-hat is larger than 1.
#' We recommend running at least four chains by default and only using the
#' sample if R-hat is less than 1.05.
#'
#' }
#'
#' \subsection{Priors}{
#'
#' A list with information about the prior distributions used and model restrictions (if there are any).
#' For more information go to \code{\link{models}}.
#'
#'
#' }
#'
#' @seealso \code{\link{summary.pandemicEstimated}}.
#'
#' @importMethodsFrom rstan summary
#'
#' @export
print.pandemicEstimated=function(x,digits=3,probs=c(0.025,0.5,0.975),info=TRUE,...){


  if(class(digits)!= "numeric" | class(info)!= "logical"){
    stop("error in 'digits' or 'info'. View ?print.pandemicEstimated")
  }

  if(x$n_waves==1 && is.null(x$seasonal_effect)){   # gen logistic
    name="static generalized logistic"
  } else if(is.null(x$seasonal_effect)==FALSE){         #gen logistic with seasonal effect
    name="static seasonal generalized logistic"
  } else {
    name=paste0("multi_waves(",x$n_waves,")")      #multiwaves
  }

  cat("pandemic_model")
  cat("\n Distribution:       ", "poisson")
  cat("\n Mean function form: ", name)
  cat("\n Type of Case:       ", x$cases.type)
  cat("\n Location:           ", x$Y$name)
  if(is.null(x$seasonal_effect)==FALSE){
  cat("\n Seasonal effect:    ", paste0(x$seasonal_effect,"(","d_",1:length(x$seasonal_effect),")" ) )
  }
  cat("\n 0bservations:       ", nrow(x$Y$data),"\n")


  # display parameters selection:
  pars=names(x$fit)[which( !grepl("mu",names(x$fit)) & !grepl("lp__",names(x$fit)) )]

  cat("\n------\n")
  cat("Parameters:\n")
  tab=rstan::summary(x$fit,pars=pars,probs=probs)$summary[,-c(2:3)] #excludes diagnostic statistic
  tab=round(tab,digits)
  print(tab)


  if(x$n_waves==1 && is.null(x$seasonal_effect)){   #gen logistic
    cat("\n------\n")
    cat("Priors:\n")
    cat("\n a ","~ ","Gamma(0.1, 0.1)")
    cat("\n b ","~ ","LogNormal(0, 20)")
    cat("\n c ","~ ","Gamma(2, 9)")
    cat("\n f ","~ ","Gamma(0.01, 0.01)\n")
  }

  if(is.null(x$seasonal_effect)==FALSE){             #gen logistic with seasonal effect
    cat("\n------\n")
    cat("Priors:\n")
    cat("\n a   ","~ ","Gamma(0.1, 0.1)")
    cat("\n b   ","~ ","LogNormal(0, 20)")
    cat("\n c   ","~ ","Gamma(2, 9)")
    cat("\n d_i ","~ ","Gamma(2,1)")
    cat("\n f   ","~ ","Gamma(0.01, 0.01)\n")
  }

  if(x$n_waves==2){
    cat("\n------\n")
    cat("Priors:\n")
    cat("\n a_i     ","~ ","Gamma(0.1, 0.1)")
    cat("\n alpha_i ","~ ","Gamma(0.01, 0.01)")
    cat("\n b_i     ","~ ","LogNormal(0, 20)")
    cat("\n c_i     ","~ ","Gamma(2, 9)")
    cat("\n delta_i ","~ ","Normal(0, 100)\n")
  }

  if(x$n_waves==1){              #gen logistic ou gen logistic with seasonal effect
    cat("\nRestrictions:")
    cat("\n 1: ", "a/b^f","<",x$config.inputs$use_inputs$p,"*population")
    cat("\n 2: ", "f > 1\n")
  }

  if(x$n_waves==2){              #multiwaves:  2waves
    cat("\nRestrictions:")
    cat("\n 1: ", "a_1/b_1","<",x$config.inputs$use_inputs$p,"*population")
    cat("\n 2: ", "a_2/b_2","<",x$config.inputs$use_inputs$p,"*population\n")
  }

  if(info){
  cat("\n------\n")
  cat("*For help interpreting the printed output see ?print.pandemicEstimated\n")
  cat("*For more information see ?'summary.pandemicEstimated\n")
  cat("*For details on the model, priors and restrictions, see the Details section in ?pandemic_model\n")
  #cat("**a/b^f represents the assymptote of the cumulative cases curve")
  }

  invisible(x)

}


