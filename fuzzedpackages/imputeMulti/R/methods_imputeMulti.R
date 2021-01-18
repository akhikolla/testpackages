
###########################################################
## Methods - mod_imputeMulti
###########################################################
## Print

#' @name show-mod_imputeMulti
#' @aliases show,mod_imputeMulti-method
#' @docType methods
#' @rdname mod_imputeMulti-class
#' @include class_imputeMulti.R
#' @title Methods for mod_imputeMulti class objects
#' @param object an object of class "mod_imputeMulti"
setMethod("show", signature= "mod_imputeMulti",
          def= function(object) {
            cat("\n Call: \n", paste(deparse(object@mle_call), sep= "\n"),
                "\n Method: ", object@method,
                "\n\n Iterations: ", object@mle_iter,
                "\n\n Log-Likelihood: ", object@mle_log_lik)
          })


## Summary
setGeneric("summary")
# @export
summary.mod_imputeMulti <- function(object, ...) {
  methods::show(object)
  cat("\n\n")
  if (object@mle_cp != "none") {
    summary(object@mle_x_y[, c("alpha", "theta_y")])  
  } 
}


#' @title Summarizing mod_imputMulti objects
#' @description summary method for class "mod_imputeMulti"
#' @param object an object of class "mod_imputeMulti"
#' @param ... further arguments passed to or from other methods.
#' @exportMethod summary
setMethod("summary", signature="mod_imputeMulti", def=summary.mod_imputeMulti)

#' @name get_parameters
#' @aliases get_parameters
#' @title Get observation level imputed values 
#' @rdname mod_imputeMulti-class
#' @exportMethod get_parameters
# @export
setGeneric("get_parameters", 
           function(object) standardGeneric("get_parameters"))

#' @rdname mod_imputeMulti-class
setMethod("get_parameters", signature= "mod_imputeMulti",
          function(object) {
            return(object@mle_x_y) 
          })

#' @title get_prior
#' @aliases get_prior
#' @rdname mod_imputeMulti-class
#' @exportMethod get_prior
# @export
setGeneric("get_prior", 
           function(object) standardGeneric("get_prior"))

#' @rdname mod_imputeMulti-class
setMethod("get_prior", signature= "mod_imputeMulti",
          function(object) {
            return(object@mle_cp) 
          })

#' @title get_iterations
#' @aliases get_iterations
#' @rdname mod_imputeMulti-class
#' @exportMethod get_iterations
# @export
setGeneric("get_iterations", 
           function(object) standardGeneric("get_iterations"))

#' @rdname mod_imputeMulti-class
setMethod("get_iterations", signature= "mod_imputeMulti",
          function(object) {
            return(object@mle_iter) 
          })

#' @title get_logLik
#' @aliases get_logLik
#' @rdname mod_imputeMulti-class
#' @exportMethod get_logLik
# @export
setGeneric("get_logLik", 
           function(object) standardGeneric("get_logLik"))

#' @rdname mod_imputeMulti-class
setMethod("get_logLik", signature= "mod_imputeMulti",
          function(object) {
            return(object@mle_log_lik) 
          })

#' @title get_method
#' @aliases get_method
#' @rdname mod_imputeMulti-class
#' @exportMethod get_method
# @export
setGeneric("get_method", 
           function(object) standardGeneric("get_method"))

#' @rdname mod_imputeMulti-class
setMethod("get_method", signature= "mod_imputeMulti",
          function(object) {
            return(object@method) 
          })


###########################################################
## Methods - mod_imputeMulti
###########################################################
## Print

#' @name show-imputeMulti
#' @aliases show,imputeMulti-method
#' @docType methods
#' @rdname imputeMulti-class
#' @include class_imputeMulti.R
#' @title Methods for mod_imputeMulti class objects
#' @param object an object of class "imputeMulti"
setMethod("show", signature= "imputeMulti",
          def= function(object) {
            cat("\n Global Call: \n", paste(deparse(object@Gcall), sep= "\n"),
                "\n Call: \n", paste(deparse(object@mle_call)),
                "\n Method: ", object@method,
                "\n\n Iterations: ", object@mle_iter,
                "\n\n Log-Likelihood: ", object@mle_log_lik,
                "\n Number Missing: ", object@nmiss)
          })

## Summary
setGeneric("summary")
# @export
summary.imputeMulti <- function(object, ...) {
  methods::show(object)
  if (object@mle_cp != "none") {
    summary(object@mle_x_y[, c("alpha", "theta_y")])  
  } 
}

#' @title Summarizing imputMulti objects
#' @description summary method for class "imputeMulti"
#' @param object an object of class "imputeMulti"
#' @param ... further arguments passed to or from other methods.
#' @exportMethod summary
setMethod("summary", signature= "imputeMulti",
          summary.imputeMulti)



#' @name get_imputations
#' @aliases get_imputations
#' @rdname imputeMulti-class
#' @exportMethod get_imputations
#' @export
setGeneric("get_imputations", 
           function(object) standardGeneric("get_imputations"))

#' @rdname imputeMulti-class
setMethod("get_imputations", signature= "imputeMulti",
          function(object) {
            return(object@data$imputed_data) 
          })

#' @name n_miss
#' @aliases n_miss
#' @rdname imputeMulti-class
#' @exportMethod n_miss
#' @export
setGeneric("n_miss", 
           function(object) standardGeneric("n_miss"))

#' @rdname mod_imputeMulti-class
setMethod("n_miss", signature= "imputeMulti",
          function(object) {
            return(object@nmiss) 
          })