#####################################
####    print.design_analysis    ####
#####################################

#----    print.design_analysis    ----

#' Print Method for design_analysis class
#'
#' @param x an object with class "design_analysis"
#' @param ... further arguments passed to or from other methods.
#'
#' @return a summary output
#' @export
#' @noRd
#'
print.design_analysis <- function(x, ...){
  da_fit <- x
  prefix <- "\t"
  output_text <- da_fit


  effect_type <- ifelse(da_fit$effect_info$effect_type=="cohen_d",
                        "cohen_d", "rho")
  sigle_value <- da_fit$effect_info$effect_function == "single_value"



  effect_info <- data.frame(effect_type = da_fit$effect_info$effect_type)

  # test_info
  test_info <- da_fit$test_info
  test_info <- as.data.frame(t(round_arg(test_info, 3)))
  critical_effect <- sign_effect(test_info$critical_effect,
                                test_info$alternative)
  test_info[["critical_effect"]] <- NULL

  # res_info
  if(da_fit$design_analysis == "retrospective"){
    res_info <- round_arg(da_fit$retrospective_res, 3)
  } else {
    res_info <- round_arg(da_fit$prospective_res, 3)
  }


  # effect_summary
  effect_summary <- round(summary(da_fit$effect_info$effect_samples),3)
  effect_summary <- c(length(da_fit$effect_info$effect_samples),
                     effect_summary)
  names(effect_summary)[1] <- "n_effect"
  effect_summary <- as.data.frame(t(round(effect_summary,3)))

  cat("\n")
  cat(strwrap("Design Analysis", prefix = prefix), sep = "\n")
  cat("\n")
  cat("Hypothesized effect: ", effect_type)
  if(sigle_value){
    cat(" =", da_fit$effect_info$effect_samples, "\n")
  } else {
    cat(" ~", deparse(da_fit$effect_info$effect_function))

    if(is.finite(da_fit$effect_info$tl) || is.finite(da_fit$effect_info$tu)){
      cat(" [tl = ", da_fit$effect_info$tl,"; tu =", da_fit$effect_info$tu, "]")
    }
    cat("\n")
    cat(capture.output(print(effect_summary,print.gap = 3, right = FALSE,
                             row.names = FALSE)), sep = "\n")

  }


  cat("\n")
  cat("Study characteristics:",
      capture.output(print(test_info,
        print.gap = 3, right = FALSE, row.names = FALSE)), sep = '\n')
  cat("\n")

  cat("Inferential risks:\n")
      if(sigle_value){
        cat(capture.output(print(res_info,
                             print.gap = 3, right = FALSE, row.names = FALSE)),
            sep = '\n')
      } else {
        error_summary <- suppressWarnings(as.data.frame(
          t(do.call(cbind, lapply(res_info, summary)))))
        cat(capture.output(print(error_summary,print.gap = 3, right = FALSE)),
            sep = "\n")

      }

  cat("\n")

  cat("Critical value(s):", effect_type, " = ", critical_effect)
  cat("\n")

  invisible(da_fit)
}

#----    summary.design_analysis    ----

#' Print Method for design_analysis class
#'
#' @param object an object with class "design_analysis"
#' @param ... further arguments passed to or from other methods.
#'
#' @return a summary output
#'
#' @export
#' @noRd
#'
summary.design_analysis <- function(object, ...){
  print(object)
  invisible(object)
}

#----









