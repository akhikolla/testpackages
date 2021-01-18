#' Names and descriptions of precompiled Stan models
#' @description By default, line 2 and 3 of comments starting
#' with # or // in Stan file are returned
#' @return A data frame with \code{model_name} and the first \code{n_lines}
#' comment lines in model as description
#' @param n_lines Number of comment lines to retrieve
#' @param skip Number of lines to skip from beginning of Stan Model file
#' @param sep separator for multiline strings
#' @importFrom stringr str_match_all
#' @export
#'
stan_model_names = function(n_lines = 2, skip = 1, sep = "\n"){
  df = unlist(lapply(stanmodels, function(model){
    code = model@model_code
    paste(str_match_all(code,"// (.*)|# (.*)")[[1]][(1 + skip):(skip + n_lines),2],
          collapse = sep)
  }))
  data.frame(model_name = names(stanmodels),
             description = df, stringsAsFactors = FALSE)
}

