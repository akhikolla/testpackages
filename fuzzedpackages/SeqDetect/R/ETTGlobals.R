merge_named_lists <- function(list1, list2, fun=NULL) {
  if(!is.list(list1) || !is.list(list2))
    stop("function merge_named_list requires two lists")
  if(is.null(fun) || !is.function(fun))
    fun <- sum
  res_list <- list1
  for(n in names(list2)) {
    if(!n %in% names(res_list)) res_list[[n]] <- list2[[n]]
    else res_list[[n]] <- fun(res_list[[n]],list2[[n]])
  }
  return(res_list)
}

c_to_string <- function(var) {
  s <- "["
  if(length(var)>0)
    for(i in 1:length(var)) {
      if(i<length(var)) s <- paste0(s,var[[i]],",")
      else s <- paste0(s,var[[i]])
    }
  return(paste0(s,"]"))
}