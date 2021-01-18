tree2indicators <-
function (fit) 
{
  if (!inherits(fit, "rpart")) 
    stop("fit have to be an rpart object")
  leaves <- path.rpart(fit, nodes = row.names(fit$frame[fit$frame$var == 
                                                          "<leaf>", ]), print.it = FALSE)
  n_leaves = length(leaves)
  list_leaves = lapply(leaves, function(u) return(u[-1]))
  list_leaves_1 = lapply(list_leaves, function(u) {
    r = sapply(u, function(uu) return(sub("=", "==", uu)))
    return(r)
  })
  list_leaves_2 = lapply(list_leaves_1, function(u) {
    r = sapply(u, function(uu) return(sub(">==", ">=", uu)))
    return(r)
  })
  list_leaves_3 = lapply(list_leaves_2, function(u) {
    r = sapply(u, function(uu) return(sub("<==", "<=", uu)))
    return(r)
  })
  list_leaves_4 = lapply(list_leaves_3, function(u) {
    r = sapply(u, function(uu) {
      if (grepl("==", uu)) {
        return(gsub("==", ",c('", gsub("$", "'))", gsub("^", 
                                                        "is.element(", gsub(",", "','", uu))), fixed = T))
      }
      else return(uu)
    })
    return(r)
  })
  indicators = lapply(list_leaves_4, function(u) {
    return(paste(u, collapse = " & "))
  })
  return(indicators)
}
