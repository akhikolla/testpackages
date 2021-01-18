
"con2df" <- function(con, Traits){
  con <- con[setdiff(names(con), c("lb", "ub", "uniform"))]
  if(length(con)==0){
    return(data.frame(dir=character(0), var=character(0), val=numeric(0), isLin=logical(0)))
  }

  dir <- str_sub(names(con),1,2)  
  var <- str_sub(names(con),4,-1)

  if(any(duplicated(var))){
    stop("Some variables appear in more than one constraint.\n")
  }  
  
  for(i in seq_along(con)){
    if((!is.numeric(con[[i]])) || length(con[[i]])>1 || is.na(con[[i]])){
      stop(paste0("Constraint ", names(con)[i], " must be a numeric value.\n"))
    }
    if((!(var[i] %in% Traits)) && con[[i]]<0){
      stop(paste0("Constraint ", names(con)[i], " must be a positive numeric value.\n"))
    }
  }
  
  condf <- data.frame(
              dir   = setNames(c("<=", ">=", "=="), c("ub", "lb", "eq"))[dir], 
              var   = var, 
              val   = unlist(con),
              isLin = var %in% Traits,
              row.names = var,
              stringsAsFactors = FALSE)
  
  return(condf)
}