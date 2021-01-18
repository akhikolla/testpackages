
# @title Pattern matching marginally missing to complete
# @description Performs pattern matching from marginally missing (z_Os_y) to potential
# complete observations (x_possible). Uses RSQLite instead of C++ selection search (marg_comp_compare)
# @param z_Os_y the marginally missing observations
# @param x_possible the set of possible fully observed observations
search_z_Os_y <- function(z_Os_y, x_possible) {
  if (is.null(names(x_possible)) | is.null(names(z_Os_y))) stop("both arguments / parameters must have names.")
  if (any(is.na(names(x_possible))) | any(is.na(names(z_Os_y)))) 
    stop("names may not be NA.")
  
  ## setup output list and SQLite database
  search_out <- vector("list", length= nrow(z_Os_y))
  
  x_p <- RSQLite::dbConnect(RSQLite::SQLite(), ":memory:")
  x_possible2 <- x_possible
  x_possible2$rownames <- 1:nrow(x_possible2)
  RSQLite::dbWriteTable(x_p, "x_possible", x_possible2)
  RSQLite::dbGetQuery(x_p, paste0("create index idx on x_possible (", paste(names(x_possible),collapse=", "), ")"))
  
  ## run
  for (s in 1:nrow(z_Os_y)) {
    search_out[[s]] <- as.integer(RSQLite::dbGetQuery(x_p, create_search_query(z_Os_y, 
                          z_Os_y[s, -ncol(z_Os_y)], var_names= names(x_possible)))$rownames)
  }
  return(search_out)
}

# helper fucntion for creating queries in search_z_Os_y
create_search_query <- function(df, row, var_names) {
  idx <- which(!is.na(row))
  n_na <- length(idx)
  q <- paste0("select rownames from x_possible where ", var_names[ idx[1] ], "= '", 
              get_level_text(get(var_names[ idx[1] ], as.environment(df)) , as.integer(row[ idx[1] ])),"'")
  if (n_na == 1) return(q)
  else {
    for (i in 2:n_na) {
      q <- paste0(q, "and ", var_names[idx[i]], "= '", 
                  get_level_text(get(var_names[ idx[i] ], as.environment(df)) , as.integer(row[ idx[i] ])), "'")
    }
    return(q)
  }
}



# nnodes <- min(nrow(dat2), parallel::detectCores() - leave_cores)
# if (grepl("Windows", utils::sessionInfo()$running)) {cl <- parallel::makeCluster(nnodes, type= "PSOCK")}
# else {cl <- parallel::makeCluster(nnodes, type= "FORK")}
# 
# comp_ind <- parallel::clusterApply(cl, x= splitRows(z_Os_y, nnodes), fun= search_z_Os_y,
#                                      x_possible= x_possible)
# 
# parallel::stopCluster(cl)


count_sumStats <- function(x_possible, dat, hasNA= c("no", "count.obs", "count.miss")) {
  # parameter checking
  hasNA <- match.arg(hasNA, several.ok= FALSE)
  if (ncol(dat) != ncol(x_possible)) stop("ncol(dat) and ncol(enum_list) must match.")
  
  ## 0. Pre-processing: convert factors to integers
  dat <- data.frame(do.call("cbind", lapply(dat, fact_to_int)))
  x_possible <- data.frame(do.call("cbind", lapply(x_possible, fact_to_int)))

  # setup database
  x_p <- RSQLite::dbConnect(RSQLite::SQLite(), ":memory:")
  RSQLite::dbWriteTable(x_p, "dat", dat)
  RSQLite::dbGetQuery(x_p, paste0("create index idx on dat (", paste(names(dat),collapse=", "), ")"))
  nm <- names(x_possible)
  x_possible$counts <- 0

  for (i in 1:nrow(x_possible)) {
    x_possible$counts[i] <-  RSQLite::dbGetQuery(x_p, create_count_query(x_possible, x_possible[i,], nm))$cnt
  }
  return(x_possible[!is.na(x_possible$counts) & x_possible$counts > 0,])
}

create_count_query <- function(df, row, var_names) {
  nx <- length(row)-1
  q <- paste0("select count(*) as cnt from dat where ", var_names[1], "= '", row[1],"'")
  if (nx == 1) return(q)
  else {
    for (i in 2:nx) {
      q <- paste0(q, "and ", var_names[i], "= '", row[i], "'")
    }
    return(q)
  }
}

