
# Solving namespace issues:
#--------------------------------------
#' @importFrom acs acs.fetch is.geo.set
#' @importFrom data.table data.table is.data.table := setnames
#' @useDynLib synthACS
#' @importFrom Rcpp sourceCpp
NULL
#--------------------------------------

##---------------------------------------------------------
## CHECK IF OBJECT IS MEMBER OF CLASS
##---------------------------------------------------------

#' @title Check macroACS class
#' @description Function that checks if the target object is a \code{macroACS} object.
#' @param x any R object.
#' @return Returns \code{TRUE} if its argument has class "macroACS" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.macroACS <- function(x) {
  inherits(x, "macroACS")
}

#' @title Check micro_synthetic class
#' @description Function that checks if the target object is a \code{micro_synthetic} object.
#' @param x any R object.
#' @return Returns \code{TRUE} if its argument has class "micro_synthetic" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.micro_synthetic <- function(x) {
  inherits(x, "micro_synthetic")
}

#' @title Check macro_micro class
#' @description Function that checks if the target object is a \code{macro_micro} object.
#' @param x any R object.
#' @return Returns \code{TRUE} if its argument has class "macro_micro" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.macro_micro <- function(x) {
  inherits(x, "macro_micro")
}


#' @title Check synthACS class
#' @description Function that checks if the target object is a \code{synthACS} object.
#' @param x any R object.
#' @return Returns \code{TRUE} if its argument has class "synthACS" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.synthACS <- function(x) {
  inherits(x, "synthACS")
}

#' @title Check smsm_set class
#' @description Function that checks if the target object is a \code{smsm_set} object.
#' @param x any R object.
#' @return Returns \code{TRUE} if its argument has class "macroACS" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.smsm_set <- function(x) {
  inherits(x, "smsm_set")
}

##---------------------------------------------------------
## Generics for class macroACS
##---------------------------------------------------------
# This set of generic fetches "estimate" or "standard error" data for each of the datasets
# contained in the 'macroACS' class. 
##---------------------------------------------------------

####  helper functions ####
# 1. save duplication of typing on rowmatching
get_rowmatch <- function(string, symbol_table) {
  ns <- length(string)
  out <- vector("list", length= ns)
  for (i in 1:ns) {
    out[[i]] <- grep(paste0("*", string[i], " *"), symbol_table)
  }
  return(Reduce("c", out))
}

# 2. save duplication on input validation
validate_get_inputs <- function(acs, geography, dataset= c("estimate", "st.err")) {
  ## check inputs
  dataset <- match.arg(dataset, several.ok= FALSE)
  if (!is.macroACS(acs)) stop("acs must be of class macroACS.")
  if (!is.character(geography)) {
    stop("geography must be specified as a character vector.")
  } else {
    if (geography != "*" && any(nchar(geography) < 4)) 
      stop("Please specify at least 4 characters for geography.")
  }
  # if no error, okay
}

#' @title Get Aggregate Data Specified Geography
#' @description Gets aggregate, macro, data, either estimate or standard error, for a specified geography
#' and specified dataset.
#' @param acs An object of class \code{"macroACS"}.
#' @param geography A character vector allowing string matching via \code{\link[base]{grep}} to 
#' a set of specified geographies. All values may be specified by \code{"*"}.
#' @param dataset Either \code{"estimate"} or \code{"st.err"}. Do you want data on estimated 
#' population counts or estimated standard errors?
#' @param choice A character vector specifying the name of one of the datasets in \code{acs} 
#' @export
fetch_data <- function(acs, geography, dataset= c("estimate", "st.err"), 
                       choice= NULL) {
  UseMethod("fetch_data", acs) 
}

#' @export
fetch_data.macroACS <- function(acs, geography, dataset= c("estimate", "st.err"), 
                                choice= NULL) {
  dataset <- match.arg(dataset, several.ok= FALSE)
  choice <- match.arg(choice, choices= names(acs[["estimates"]]), several.ok= FALSE)
  
  ## check inputs
  validate_get_inputs(acs, geography, dataset)
  ## execute return
  if (length(geography) > 1L || geography != "*") {
    rowid <- get_rowmatch(geography, acs$geography$NAME)
    if (dataset == "estimate") { return(acs$estimates[[choice]][rowid, ]) } 
    else if (dataset == "st.err") { return(acs$standard_error[[choice]][rowid, ]) }
    else { stop("input 'dataset' is not valid.") }
  } else {
    if (dataset == "estimate") { return(acs$estimates[[choice]]) } 
    else if (dataset == "st.err") { return(acs$standard_error[[choice]]) }
    else { stop("input 'dataset' is not valid.") }
  }
}

##---------------------------------------------------------
## generics for getting the span, endyear and geography -- macroACS
##---------------------------------------------------------
# @title Constructor function for the "macroACS" class
# @description Constructor function for the "macroACS" class
# @param endyear An integer specifying the eendyear
# @param span An integer specifying the span of the data collection period. One of \code{c(1,3,5)}
# @param estimates a \code{list} of \code{data.frame}s
# @param standard_error a \code{list} of \code{data.frame}s
# @param geography A \code{data.frame} specifying the geography to which the data corrresponds
# @param geo_title An object of class 'geo'.
new_macroACS <- function(endyear, span, estimates, standard_error, geography,
                         geo_title) {
  # validate inputs
  if (! span %in% c(1,3,5)) stop("The ACS API only supports data spans of 1, 3, and 5 years.")
  if (endyear %% 1 != 0 | endyear < 2009) stop("endyear must be an integer >= 2009 (when ACS data begins).")
  if (!is.data.frame(geography)) stop("geography must be a data.frame")
  #if (!acs::is.geo(geo_title[[1]])) stop("geo_title must be a 'geo' object.")
  if (!is.list(estimates) & all(!unlist(lapply(estimates, is.data.frame))))
    stop("estimates must be a list of data.frames.")
  if (!is.list(standard_error) & all(!unlist(lapply(standard_error, is.data.frame))))
    stop("standard_error must be a list of data.frames.")
  
  # create new macroACS object
  new_macro <- list(endyear= endyear, span= span,
                    estimates= estimates,
                    standard_error= standard_error,
                    geography= geography,
                    geo_title= geo_title)
  class(new_macro) <- "macroACS"
  return(new_macro)
}


#' @title Get the span from a "macroACS" object.
#' @description Get the data collection span from a "macroACS" object
#' @param acs An object of class \code{"macroACS"}.
#' @export
get_span <- function(acs) {
  UseMethod("get_span", acs)
}

#' @export
get_span.macroACS <- function(acs) {
  return(acs$span)
}

#' @title Get the endyear from a "macroACS" object.
#' @description Get the data collection endyear from a "macroACS" object
#' @param acs An object of class \code{"macroACS"}.
#' @export
get_endyear <- function(acs) {
  UseMethod("get_endyear", acs)
}

#' @export
get_endyear.macroACS <- function(acs) {
  return(acs$endyear)
}

#' @title Get the geography title from a "macroACS" object.
#' @description Get the summary information of the geography selected from a "macroACS" object
#' @param acs An object of class \code{"macroACS"}.
#' @export
get_geography <- function(acs) {
  UseMethod("get_geography", acs)
}

#' @export
get_geography.macroACS <- function(acs) {
  return(acs$geo_title)
}

#' @title Get dataset names from a "macroACS" object.
#' @description Get the names of the datasets in a given "macroACS" object. 
#' @param acs An object of class \code{"macroACS"}.
#' @seealso \code{\link{fetch_data}}
#' @export
get_dataset_names <- function(acs) {
  UseMethod("get_dataset_names", acs)
}

#' @export
get_dataset_names.macroACS <- function(acs) {
  return(names(acs$estimates))
}

#' @title Split a "macroACS" object
#' @description Split a "macroACS" object into subsets. This may be helpful for users who have
#' limited memory available on their machines before proceding to derive sample synthetic micro data.
#' @param acs An object of class \code{"macroACS"}.
#' @param n_splits An integer for the number of splits you wish to create.
#' @seealso \code{\link{derive_synth_datasets}}
#' @export
split <- function(acs, n_splits) {
  UseMethod("split", acs)
}

#' @export
split.macroACS <- function(acs, n_splits) {
  # keep meta data
  sp <- get_span(acs)
  ey <- get_endyear(acs)
  orig_geo <- acs$geography
  geo_title <- NULL
  
  # split
  nx <- nrow(acs$geography)
  split_idx <- parallel::splitIndices(nx, ncl= n_splits)
  
  split_macroACS <- vector("list", length= n_splits)
  for (i in 1:n_splits) {
    geo <- orig_geo[ split_idx[[i]], ]
    est <- lapply(acs$estimates, function(l, idx) {return(l[idx, ])}, idx= split_idx[[i]])
    se <- lapply(acs$standard_error, function(l, idx) {return(l[idx, ])}, idx= split_idx[[i]])
    
    split_macroACS[[i]] <- new_macroACS(endyear= ey, span= sp, 
                                        estimates= est, standard_error= se,
                                        geography= geo,
                                        geo_title= geo_title)
  }
  return(split_macroACS)
}

#' @title Generate attribute vectors
#' @description Generate a list of attribute vectors for new synthetic attribute creation from a 
#' "macroACS" object.
#' @param acs An object of class \code{"macroACS"}.
#' @param choice A character vector specifying the name of one of the datasets in \code{acs} 
#' @seealso \code{\link{all_geog_synthetic_new_attribute}}, \code{\link{synthetic_new_attribute}}
#' @export
gen_attr_vectors <- function(acs, choice) {
  UseMethod("gen_attr_vectors", acs)
}

#' @export
gen_attr_vectors.macroACS <- function(acs, choice) {
  df <- fetch_data(acs, geography= "*", dataset= "estimate", choice= choice)
  return(disaggregate_md(df))
}

##---------------------------------------------------------
## GENERICS FOR CLASS "macro_micro" -- SPECIFICALLY ADDING CONSTRAINT LISTS
## FOR THE 10 DEFAULT SYNTHETIC ATTRIBUTES / VARIABLES
##---------------------------------------------------------

# helper function to save typing
# unsures that constraint population equals macro-population (for a geography)
# is used in synthACS-methods when method= "synthetic" (see below functions)
equal_constraint_populations <- function(constr_vec, geo_pop) {
  if (sum(constr_vec, na.rm=TRUE) == geo_pop) {return(constr_vec)}
  else if (sum(constr_vec, na.rm=TRUE) > geo_pop){ # equality from max
    idx <- which.max(constr_vec)
    constr_vec[idx] <- constr_vec[idx] + (geo_pop - sum(constr_vec, na.rm=TRUE))
    return(constr_vec)
  } else if (sum(constr_vec, na.rm=TRUE) < geo_pop) { # equality from min
    idx <- which.min(constr_vec)
    constr_vec[idx] <- constr_vec[idx] + (geo_pop - sum(constr_vec, na.rm=TRUE))
    return(constr_vec)
  }
}


#' @title Create gender constraint list to a set of geographies
#' @description Create a new gender constraint list to the mapping between a a set of macro datasets and a 
#' matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  g1 <- all_geog_constraint_gender(obj, "synthetic")
#'  g2 <- all_geog_constraint_gender(obj, "macro_table")
#' }
#' @export
all_geog_constraint_gender <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_gender", obj)
}

#' @export
all_geog_constraint_gender.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$gender, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- c("Male", "Female")
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      constr_vec <- l[[1]]$age_by_sex[2:3]
      names(constr_vec) <- c("Male", "Female")
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}


#' @title Create age constraint list to a set of geographies
#' @description Create a new age constraint list to the mapping between a a set of macro datasets and a 
#' matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  a1 <- all_geog_constraint_age(obj, "synthetic")
#'  a2 <- all_geog_constraint_age(obj, "macro_table")
#' }
#' @export
all_geog_constraint_age <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_age", obj)
}

#' @export
all_geog_constraint_age.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$age, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$age)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      constr_vec <- l[[1]]$age_by_sex[-c(1:3)]
      constr_vec <- apply(cbind(constr_vec[1:16], constr_vec[17:32]), 1, sum)
      names(constr_vec) <- levels(l[[2]]$age)
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}


#' @title Create marital status constraint list to a set of geographies
#' @description Create a new marital status constraint list to the mapping between a a set of macro 
#' datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  m1 <- all_geog_constraint_marital_status(obj, "synthetic")
#'  m2 <- all_geog_constraint_marital_status(obj, "macro_table")
#' }
#' @export
all_geog_constraint_marital_status <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_marital_status", obj)
}

#' @export
all_geog_constraint_marital_status.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$marital_status, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$marital_status)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) { 
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$marital_status[-1]
      nv_mar <- sum(v[which(grepl("nvr_married", names(v)))])
      mar <- v[which(grepl("married", names(v)))]
      mar <- sum(mar[which(!grepl("nvr_married", names(mar)))])
      mar_apart <- sum(v[which(grepl("mar_apart", names(v)))])
      div <- sum(v[which(grepl("divor", names(v)))])
      wid <- sum(v[which(grepl("widow", names(v)))])
      # under 15 not accounted for in base table; assume all never married 
      # (matches assumptions from synth_data_mar)
      pop_u15 <- sum(l$macro_constraints$age_by_sex[c(4,20)])
      
      constr_vec <- c("never_mar"= nv_mar + pop_u15, "married"= mar, "mar_apart"=mar_apart, 
                      "widowed"= wid, "divorced"= div)
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}


#' @title Create educational attainment constraint list to a set of geographies
#' @description Create a new educational attainment constraint list to the mapping between a a set 
#' of macro datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  e1 <- all_geog_constraint_edu(obj, "synthetic")
#'  e2 <- all_geog_constraint_edu(obj, "macro_table")
#' }
#' @export
all_geog_constraint_edu <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_edu", obj)
}

#' @export
all_geog_constraint_edu.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$edu_attain, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$edu_attain)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$edu[-c(1:3)]
      lt_hs <- sum(v[which(grepl("lt_hs", names(v)))])
      some_hs <- sum(v[which(grepl("some_hs", names(v)))])
      hs_grad <- sum(v[which(grepl("hs_grad", names(v)))])
      some_col <- sum(v[which(grepl("some_col", names(v)))])
      assoc_dec <- sum(v[which(grepl("assoc_dec", names(v)))])
      ba_deg <- sum(v[which(grepl("ba_deg", names(v)))])
      grad_deg <- sum(v[which(grepl("grad_deg", names(v)))])
      # under 18 not accounted for in base table 
      # assumptions from synth_data_edu:
        # under 15 == less than high school education
        # 15-17 == some high school education
      pop_u15 <- sum(l$macro_constraints$age_by_sex[c(4,20)])
      pop_15_17 <- sum(l$macro_constraints$age_by_sex[c(5,21)])
      
      constr_vec <- c("lt_hs"=  lt_hs + pop_u15, "some_hs"= some_hs + pop_15_17,
                      "hs_grad"= hs_grad, "some_col"= some_col, "assoc_dec"= assoc_dec,
                      "ba_deg"= ba_deg, "grad_deg"= grad_deg)
      
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}


#' @title Create employment status constraint list to a set of geographies
#' @description Create a new employment status constraint list to the mapping between a a set 
#' of macro datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  e1 <- all_geog_constraint_employment(obj, "synthetic")
#'  e2 <- all_geog_constraint_employment(obj, "macro_table")
#' }
#' @export
all_geog_constraint_employment <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_employment", obj)
}

#' @export
all_geog_constraint_employment.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$emp_status, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$emp_status)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$emp_status[-1]
      employed <- sum(v[which(grepl("employed", names(v)))])
      unemp <- sum(v[which(grepl("unemp", names(v)))])
      no_labor <- sum(v[which(grepl("no_labor", names(v)))])
      # under 16 not accounted for in base table 
      # assumptions from synth_data_emp:
      # assumptions: under 15 == not in labor force
      # assumptions: 15-17 best represented by 16-19 (will ignore 16 yr olds and let pop-match fix it)
      pop_u15 <- sum(l$macro_constraints$age_by_sex[c(4,20)])
      
      constr_vec <- c("not_in_labor_force"= no_labor + pop_u15, "employed"= employed, "unemployed"= unemp)
      
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}


#' @title Create nativity status constraint list to a set of geographies
#' @description Create a new nativity status constraint list to the mapping between a a set 
#' of macro datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  n1 <- all_geog_constraint_nativity(obj, "synthetic")
#'  n2 <- all_geog_constraint_nativity(obj, "macro_table")
#' }
#' @export
all_geog_constraint_nativity <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_nativity", obj)
}

#' @export
all_geog_constraint_nativity.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$nativity, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$nativity)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$nativity[-c(1:10)]
      born_other_state <- sum(v[which(grepl("born_out_state", names(v)))])
      born_out_us <- sum(v[which(grepl("born_st_res", names(v)))])
      born_state_residence <- sum(v[which(grepl("no_labor", names(v)))])
      foreigner <- sum(v[which(grepl("foreigner", names(v)))])
      
      constr_vec <- c("born_other_state"= born_other_state, "born_out_us"= born_out_us, 
                      "born_state_residence"= born_state_residence, "foreigner"= foreigner)
      
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}

#' @title Create poverty status constraint list to a set of geographies
#' @description Create a new poverty status constraint list to the mapping between a a set 
#' of macro datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  p1 <- all_geog_constraint_poverty(obj, "synthetic")
#'  p2 <- all_geog_constraint_poverty(obj, "macro_table")
#' }
#' @export
all_geog_constraint_poverty <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_poverty", obj)
}

#' @export
all_geog_constraint_poverty.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$pov_status, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$pov_status)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$pov_status1[-1]
      at_above <- sum(v[which(grepl("gt_eq_pov", names(v)))])
      below <- sum(v[which(grepl("lt_pov", names(v)))])
      # synth_pov: (individuals 15 years and over for whom poverty status is determined)
      # proportionally allocate u15
      pop_u15 <- sum(l$macro_constraints$age_by_sex[c(4,20)])
      
      constr_vec <- c("at_above_pov_level"= at_above + round(pop_u15 * at_above / (at_above + below) ,0), 
                      "below_pov_level"= below + round(pop_u15 * below / (at_above + below) ,0))
      
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}

#' @title Create geographic mobility constraint list to a set of geographies
#' @description Create a new geographic mobility constraint list to the mapping between a a set 
#' of macro datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  gm1 <- all_geog_constraint_geog_mob(obj, "synthetic")
#'  gm2 <- all_geog_constraint_geog_mob(obj, "macro_table")
#' }
#' @export
all_geog_constraint_geog_mob <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_geog_mob", obj)
}

#' @export
all_geog_constraint_geog_mob.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$geog_mobility, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$geog_mobility)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$geo_mob_edu
      CNT <- v[1]
      same_house <- sum(v[which(grepl("same_house_all", names(v)))])
      same_cty   <- sum(v[which(grepl("same_cnty_all", names(v)))])
      same_st    <- sum(v[which(grepl("same_st_all", names(v)))])
      diff_st    <- sum(v[which(grepl("diff_st_all", names(v)))])
      abroad     <- sum(v[which(grepl("abroad_all", names(v)))])
      
      # v <- l[[1]]$edu[4:10] + v <- l[[1]]$edu[39:45] #18-24 only
      # lt_hs <- sum(v[which(grepl("lt_hs", names(v)))])
      # some_hs <- sum(v[which(grepl("some_hs", names(v)))])
      # hs_grad <- sum(v[which(grepl("hs_grad", names(v)))])
      # some_col <- sum(v[which(grepl("some_col", names(v)))])
      # assoc_dec <- sum(v[which(grepl("assoc_dec", names(v)))])
      # ba_deg <- sum(v[which(grepl("ba_deg", names(v)))])
      # grad_deg <- sum(v[which(grepl("grad_deg", names(v)))])
      
      # synth_geomob: (Universe:  Population 25 years and over in the United States)
      # ALLOCATE ALL PROPORTIONALLY BASED ON EDU
      # u15 = less than high school 
      # 15_17 = some HS 
      # 18_24 = have data; but just proportionally allocate equally
      pop_u15 <- sum(l$macro_constraints$age_by_sex[c(4,20)])
      pop_15_17 <- sum(l$macro_constraints$age_by_sex[c(5,21)])
      pop_18_24 <- sum(l$macro_constraints$age_by_sex[c(6,22)])
      
      constr_vec <- c("diff state"= diff_st + round(pop_u15 * v[21] / same_house, 0) + 
                        round(pop_15_17 * v[22] / diff_st, 0) + round(pop_18_24 * diff_st / CNT, 0), 
                      "moved from abroad"= abroad + round(pop_u15 * v[27] / same_house, 0) + 
                        round(pop_15_17 * v[28] / abroad, 0) + round(pop_18_24 * abroad / CNT, 0), 
                      "same county"= same_cty + round(pop_u15 * v[9] / same_house, 0) + 
                        round(pop_15_17 * v[10] / same_cty, 0) + round(pop_18_24 * same_cty / CNT, 0),
                      "same house"= same_house + round(pop_u15 * v[3] / same_house, 0) + 
                        round(pop_15_17 * v[4] / same_house, 0) + round(pop_18_24 * same_house / CNT, 0), 
                      "same state"= same_st + round(pop_u15 * v[15] / same_house, 0) + 
                        round(pop_15_17 * v[16] / same_st, 0) + round(pop_18_24 * same_st / CNT, 0)) 
      
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$geog_mobility)
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}

#' @title Create individual income constraint list to a set of geographies
#' @description Create a new individual income constraint list to the mapping between a a set 
#' of macro datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  i1 <- all_geog_constraint_income(obj, "synthetic")
#'  i2 <- all_geog_constraint_income(obj, "macro_table")
#' }
#' @export
all_geog_constraint_income <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_income", obj)
}

#' @export
all_geog_constraint_income.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$ind_income, sum) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$ind_income)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$by_inc_12mo
      # synth_inc: (Universe:  Population 15 years and over in the United States)
      # proportionally allocate u15 to under 35k
      pop_u15 <- sum(l$macro_constraints$age_by_sex[c(4,20)])
      
      constr_vec <- c("1_lt10k"= v[3] + round(pop_u15 * v[3] / sum(v[2:6]), 0),
                      "10k_lt15k"= v[4] + round(pop_u15 * v[4] / sum(v[2:6]), 0),
                      "15k_lt25k"= v[5] + round(pop_u15 * v[5] / sum(v[2:6]), 0),
                      "25k_lt35k"= v[6] + round(pop_u15 * v[6] / sum(v[2:6]), 0),
                      "35k_lt50k"= v[7],
                      "50k_lt65k"= v[8],
                      "65k_lt75k"= v[9],
                      "gt75k"= v[10],
                      "no_income"= v[2] + round(pop_u15 * v[2] / sum(v[2:6]) , 0))
      
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$ind_income)
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}

#' @title Create race constraint list to a set of geographies
#' @description Create a new race constraint list to the mapping between a a set 
#' of macro datasets and a matching set of micro dataset (supplied as class 'synthACS').
#' @param obj An object of class \code{"synthACS"}.
#' @param method One of \code{c("synthetic", "macro.table")}. Specifying \code{"synthetic"} indicates
#' that constraints are built by marginalizing the synthetic micro datasets. Specifying 
#' \code{"macro.table"} indicates that the constraints are build from the data in the base ACS tables.
#' @seealso \code{\link{all_geogs_add_constraint}}
#' @examples \dontrun{
#'  # assumes that obj of class 'synthACS' already exists in your environment
#'  r1 <- all_geog_constraint_race(obj, "synthetic")
#'  r2 <- all_geog_constraint_race(obj, "macro_table")
#' }
#' @export
all_geog_constraint_race <- function(obj, method= c("synthetic", "macro.table")) {
  UseMethod("all_geog_constraint_race", obj)
}

#' @export
all_geog_constraint_race.synthACS <- function(obj, method= c("synthetic", "macro.table")) {
  method= match.arg(method, several.ok= FALSE)
  
  # A - synthetic
  if (method == "synthetic") {
    constraint <- lapply(obj, function(l) {
      # marginalize probability vector by attribute, normalize by total_pop
      constr_vec <- as.numeric(round(tapply(l[[2]]$p, l[[2]]$race, sum, na.rm=T) * l[[1]]$age_by_sex[1], 0))
      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$race)
      # check that population matches macro pop and return
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  # B - macro.table
  else if (method == "macro.table") {
    constraint <- lapply(obj, function(l) {
      # done via name pattern matching and summing across age brackets
      v <- l[[1]]$pop_by_race
      # race includes breakout by white alone v hispanic white... need to reduce counts to match
      # total (v[1] vs sum(v[3:10]))
      white   <- round(sum(v[which(grepl("total_white_alone", names(v)))]) * v[1] / sum(v[3:10]), 0)
      black   <- round(sum(v[which(grepl("total_black", names(v)))]) * v[1] / sum(v[3:10]), 0)
      hisp    <- round(sum(v[which(grepl("total_hisp_lat", names(v)))]) * v[1] / sum(v[3:10]), 0)
      asian   <- round(sum(v[which(grepl("total_asian", names(v)))]) * v[1] / sum(v[3:10]), 0)
      native  <- round(sum(v[which(grepl("total_nat_amer", names(v)))]) * v[1] / sum(v[3:10]), 0)
      pac     <- round(sum(v[which(grepl("total_pac_isl", names(v)))]) * v[1] / sum(v[3:10]), 0)
      more2   <- round(sum(v[which(grepl("total_2p_races", names(v)))]) * v[1] / sum(v[3:10]), 0)

      constr_vec <- c("asian"= asian,
                      "black, afr amer"= black,
                      "hispanic, latino"= hisp,
                      "native amer"= native,
                      "pacific islander"= pac,
                      "two or more races"= more2,
                      "white alone"=  white)

      constr_vec <- ifelse(is.na(constr_vec), 0, constr_vec)
      names(constr_vec) <- levels(l[[2]]$race)
      return(equal_constraint_populations(constr_vec, l[[1]]$age_by_sex[1]))
    })
  }
  return(constraint)
}

##---------------------------------------------------------
## Generics for class micro_synthetic
##---------------------------------------------------------
#' @title Marginalize synthetic attributes
#' @description Marginalize, (ie- reduce in number), attributes of a synthetic dataset of class
#' 'micro_synthetic' or a list of synthetic datasets of class 'synthACS'. This is done
#' by marginalizing the joint distribution based on a set of specified attributes (see Arguments below).
#' @param obj An object of class \code{"micro_synthetic"}.
#' @param varlist A character vector of variable, or attribute, names in \code{obj}. 
#' @param marginalize_out Logical. Do you wish to *remove* the variables in \code{varlist} 
#' instead of keeping them? Defaults to \code{FALSE}
#' @examples {
#' # dummy data setup
#' set.seed(567L)
#' df <- data.frame(gender= factor(sample(c("male", "female"), size= 100, replace= TRUE)),
#'                  age= factor(sample(1:5, size= 100, replace= TRUE)),
#'                  pov= factor(sample(c("below poverty", "at above poverty"), 
#'                                    size= 100, replace= TRUE, prob= c(.15,.85))),
#'                  p= runif(100))
#' df$p <- df$p / sum(df$p)
#' class(df) <- c("data.frame", "micro_synthetic")
#' 
#' df2 <- marginalize_attr(df, varlist= "gender")
#' df3 <- marginalize_attr(df, varlist= c("gender", "age"))
#' df4 <- marginalize_attr(df, varlist= c("gender", "age"), marginalize_out= TRUE)
#' 
#' df_list <- replicate(10, df, simplify= FALSE)
#' dummy_list <- replicate(10, list(NULL), simplify= FALSE)
#' df_list <- mapply(function(a,b) {return(list(a, b))}, a= dummy_list, b= df_list, SIMPLIFY = FALSE)
#' class(df_list) <- c("list", "synthACS")
#' 
#' # run the function
#' df_list2 <- marginalize_attr(df_list, varlist= c("gender", "age"))
#' }
#' @export
#' 
marginalize_attr <- function(obj, varlist, marginalize_out= FALSE) {
  UseMethod("marginalize_attr", obj)
}

#' @export
marginalize_attr.micro_synthetic <- function(obj, varlist, marginalize_out= FALSE) {
  if (!is.character(varlist) || !is.vector(varlist)) stop("varlist must be specified as a character vector.")
  if (!all(sapply(varlist, function(l) exists(l, as.environment(obj)))))
    stop("at least one variable in varlist is not in obj.")
  if (!("p" %in% names(obj))) {
    p_idx <- which(apply(obj, 2, is.numeric))
    if (length(p_idx) > 1L) stop("obj appears to have more than 1 probability vector.")
    p_name <- names(obj)[p_idx]
    names(obj)[p_idx] <- "p"
  }

  data.table::setDT(obj)
  if (!marginalize_out) {
    obj <- obj[,sum(p), by= varlist]
    data.table::setnames(obj, "V1", ifelse(exists("p_name"), p_name, "p"))
    class(obj) <- c("data.table", "data.frame", "micro_synthetic")
    return(obj)
  } else {
    p <- NULL # needed for R CMD Check -- see stackoverflow./com/questions/8096313
    vlist2 <- names(obj)[which(!names(obj) %in% c(varlist, "p"))]
    obj <- obj[,sum(p), by= vlist2]
    data.table::setnames(obj, "V1", ifelse(exists("p_name"), p_name, "p"))
    class(obj) <- c("data.table", "data.frame", "micro_synthetic")
    return(obj)
  }
}

#' @export
marginalize_attr.synthACS <- function(obj, varlist, marginalize_out= FALSE) {
  obj <- lapply(obj, function(l, v, m) {
    l[[2]] <- marginalize_attr.micro_synthetic(l[[2]], v, m)
    return(l)
  }, v= varlist, m= marginalize_out)
  class(obj) <- c("synthACS", "list")
  return(obj)
}



##---------------------------------------------------------
## Generics for class smsm_set
##---------------------------------------------------------
#' @title Summarizing SMSM fits
#' @description \code{summary} method for class 'smsm_set'. 
#' @param object An object of class \code{'smsm_set'}, typically a result of call to 
#' \code{\link{all_geog_optimize_microdata}}
#' @param ... additional arguments affecting the summary produced.
#' @export
summary.smsm_set <- function(object, ...) {
  
  tae_pct <- unlist(object$tae) / unlist(lapply(object$best_fit, nrow)) / object$D
  tae_q <- round(stats::quantile(tae_pct, c(0,.25,.5,.75,.9,.95)), 6)
  names(tae_q) <- c("0%", "25%", "50%", "75%", "90%", "95%")
  n_early <- sum(unlist(object$iter) < object$max_iter)
  
  cat("\n Call: \n", paste(deparse(object$call), collapse= "\n"),
      "\n \n Seed: ", object$seed,
      "\n n-Constraints: ", object$D,
      "\n \n Maximum Iterations: ", object$max_iter,
      "\n %-Early Stop: ", round(n_early / length(object$iter), 4),
      "\n \n Mean %-TAE: ", round(mean(tae_pct), 6),
      "\n Median %-TAE: ", round(stats::median(tae_pct), 6),
      "\n Max %-TAE: ", round(max(tae_pct), 6),
      "\n %-TAE quantiles: \n ")
  print(tae_q)
}

#' @title Plot simulated annealing path
#' @description Plot the path TAE in the simulated annealing algorithm for a given geography
#' @param object An object of class \code{'smsm_set'}, typically a result of call to 
#' \code{\link{all_geog_optimize_microdata}}
#' @param geography A string allowing string matching via \code{\link[base]{grep}} to 
#' a specified geography.
#' @param ... additional arguments passed to other methods
#' @export
plot_TAEpath <- function(object, geography,  ...) {
  UseMethod("plot_TAEpath", object)
}

#' @export
plot_TAEpath.smsm_set <- function(object, geography, ...) {
  
  if (length(geography) != 1) stop("Please specify a single geography")
  idx <- get_rowmatch(geography, names(object$tae))
  if (length(idx) > 1) stop("geography specification returns multiple results. Please be more specific.")
  tae_path <- object$tae_paths[[ idx ]]
  if (!is.matrix(tae_path)) {
    tae_path <- matrix(tae_path, ncol= 2)
  }
  
  y_min <- max(min(tae_path) - 10, 0)
  y_max <- round(max(tae_path) * 1.1, 0)
  graphics::plot(x= 1:nrow(tae_path), y= tae_path[,2], type= "b", col= "black", pch= 18,
       xlab= "Iteration",ylab= "TAE", ylim= c(y_min, y_max), ...)
  graphics::points(x= 1:nrow(tae_path), y= tae_path[,1], pch= 1, col= "red")
  graphics::title(main= paste(names(object$best_fit)[idx], "\nPath of Simulated Annealing Fit"))
  graphics::legend(x= nrow(tae_path) * .76, y = y_max * .95, legend= c("Proposal", "Current"), title= "TAE",
         pch= c(1, 18), col= c("red", "black"), cex= 0.8, bty= "n")
}

# @description save some duplicate typing in checking global parameter matches in combine_smsm() below
# @param l A list of objects of class 'smsm_set'.
# @param param_name a string specifying one of the elements of an 'smsm_set' object
# @param warning_msg A string specifying the warning messsage if equality is not found.
global_param_check <- function(l, param_name, warning_msg) {
  param <- unlist(lapply(l, "[[", param_name))
  if (!all(param == param[1])) warning(warning_msg)
  return(param[1])
}

#' @title Combine separate SMSM optimizations
#' @description Combine objects of class "smsm_set" into a single object of class "smsm_set"
#' @param ... A list of objects of class 'smsm_set'.
#' @seealso \code{\link[synthACS]{split}}, \code{\link[synthACS]{all_geog_optimize_microdata}}
#' @examples \dontrun{
#'  combined <- combine_smsm(smsm1, smsm2, smsm3)
#' }
#' @export
combine_smsm <- function(...) {
   smsm <- list(...)
   if (!all(unlist(lapply(smsm, is.smsm_set)))) stop("All items supplied via '...' must be of class 'smsm_set'.")
   
   # 01. check if global parameter values are the same for all objects
   max_iter <- global_param_check(smsm, "max_iter",
                 warning_msg= "'smsm_set' objects do not have the equal max_iter values. Setting max_iter to first value.")
     
   D <- global_param_check(smsm, "D",
                 warning_msg= "'smsm_set' objects do not have the equal D values. Setting D to first value.")
   
   p_accept <- global_param_check(smsm, "p_accept",
                 warning_msg= "'smsm_set' objects do not have the equal p_accept values. Setting p_accept to first value.")
   seed <- unlist(lapply(smsm, "[[", "seed"))
   if (!all(seed == seed[1])) warning("'smsm_set' objects do not have equal seed values. Returning all seed values.")
   else {seed <- seed[1]}
   # [remove, message is confusing] message("call will be set to NULL.")
   
   # 02. create return structure / combine list of smsm_set objects
   ret <- list(best_fit=  do.call("c", lapply(smsm, "[[", "best_fit")),
               tae=       do.call("c", lapply(smsm, "[[", "tae")),
               call=      NULL,
               p_accept=  p_accept,
               iter=      do.call("c", lapply(smsm, "[[", "iter")),
               max_iter=  max_iter,
               tae_paths= do.call("c", lapply(smsm, "[[", "tae_paths")),
               seed=      seed,
               D=         D)
   class(ret) <- "smsm_set"
   return(ret)
}

#' @title Extract best fit for a specified geogrpahy from an 'smsm_set' object
#' @description Extract the best fit micro population (resulting from the simulated annealing 
#' algorithm) for a given geography.
#' @param obj An object of class \code{'smsm_set'}, typically a result of call to 
#' \code{\link{all_geog_optimize_microdata}}
#' @param geography A string allowing string matching via \code{\link[base]{grep}} to 
#' a specified geography.
#' @export
get_best_fit <- function(obj, geography) {
  UseMethod("get_best_fit", obj)
}

#' @export
get_best_fit.smsm_set <- function(obj, geography) {
  if (length(geography) != 1) stop("Please specify a single geography")
  idx <- get_rowmatch(geography, names(obj$best_fit))
  if (length(idx) > 1) stop("geography specification returns multiple results. Please be more specific.")
  
  return(obj$best_fit[[ idx ]])
}

#' @title Extract the final TAE for a specified geogrpahy from an 'smsm_set' object
#' @description Extract the final TAE (resulting from the simulated annealing 
#' algorithm) for a given geography.
#' @param obj An object of class \code{'smsm_set'}, typically a result of call to 
#' \code{\link{all_geog_optimize_microdata}}
#' @param geography A string allowing string matching via \code{\link[base]{grep}} to 
#' a specified geography.
#' @export
get_final_tae <- function(obj, geography) {
  UseMethod("get_final_tae", obj)
}

#' @export
get_final_tae.smsm_set <- function(obj, geography) {
  if (length(geography) != 1) stop("Please specify a single geography")
  idx <- get_rowmatch(geography, names(obj$best_fit))
  if (length(idx) > 1) stop("geography specification returns multiple results. Please be more specific.")
  
  return(obj$tae[[ idx ]])
}

## Z-statistics???

