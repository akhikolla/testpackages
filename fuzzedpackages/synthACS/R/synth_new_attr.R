
#' @title Add a new attribute to a synthetic_micro dataset
#' @description Add a new attribute to a synthetic_micro dataset using conditional relationships
#' between the new attribute and existing attributes (eg. wage rate conditioned on age and education 
#' level).  
#' @section Details:
#' New synthetic variables are introduced to the existing data via conditional probability. Similar 
#' to \code{\link{derive_synth_datasets}}, the goal with this function is to generate a joint 
#' probability distribution for an attribute vector; and, to create synthetic individuals from 
#' this distribution. Although no limit is placed on the number of variables on which to condition, 
#' in practice, data rarely exists which allows more than two or three conditioning variables. Other 
#' variables are assumed to be independent from the new attribute. 
#' 
#' ** There are four different types of conditional/marginal probability models which may be considered
#' for a given new attribute:
#'  (1) Independence: it is assumed that each of the variables is independent of the others
#'  (2) Pairwise conditional independence: it is assumed that attributes are related to 
#'  only one other attribute and independent of all others.
#'  (3) Conditional independence: Attributes can be depedent on some subset of other attributes and 
#'  independent of the rest.
#'  (4) In the most general case, all attributes are jointly interrelated.
#'  
#' Conditioning is implemented via symbol-tables (\code{sym_tbl}) to ensure accurate matching between
#' conditioning variables, new attribute levels, and new attribute probabilities. The symbol table
#' is constructed such that the key in the symbol-table's key-value pair is the specific values for 
#' the set of conditioning variables. This key is the first N columns of \code{sym_tbl}.  A 
#' recursive approach is employed to conditionally partition \code{sym_tbl}. In this sense, the 
#' *order* in which the conditional variables are supplied matters.
#' 
#' The value is final 2 columns of \code{sym_tbl} which are a pair of (A) either counts or percentages 
#' used to specify the probability for the new attribute and (B) the level that the new attribute takes on.
#' 
#' @param df An R object of class "synthetic_micro". 
#' @param prob_name A string specifying the column name of the \code{df} containing the
#' probabilities for each synthetic observation.
#' @param attr_name A string specifying the desired name of the new attribute to be added to the data.
#' @param conditional_vars An character vector specifying the existing variables, if any, on which 
#' the new attribute (variable) is to be conditioned on. Variables must be specified in order. 
#' Defaults to \code{NULL} ie- an unconditional new attribute.
#' @param sym_tbl sym_tbl A \code{data.frame} symbol table with N + 2 columns. The last two columns must be:
#' 1. A vector containing the new attribute counts or percentages; 2. is a vector of the new attribute 
#' levels. The first N columns must match the conditioning scheme imposed by the variables in 
#' \code{conditional_vars}. See details and examples. 
#' @return A new synthetic_micro dataset with class "synthetic_micro".
#' @examples {
#' set.seed(567L)
#' df <- data.frame(gender= factor(sample(c("male", "female"), size= 100, replace= TRUE)),
#'                 edu= factor(sample(c("LT_college", "BA_degree"), size= 100, replace= TRUE)),
#'                 p= runif(100))
#' df$p <- df$p / sum(df$p)
#' class(df) <- c("data.frame", "micro_synthetic")
#' ST <- data.frame(gender= c(rep("male", 3), rep("female", 3)),
#'                  attr_pct= c(0.1, 0.8, 0.1, 0.05, 0.7, 0.25),
#'                  levels= rep(c("low", "middle", "high"), 2))
#' df2 <- synthetic_new_attribute(df, prob_name= "p", attr_name= "SES", conditional_vars= "gender",
#'          sym_tbl= ST)
#' 
#' ST2 <- data.frame(gender= c(rep("male", 3), rep("female", 6)),
#'                   edu= c(rep(NA, 3), rep(c("LT_college", "BA_degree"), each= 3)),
#'                   attr_pct= c(0.1, 0.8, 0.1, 10, 80, 10, 5, 70, 25),
#'                   levels= rep(c("low", "middle", "high"), 3))
#' df2 <- synthetic_new_attribute(df, prob_name= "p", attr_name= "SES",
#'          conditional_vars= c("gender", "edu"),
#'          sym_tbl= ST2)
#' }
#' @export
synthetic_new_attribute <- function(df, prob_name= "p",
                                    attr_name= "variable",
                                    conditional_vars= NULL,
                                    sym_tbl= NULL) {
  # 01. Error checking
  #------------------------------------
  n_st <- ncol(sym_tbl); n_cv <- length(conditional_vars)
  
  if (!is.micro_synthetic(df)) stop("Input appropriate df -- use class 'micro_synthetic'.")
  if (!is.character(prob_name) || length(prob_name) != 1L) stop("prob_name must be a string.")
  if (!exists(prob_name, as.environment(df))) stop("prob_name is not in df.")
  if (!is.character(attr_name) || length(attr_name) != 1L) stop ("attr_name must be specified as a character string")
  if (!is.numeric(sym_tbl[, n_st - 1])) 
    stop("The second to last column of sym_tbl must be numeric.")
  if (!is.null(conditional_vars)) {
    # error check conditional_vars
    if (!all(sapply(conditional_vars, is.character)))
      stop("conditional_vars must be specified as strings.")
    if (!all(sapply(conditional_vars, function(l) exists(l, as.environment(df)))))
      stop("at least one conditional_var is not in df.")
    # error check ht_list, must be co-specified.
    if (is.null(sym_tbl)) stop("sym_tbl must be specified")
    else {
      if (!is.data.frame(sym_tbl) | (n_st - 2 != n_cv))
        stop("sym_tbl must contain conditioning variables of equal length as conditional_vars.")
      if (is.data.table(sym_tbl)) { class(sym_tbl) <- "data.frame" } 
      if ( !all(unlist(lapply(sym_tbl[, 1:(n_st - 2)], 
            function(l) { sapply(l, function(i) is.character(i) || is.factor(i)) }))) )
        stop("all conditioning elements of sym_tbl must be strings or factors.")
      if (any(names(sym_tbl)[1:(n_st - 2)] != conditional_vars))
        stop("Variable names in sym_tbl must match conditional_vars.")
    }
  }
  
  # 02. Variable conditioning and Apply new synthetic variable
  #------------------------------------
  if (!is.null(conditional_vars)) {
    dat <- cond_var_split(df= df, prob_name= prob_name,
                          attr_name= attr_name,
                          conditional_vars= conditional_vars,
                          sym_tbl= sym_tbl)
  } else { # apply new attribute unconditionally
    n_lvl <- length(table(sym_tbl[, n_st]))
    dat <- replicate(n_lvl, df, simplify = FALSE)
    
    sym_tbl[, n_st - 1] <- sym_tbl[, n_st - 1] / sum(sym_tbl[, n_st - 1])
    sym_tbl <- base::split(sym_tbl, 1:nrow(sym_tbl))
    
    dat <- do.call("rbind", mapply(add_synth_attr_level, 
                                   dat= dat, prob_name= prob_name, attr_name= attr_name,
                                   attr= sym_tbl,
                                   SIMPLIFY = FALSE)) 
  }
  
  # 03. return
  #------------------------------------
  if  (!is.data.table(dat)) { dat <- dat[dat[prob_name] > 0,] }
  else { dat <- dat[get(prob_name, as.environment(dat)) > 0,] }
  
  if (!is.micro_synthetic(dat))   class(dat) <- c(class(dat), "micro_synthetic")
  return(dat)
}

# simple helper function to reduce typing.
split_df <- function(d, var) {
  base::split(d, get(var, as.environment(d)))
}



# @param df An R object of class "synthetic_micro". 
# @param prob_name A string specifying the column name of the \code{df} containing the
# probabilities for each synthetic observation.
# @param attr_name A string specifying the desired name of the new attribute to be added to the data.
# @param conditional_vars An character vector specifying the existing variables, if any, on which 
# the new attribute (variable) is to be conditioned on. Variables must be specified in order. 
# Defaults to \code{NULL} ie- an unconditional new attribute.
# @param sym_tbl A \code{data.frame} symbol table with N + 2 columns. The last two columns must be:
# 1. A vector containing the new attribute counts or percentages; 2. is a vector of the new attribute 
# levels. The first N columns must match the conditioning scheme imposed by 
# the variables in \code{conditional_vars}....
cond_var_split <- function(df, prob_name, attr_name= "variable", 
                           conditional_vars, sym_tbl) {
  if (nrow(df) < 1L) return(df)
  if (is.data.table(sym_tbl)) { class(sym_tbl) <- "data.frame" } 
  
  cv_n <- length(conditional_vars)
  st_n <- ncol(sym_tbl) - 2
  
  if (cv_n == 0 & st_n == 0) { # if lengths == 1, bottom of tree. Apply and return
    return(add_synth_attr(l= df, prob_name= prob_name, sym_tbl= sym_tbl, attr_name= attr_name))
  } else { 
    if (any(names(sym_tbl)[1:st_n] != conditional_vars))
      stop("Variable names in sym_tbl must match conditional_vars.") # quick error check
    
    # else: split data & ST conditionally, then recurse
    df <- split_df(df, conditional_vars[1])
    df <- df[which(unlist(lapply(df, nrow)) > 0)]
    nm_df <- names(df)
    df <- df[order(names(df))]
    
    if (!all(is.na(sym_tbl[,1]))) {
      sym_tbl <- base::split(sym_tbl[,-1], sym_tbl[,1])
      # validate ST only includes values from DF
      sym_tbl <- sym_tbl[ which(names(sym_tbl) %in% nm_df) ]
      sym_tbl <- sym_tbl[order(names(sym_tbl))]
      
      return(do.call("rbind", mapply(cond_var_split, 
             df= df, prob_name= prob_name, attr_name= attr_name, 
             conditional_vars= ifelse(cv_n == 1, replicate(st_n - 1, NULL), conditional_vars[-1]), 
             sym_tbl= sym_tbl, SIMPLIFY= FALSE)))
    } else {
      return(do.call("rbind", lapply(df, cond_var_split,
             prob_name= prob_name, attr_name= attr_name, 
             conditional_vars= conditional_vars[-1], 
             sym_tbl= sym_tbl[,-1])))
    }
  }
}



# @param l A data.frame -- a conditional subset of df 
# @param sym_tbl A \code{data.frame} symbol table with two columns: Column 1 is a vector containing the new 
# attribute counts or percentages; Column 2 is a vector of the new attribute levels.
# @param attr_name the name of the new variable 
# @ prob_name A string specifying the column name within \code{l} containing the
# probabilities for each synthetic observation.
add_synth_attr <- function(l, prob_name, sym_tbl, attr_name= "variable") {
  if (nrow(l) < 1L) return(l)
  if (is.data.table(sym_tbl)) { class(sym_tbl) <- "data.frame" } 
  
  if (ncol(sym_tbl) != 2) stop("incorrect dimensions for sym_tbl.")
  if (!(is.factor(sym_tbl[,2]) | is.character(sym_tbl[,2]))) stop("sym_tbl new-levels incorrectly specified.")
  if (!is.numeric(sym_tbl[,1])) 
    stop("sym_tbl attribute counts incorrectly specified. Must be specified as  either percentages summing to 1 or counts > 1.")
  # unit norm
  if (sum(sym_tbl[,1]) > 1) { sym_tbl[,1] <- sym_tbl[,1] / sum(sym_tbl[,1])} 
  
  # split by row
  sym_tbl <- base::split(sym_tbl, f= 1:nrow(sym_tbl))
  
  # replicate data and apply new levels/probabilities
  dat <- replicate(length(sym_tbl), l, simplify = FALSE)
  dat <- do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= prob_name, 
                                 attr= sym_tbl, attr_name= attr_name,
                                 SIMPLIFY = FALSE))
  
  return(dat)
}


# helper function for bottom of recursion -- apply new level to 
# smallest subset of data and update probabilities appropriately

# @title add new synthetic attribute to a dataset
# @param dat A current dataset on which to add an attribute
# @param prob_name A string specifying the name of the probability vector within \code{dat}
# @param attr A vector of length two. The first element must be the attribute probability and
# the second element must be the attribute level.
# @param attr_name A string specifying the name of the new attribute. (eg. "gender")
add_synth_attr_level <- function(dat, prob_name, attr, attr_name= "variable") {
  p <- get(prob_name, as.environment(dat))
  if (is.data.table(dat)) {
    d_temp <- dat[, which(names(dat) != prob_name), with= FALSE]
    d_temp[, `:=` (V_NEW= attr[[2]],
                   p= p * attr[[1]])]
    data.table::setnames(d_temp, old= c("V_NEW", "p"), c(attr_name, prob_name))
    return(d_temp)
  } else {
    d_temp <- dat[, which(names(dat) != prob_name)]  
    d_temp[attr_name] <- attr[[2]]
    d_temp[prob_name] <- p * attr[[1]]
    return(d_temp)
  }
}
