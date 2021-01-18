## \code{$free_coefficient()} / \code{$penalize_coefficient()} / \code{$fix_coefficient()} sets the coefficient named \code{name} as FREE / PENALIZED / FIXED with starting value \code{start}. ##
lslx$set("public",
         "free_coefficient",
         function(name,
                  start,
                  verbose = TRUE) {
           private$set_coefficient(
             name = name,
             start = start,
             action = "free",
             verbose = verbose
           )
         })

lslx$set("public",
         "fix_coefficient",
         function(name,
                  start,
                  verbose = TRUE) {
           private$set_coefficient(
             name = name,
             start = start,
             action = "fix",
             verbose = verbose
           )
         })

lslx$set("public",
         "penalize_coefficient",
         function(name,
                  start,
                  penalty,
                  set,
                  weight,
                  verbose = TRUE) {
           private$set_coefficient(
             name = name,
             start = start,
             penalty = penalty,
             set = set,
             weight = weight,
             action = "penalize",
             verbose = verbose
           )
         })

## \code{$set_coefficient_type()} reset the types of coefficients in model specification. ##
lslx$set("public",
         "set_coefficient_type",
         function(name, type) {
           private$set_coefficient(
             name = name,
             type = type,
             action = "type",
             verbose = FALSE
           )
         })

## \code{$set_coefficient_start()} reset the starting values of coefficients in model specification. ##
lslx$set("public",
         "set_coefficient_start",
         function(name, start) {
           private$set_coefficient(
             name = name,
             start = start,
             action = "start",
             verbose = FALSE
           )
         })

lslx$set("private",
         "set_coefficient",
         function(name,
                  start,
                  penalty,
                  set,
                  weight,
                  type,
                  action,
                  verbose = TRUE) {
           if (missing(name)) {
             stop("Argument 'name' must be given.")
           }
           name <-
             gsub(pattern = "[[:blank:]]",
                  replacement = "",
                  x = name)
           name <-
             sapply(
               X = name,
               FUN = function(name_i) {
                 if (!grepl(pattern = "<-[^>]|<->|\\||\\*\\*", x = name_i)) {
                   stop("Some coefficient name contains invalid operator.")
                 }
                 if (!grepl(pattern = "/", x = name_i)) {
                   if (length(private$model$level_group) == 1) {
                     name_i <- paste0(name_i, "/", private$model$level_group)
                   } else {
                     stop("Some coefficient name doesn't contain group name.")
                   }
                 }
                 name_i_split <-
                   strsplit(x = name_i, split = c("<-|/|<->|\\||\\*\\*"))[[1]]
                 if (length(name_i_split) != 3) {
                   stop("The format of some coefficient name is incorrect.")
                 }
                 left_i <- name_i_split[1]
                 right_i <- name_i_split[2]
                 group_i <- name_i_split[3]
                 
                 if (left_i == 1) {
                   stop("Intercept term '1' cannot be presented at the left-hand side.")
                 }
                 
                 if (!(left_i %in% private$model$name_eta)) {
                   stop(
                     "Some specified left response or factor name is unrecognized.",
                     "\n  Response or factor name(s) currently recognized by 'lslx' is ",
                     do.call(paste, as.list(private$model$name_eta)),
                     ".",
                     "\n  The unrecognized variable or factor name is ",
                     left_i,
                     "."
                   )
                 }
                 
                 if (!(right_i %in% c(private$model$name_eta, private$model$name_threshold, 1))) {
                   stop(
                     "Some specified right response or factor name is unrecognized.",
                     "\n  Response or factor name(s) currently recognized by 'lslx' is ",
                     do.call(paste, as.list(private$model$name_eta)),
                     ".",
                     "\n  The unrecognized response or factor name is ",
                     right_i,
                     "."
                   )
                 }
                 
                 if (!(group_i %in% private$model$level_group)) {
                   stop(
                     "Some specified group name is unrecognized.",
                     "\n  Group name(s) currently recognized by 'lslx' is ",
                     do.call(paste, as.list(private$model$level_group)),
                     ".",
                     "\n  The unrecognized group name is ",
                     group_i,
                     "."
                   )
                 }
                 if (grepl(pattern = "<->",
                           x = name_i)) {
                   if (right_i == 1) {
                     stop("Intercept term '1' cannot be presented at the right-hand side of '<->'.")
                   }
                   if ((
                     match(left_i, private$model$name_eta) <
                     match(right_i, private$model$name_eta)
                   )) {
                     name_i <-
                       paste0(right_i, "<->", left_i, "/", group_i)
                   }
                 }
                 return(name_i)
               },
               simplify = TRUE,
               USE.NAMES = FALSE
             )
           
           if (action == "start") {
             if (missing(start)) {
               stop("Argument 'start' must be given if 'action' == 'start'.")
             }
             if (!is.numeric(start)) {
               stop("The argument 'start' must be a numeric vector.")
             }
             if (length(start) != length(name)) {
               stop(
                 "Argument 'start' has ambiguous length.",
                 "\n  The length of 'start' must be the same with the length of 'name' if 'action' == 'start'."
               )
             }
             private$model$specification[name, "start"] <- start
           } else if (action == "type") {
             if (missing(type)) {
               stop("Argument 'type' must be given if 'action' == 'type'.")
             }
             if (!is.character(type)) {
               stop("Argument 'type' must be a character vector.")
             }
             if (length(type) == 1) {
               type <- rep(type, length(name))
             }
             if (length(type) != length(name)) {
               stop(
                 "Argument 'type' has ambiguous length.",
                 "\n  The length of 'type' must be the same with the length of 'name' if 'action' == 'type'."
               )
             }
             private$model$specification[name, "type"] <- type
           } else {
             if (missing(start)) {
               if (action == "free") {
                 if (any(private$model$level_group %in% private$model$reference_group)) {
                   start <- rep(0, length(name))
                 } else {
                   start <- rep(NA_real_, length(name))
                 }
               } else {
                 start <- rep(0, length(name))
               }
             } else {
               if (!is.numeric(start)) {
                 stop("The argument 'start' must be a numeric vector.")
               }
               if (length(start) == 1) {
                 start <- rep(start, length(name))
               }
               if (length(start) != length(name)) {
                 stop(
                   "Argument 'start' has ambiguous length.",
                   "\n  The length of 'start' must be the same with the length of 'name' or just one."
                 )
               }
             }
             if (missing(penalty)) {
               if (action == "penalize") {
                 penalty <- rep("default", length(name))
               } else {
                 penalty <- rep("none", length(name))
               }
             } else {
               if (!is.character(penalty)) {
                 stop("Argument 'penalty' must be a character vector.")
               } else {
                 if (!(penalty %in% c("lasso", "ridge", "mcp", "elastic_net"))) {
                   stop("Elements in argument 'penalty' must be 'lasso', 'ridge', 'mcp', or 'elastic_net'.")      
                 }
               }
               if (length(penalty) == 1) {
                 penalty <- rep(penalty, length(name))
               } else {
                 if (length(penalty) != length(name)) {
                   stop(
                     "Argument 'penalty' has ambiguous length.",
                     "\n  The length of 'penalty' must be the same with the length of 'name' or just one."
                   )
                 }
               }
             }
             
             if (missing(set)) {
               if (action == "penalize") {
                 set <- rep(1, length(name))
               } else {
                 set <- rep(0, length(name))
               }
             } else {
               if (!is.numeric(set)) {
                 stop("Argument 'set' must be a numeric vector.")
               } else {
                 if (any(!(set %in% c(1, 2)))) {
                   stop("Elements in argument 'set' must be 1 or 2.")      
                 }
                 if (length(set) == 1) {
                   set <- rep(set, length(name))
                 } else {
                   if (length(set) != length(name)) {
                     stop(
                       "Argument 'set' has ambiguous length.",
                       "\n  The length of 'set' must be the same with the length of 'name' or just one."
                     )
                   }
                 }
               }
               set <- rep(set, length(name))
             }
             if (missing(weight)) {
               if (action == "penalize") {
                 weight <- rep(1, length(name))
               } else {
                 weight <- rep(0, length(name))
               }
             } else {
               if (!is.numeric(weight)) {
                 stop("Argument 'weight' must be a numeric vector.")
               } 
               if (length(weight) == 1) {
                 weight <- rep(weight, length(name))
               } else {
                 if (length(weight) != length(name)) {
                   stop(
                     "Argument 'weight' has ambiguous length.",
                     "\n  The length of 'weight' must be the same with the length of 'name' or just one."
                   )
                 }
               }
             }
             
             if (action == "free") {
               type <- "free"
             } else if (action == "fix") {
               type <- "fixed"
             } else if (action == "penalize") {
               type <- "pen"
             } else {
             }

             for (i in seq_len(length(name))) {
               if (name[i] %in% rownames(private$model$specification)) {
                 private$model$specification[name[i], "type"] <- type
                 private$model$specification[name[i], "start"] <-
                   start[i]
                 private$model$specification[name[i], "penalty"] <-
                   penalty[i]
                 private$model$specification[name[i], "set"] <-
                   set[i]
                 private$model$specification[name[i], "weight"] <-
                   weight[i]
                 specification_i <-
                   private$model$specification[name[i], , drop = FALSE]
               } else {
                 relation_i <-
                   substr(name[i],
                          start = 1,
                          stop = regexpr("/", name[i]) - 1)
                 name_i_split <-
                   strsplit(x = name[i], split = c("<-|/|<->|\\||\\*\\*"))[[1]]
                 left_i <- name_i_split[1]
                 right_i <- name_i_split[2]
                 group_i <- name_i_split[3]
                 
                 matrix_i <-
                   ifelse(grepl(pattern = "\\*\\*",
                                x = name[i]),
                          "psi",
                          ifelse(grepl(pattern = "\\|",
                                       x = name[i]),
                                 "gamma",
                                 ifelse(grepl(pattern = "<->",
                                              x = name[i]),
                                        "phi",
                                        ifelse(
                                          grepl(pattern = "<-[^>]",
                                                x = name[i]),
                                          ifelse(right_i == "1",
                                                 "alpha",
                                                 "beta")
                                        ))))
                 block_i_left <-
                   ifelse(left_i %in% private$model$name_response,
                          "y",
                          "f")
                 block_i_right <-
                   ifelse(
                     right_i %in% private$model$name_response,
                     "y",
                     ifelse(
                       right_i %in% private$model$name_factor,
                       "f",
                       ifelse(right_i == "1",
                              "1",
                              "t")
                     )
                   )
                 block_i_middle <-
                   ifelse(
                     matrix_i %in% c("alpha", "beta"),
                     "<-",
                     ifelse(matrix_i == "phi",
                            "<->",
                            ifelse("gamma",
                                   "|",
                                   "**"))
                   )
                 block_i <-
                   paste0(block_i_left, block_i_middle, block_i_right)
                 
                 specification_i <-
                   data.frame(
                     relation = relation_i,
                     left = left_i,
                     right = right_i,
                     group = group_i,
                     reference =
                       ifelse(
                         is.null(private$model$reference_group),
                         FALSE,
                         ifelse(group_i == private$model$reference_group,
                                TRUE,
                                FALSE)
                       ),
                     matrix = matrix_i,
                     block = block_i,
                     type = type,
                     start = start[i],
                     label = NA_character_,
                     penalty = penalty[i],
                     set = set[i],
                     weight = weight[i],
                     stringsAsFactors = FALSE
                   )
                 rownames(specification_i) <- name[i]
                 private$model$specification <-
                   rbind(private$model$specification,
                         specification_i)
               }
               if (verbose) {
                 cat(
                   "The relation",
                   specification_i$relation,
                   "under",
                   specification_i$group,
                   "is set as",
                   ifelse(
                     type == "pen",
                     "PENALIZED",
                     ifelse(type == "free",
                            "FREE",
                            "FIXED")
                   ),
                   "with starting value =",
                   paste0(specification_i$start, "."),
                   "\n"
                 )
               }
             }
             
             if (verbose & !is.null(private$model$reference_group)) {
               cat(
                 "NOTE: Because",
                 private$model$reference_group,
                 "is set as reference,",
                 "a relation under other group",
                 "actually represents an increment. \n"
               )
               cat(
                 "NOTE: Please check whether the starting value for the increment represents a difference. \n"
               )
             }
             
             if (length(private$model$ordered_variable) > 0) {
               specification_gamma <- private$model$specification[private$model$specification$matrix == "gamma", ]
               specification_non_gamma <- private$model$specification[private$model$specification$matrix != "gamma", ]
               specification_gamma <- 
                 specification_gamma[order(
                   specification_gamma$reference,
                   specification_gamma$group,
                   specification_gamma$matrix,
                   specification_gamma$block,
                   match(specification_gamma$left, self$name_eta),
                   match(specification_gamma$right, c("1", self$name_threshold, self$name_eta)),
                   method = "radix"
                 ),]
               specification_non_gamma <- 
                 specification_non_gamma[order(
                   specification_non_gamma$reference,
                   specification_non_gamma$group,
                   specification_non_gamma$matrix,
                   specification_non_gamma$block,
                   match(specification_non_gamma$right, c("1", private$model$name_threshold, private$model$name_eta)),
                   match(specification_non_gamma$left, private$model$name_eta),
                   method = "radix"
                 ),]
               private$model$specification <- rbind(specification_gamma, 
                                                    specification_non_gamma)
               private$model$specification <-
                 private$model$specification[order(
                   private$model$specification$reference,
                   private$model$specification$group,
                   private$model$specification$matrix,
                   private$model$specification$block,
                   method = "radix"
                 ),]
             } else {
               private$model$specification <-
                 private$model$specification[order(
                   private$model$specification$reference,
                   private$model$specification$group,
                   private$model$specification$matrix,
                   private$model$specification$block,
                   match(private$model$specification$right, c("1", private$model$name_threshold, private$model$name_eta)),
                   match(private$model$specification$left, private$model$name_eta),
                   method = "radix"
                 ), ]
             }
             private$model$specification <-
               private$model$specification[order(
                 private$model$specification$reference,
                 decreasing = TRUE,
                 method = "radix"
               ), ]
             
             private$model$name_endogenous <-
               unique(private$model$specification$left[private$model$specification$matrix == "beta"])
             private$model$name_exogenous <-
               setdiff(x = private$model$name_eta,
                       y = private$model$name_endogenous)
           }
           private$fitting <- NULL
         })
