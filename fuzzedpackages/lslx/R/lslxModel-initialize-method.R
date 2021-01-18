## \code{$new()} initializes a new \code{lslxModel} object. ##
lslxModel$set("public",
              "initialize",
              function(model,
                       numeric_variable,
                       ordered_variable,
                       weight_variable,
                       auxiliary_variable,
                       group_variable,
                       reference_group,
                       level_group,
                       nlevel_ordered) {
                private$initialize_model(model = model)
                private$initialize_specification()
                private$initialize_tag(
                  numeric_variable = numeric_variable,
                  ordered_variable = ordered_variable,
                  weight_variable = weight_variable,
                  auxiliary_variable = auxiliary_variable,
                  group_variable = group_variable,
                  reference_group = reference_group,
                  level_group = level_group,
                  nlevel_ordered = nlevel_ordered
                )
                private$organize_specification()
                private$expand_specification()
              })

## \code{$initialize_model()} initializes a cleaned model. ##
lslxModel$set("private",
              "initialize_model",
              function(model) {
                model <-
                  gsub(pattern = "[[:blank:]]",
                       replacement = "",
                       x = model)
                model <-
                  gsub(pattern = "\\$|\\?|\\\\|\\^|%|&|#|\\[|\\]|\\{|\\}",
                       replacement = "",
                       x = model)
                model <-
                  gsub(pattern = ";",
                       replacement = "\n",
                       x = model)
                model <-
                  gsub(pattern = "\n{2,}",
                       replacement = "\n",
                       x = model)
                model <-
                  unlist(x = strsplit(x = model,
                                      split = "\n"),
                         use.names = FALSE)
                model <-
                  gsub(pattern = "=~",
                       replacement = ":=>",
                       x = model)
                model <-
                  gsub(pattern = "~~",
                       replacement = "<=>",
                       x = model)
                model <-
                  gsub(pattern = "~\\*~",
                       replacement = "\\*=\\*",
                       x = model)
                model <-
                  ifelse(
                    !grepl(pattern = "\\|=|=\\||\\|~|~\\|",
                           x = model),
                    gsub(
                      pattern = "\\|",
                      replacement = "|=",
                      x = model
                    ),
                    model
                  )
                model <-
                  ifelse(
                    !grepl(pattern = "<~|<~:|~>|:~>|<~>|\\|~|~\\|",
                           x = model),
                    gsub(
                      pattern = "~",
                      replacement = "<=",
                      x = model
                    ),
                    model
                  )
                self$model <- model
              })


## \code{$initialize_specification()} initializes a specification table. ##
lslxModel$set("private",
              "initialize_specification",
              function() {
                operator <- c("|=",
                              "=|",
                              "|~",
                              "~|",
                              "*=*",
                              "*~*",
                              "<=:",
                              "<~:",
                              ":=>",
                              ":~>",
                              "<=",
                              "<~",
                              "=>",
                              "~>",
                              "<=>",
                              "<~>")
                self$specification <-
                  do.call(what = rbind,
                          args = lapply(
                            X = self$model,
                            FUN = function(model_i) {
                              operator_i <- operator[sapply(
                                X = c(
                                  "\\|=",
                                  "=\\|",
                                  "\\|~",
                                  "~\\|",
                                  "\\*=\\*",
                                  "\\*~\\*",
                                  "<=:",
                                  "<~:",
                                  ":=>",
                                  ":~>",
                                  "<=[^:>]",
                                  "<~[^:>]",
                                  "[^:<]=>",
                                  "[^:<]~>",
                                  "<=>",
                                  "<~>"
                                ),
                                FUN = function(pattern) {
                                  grepl(pattern, model_i)
                                }
                              )]
                              if (length(operator_i) > 0) {
                                model_i_split <-
                                  strsplit(x = model_i, split = operator_i, fixed = TRUE)[[1]]
                                if (operator_i %in%
                                    c(":=>", ":~>", "=>", "~>", "=|", "~|")) {
                                  if (operator_i %in%
                                      c(":=>", ":~>", "=>", "~>")) {
                                    operator_i <-
                                      paste0(rev(gsub(
                                        pattern = ">",
                                        replacement = "<",
                                        x = substring(operator_i,
                                                      1:nchar(operator_i),
                                                      1:nchar(operator_i))
                                      )),
                                      collapse = "")
                                  } else {
                                    operator_i <-
                                      paste0(rev(substring(
                                        operator_i,
                                        1:nchar(operator_i),
                                        1:nchar(operator_i)
                                      )),
                                      collapse = "")
                                  }
                                  model_i <-
                                    c(left = model_i_split[2],
                                      operator = operator_i,
                                      right = model_i_split[1])
                                } else {
                                  model_i <-
                                    c(left = model_i_split[1],
                                      operator = operator_i,
                                      right = model_i_split[2])
                                }
                                left_i <-
                                  unlist(strsplit(x = model_i[["left"]],
                                                  split = "\\+"),
                                         use.names = FALSE)
                                right_i <-
                                  unlist(strsplit(x = model_i[["right"]],
                                                  split = "\\+"),
                                         use.names = FALSE)
                                model_i <-
                                  expand.grid(
                                    relation = NA_character_,
                                    left = left_i,
                                    right = right_i,
                                    operator =  operator_i,
                                    KEEP.OUT.ATTRS = FALSE,
                                    stringsAsFactors = FALSE
                                  )
                                left_i_split <-
                                  strsplit(model_i$left, split = "\\*")
                                right_i_split <-
                                  strsplit(model_i$right, split = "\\*")
                                model_i$left <-
                                  sapply(
                                    left_i_split,
                                    FUN = function(i)
                                      getElement(i, length(i))
                                  )
                                model_i$right <-
                                  sapply(
                                    right_i_split,
                                    FUN = function(i)
                                      getElement(i, length(i))
                                  )
                                if (any(model_i$right[model_i$operator %in%
                                                      c("<=>", "<~>")] == "1") |
                                    any(model_i$left[model_i$operator %in%
                                                     c("<=>", "<~>")] == "1")) {
                                  stop(
                                    "Intercept term '1' cannot present at any side of expression for covariance specification."
                                  )
                                }
                                if (any(model_i$left[!(model_i$operator %in%
                                                       c("<=>", "<~>"))] == "1")) {
                                  stop("Intercept term '1' cannot present at the arrow side of expression.")
                                }
                                if (any(model_i$right[model_i$operator %in%
                                                      c("|=", "|~")] == "1") |
                                    any(model_i$left[model_i$operator %in%
                                                     c("|=", "|~")] == "1")) {
                                  stop(
                                    "Intercept term '1' cannot present at any side of expression for threshold specification."
                                  )
                                }
                                if (any(model_i$right[model_i$operator %in%
                                                      c("*=*", "*~*")] == "1") |
                                    any(model_i$left[model_i$operator %in%
                                                     c("*=*", "*~*")] == "1")) {
                                  stop(
                                    "Intercept term '1' cannot present at any side of expression for scale specification."
                                  )
                                }

                                
                                model_i$left_prefix <-
                                  sapply(
                                    left_i_split,
                                    FUN = function(i) {
                                      ifelse(length(i) == 1L, NA_character_, i[1])
                                    }
                                  )
                                model_i$right_prefix <-
                                  sapply(
                                    right_i_split,
                                    FUN = function(i) {
                                      ifelse(length(i) == 1L, NA_character_, i[1])
                                    }
                                  )
                                if (any(!(is.na(model_i$left_prefix) |
                                          is.na(model_i$right_prefix)))) {
                                  stop("Prefix before '*' cannot simultaneously present at both side of expression.")
                                } else if (any(!is.na(model_i$left_prefix))) {
                                  model_i$prefix <-
                                    model_i$left_prefix
                                } else if (any(!is.na(model_i$right_prefix))) {
                                  model_i$prefix <-
                                    model_i$right_prefix
                                } else {
                                  model_i$prefix <- NA_character_
                                }
                                model_i$left_prefix <- NULL
                                model_i$right_prefix <- NULL
                                model_i$relation <-
                                  paste0(model_i$left,
                                         ifelse(
                                           model_i$operator %in%
                                             c("<=:", "<~:", "<=", "<~"),
                                           "<-",
                                           ifelse(model_i$operator %in%
                                                    c("<=>", "<~>"),
                                                  "<->",
                                                  ifelse(model_i$operator %in% 
                                                           c("|=", "|~"),
                                                         "|",
                                                         "**"))),
                                         model_i$right)
                                model_i$operator <-
                                  ifelse(
                                    model_i$right == "1",
                                    ifelse(
                                      model_i$operator == "<=:",
                                      "<=",
                                      ifelse(model_i$operator == "<~:",
                                             "<~",
                                             model_i$operator)
                                    ),
                                    model_i$operator
                                  )
                              } else {
                                model_i <- data.frame()
                              }
                              return(model_i)
                            }
                          ))
                if (any(grepl(
                  pattern = "[[:digit:]]",
                  x = substr(
                    x = setdiff(x = unique(
                      c(self$specification$left,
                        self$specification$right)
                    ),
                    y = c("1")),
                    start = 1,
                    stop = 1
                  )
                ))) {
                  stop(
                    "Names of variable(s) or factor(s) cannot start with numbers.",
                    "\n  Please check the specified 'model'."
                  )
                }
              })


## \code{$initialize_tag()} initializes tags for variables. ##
lslxModel$set("private",
              "initialize_tag",
              function(numeric_variable,
                       ordered_variable,
                       weight_variable,
                       auxiliary_variable,
                       group_variable,
                       reference_group,
                       level_group,
                       nlevel_ordered) {
                self$numeric_variable <- numeric_variable
                self$ordered_variable <- ordered_variable
                self$weight_variable <- weight_variable
                self$auxiliary_variable <- auxiliary_variable
                self$group_variable <- group_variable
                self$reference_group <- reference_group
                self$level_group <- level_group
                self$name_factor <-
                  unique(self$specification[self$specification$operator %in%
                                              c("<=:", "<~:"),
                                            "right"])
                self$name_response <-
                  setdiff(x = unique(unlist(self$specification[!(self$specification$operator %in%
                                                                   c("|=", "|~")),
                                                               c("left", "right")])),
                          y = c(self$name_factor, "1"))
                if (!all(self$name_response %in% 
                         union(x = self$numeric_variable,
                               y = self$ordered_variable))) {
                  stop(
                    "Some response variable in 'model' cannot be found in 'data' or 'sample_cov'.",
                    "\n  Response variables specified by 'model' are ",
                    do.call(paste, as.list(self$name_response)),
                    ".",
                    "\n  Column names of 'data' or 'sample_cov' are ",
                    do.call(paste, as.list(union(x = self$numeric_variable,
                                                 y = self$ordered_variable))),
                    "."
                  )
                } 
                self$numeric_variable <- 
                  intersect(x = self$name_response,
                            y = self$numeric_variable)
                self$ordered_variable <- 
                  intersect(x = self$name_response,
                            y = self$ordered_variable)
                if (length(self$ordered_variable) > 0) {
                  self$nlevel_ordered <- 
                    nlevel_ordered[self$ordered_variable] 
                  self$name_threshold <- paste0("t", 1:(max(nlevel_ordered) - 1))
                } else {
                  self$nlevel_ordered <- numeric(0)
                  self$name_threshold <- character(0)
                }
                self$auxiliary_variable <-
                  setdiff(x = self$auxiliary_variable,
                          y = self$name_response)
                self$name_eta <-
                  c(self$name_response, self$name_factor)
                self$name_endogenous <-
                  unique(self$specification[(self$specification$operator %in%
                                               c("<=:", "<~:", "<=", "<~")) &
                                              (self$specification$right != "1"),
                                            "left"])
                self$name_exogenous <-
                  unique(setdiff(x = self$name_eta,
                                 y = self$name_endogenous))
              })


## \code{$organize_specification()} organizes the specification table. ##
lslxModel$set("private",
              "organize_specification",
              function() {
                self$specification <-
                  do.call(what = rbind,
                          args = lapply(
                            split(self$specification, 1:nrow(self$specification)),
                            FUN = function(specification_i) {
                              if ((specification_i[["operator"]] %in% c("<=>", "<~>")) &
                                  (
                                    match(specification_i[["left"]], self$name_eta) <
                                    match(specification_i[["right"]], self$name_eta)
                                  )) {
                                specification_i <-
                                  data.frame(
                                    relation = paste0(specification_i[["right"]],
                                                      "<->",
                                                      specification_i[["left"]]),
                                    left = specification_i[["right"]],
                                    right = specification_i[["left"]],
                                    operator = specification_i[["operator"]],
                                    prefix = specification_i[["prefix"]]
                                  )
                              }
                              return(specification_i)
                            }
                          ))
                self$specification <-
                  self$specification[!duplicated(self$specification$relation,
                                                 fromLast = TRUE),]
                if (length(self$ordered_variable) > 0) {
                  relation_gamma <-
                    unlist(mapply(
                      FUN = function(ordered_variable_i, 
                                     nlevel_ordered_i) {
                        paste0(ordered_variable_i,
                               "|",
                               paste0("t", 1:(nlevel_ordered_i - 1)))
                      },
                      self$ordered_variable,
                      self$nlevel_ordered),
                      use.names = FALSE)
                  if (!all(self$specification[(self$specification$operator %in%
                                               c("|=", "|~")),
                                              "relation"] %in% (relation_gamma))) {
                    stop("Some specification for thresholds are not recognized.")
                  }
                  relation_gamma <-
                    setdiff(x = relation_gamma, y = self$specification$relation)
                } else {
                  relation_gamma <- character(0)
                }
                if (length(relation_gamma) > 0) {
                  specification_gamma <-
                    data.frame(
                      relation = relation_gamma,
                      left = substr(
                        relation_gamma,
                        start = 1,
                        stop = regexpr("\\|", relation_gamma) - 1
                      ),
                      right = substr(
                        relation_gamma,
                        start = regexpr("\\|", relation_gamma) + 1,
                        stop = nchar(relation_gamma)
                      ),
                      operator = "|=",
                      prefix = NA_character_,
                      stringsAsFactors = FALSE
                    )
                } else {
                  specification_gamma = data.frame()
                }
                if (length(self$ordered_variable) > 0) {
                  relation_psi <- paste0(self$ordered_variable, 
                                         "**", 
                                         self$ordered_variable)
                  
                  if (!all(self$specification[(self$specification$operator %in%
                                               c("*=*", "*~*")),
                                              "relation"] %in% (relation_psi))) {
                    stop("Some specification for scale are not recognized.")
                  }
                  relation_psi <-
                    setdiff(x = relation_psi, y = self$specification$relation)
                } else {
                  relation_psi <- character(0)
                }
                if (length(relation_psi) > 0) {
                  specification_psi <-
                    data.frame(
                      relation = relation_psi,
                      left = substr(
                        relation_psi,
                        start = 1,
                        stop = regexpr("\\*\\*", relation_psi) - 1
                      ),
                      right = substr(
                        relation_psi,
                        start = regexpr("\\*\\*", relation_psi) + 2,
                        stop = nchar(relation_psi)
                      ),
                      operator = "*=*",
                      prefix = 1,
                      stringsAsFactors = FALSE
                    )
                } else {
                  specification_psi = data.frame()
                }
                
                if (any(self$specification$right == "1")) {
                  if (length(intersect(x = self$numeric_variable,
                                       y = self$name_exogenous)) > 0) {
                    relation_alpha <-
                      setdiff(
                        x = paste(
                          intersect(
                            x = self$numeric_variable,
                            y = self$name_exogenous
                          ),
                          "1",
                          sep = "<-"
                        ),
                        y = self$specification$relation
                      )
                  } else {
                    relation_alpha <- character()
                  }
                } else {
                  if (length(self$numeric_variable) > 0) {
                    relation_alpha <-
                      setdiff(
                        x = paste(self$numeric_variable,
                                  "1",
                                  sep = "<-"),
                        y = self$specification$relation
                      )
                  } else {
                    relation_alpha <- character()
                  }
                }
                if (length(relation_alpha) > 0) {
                  specification_alpha <- data.frame(
                    relation = relation_alpha,
                    left = substr(
                      relation_alpha,
                      start = 1,
                      stop = regexpr("<-", relation_alpha) - 1
                    ),
                    right = "1",
                    operator = "<=",
                    prefix = NA_character_,
                    stringsAsFactors = FALSE
                  )
                } else {
                  specification_alpha = data.frame()
                }
                if (length(self$name_exogenous) > 1) {
                  relation_phi <-
                    setdiff(x = c(
                      paste0(self$name_eta,
                             "<->",
                             self$name_eta),
                      paste0(apply(
                        combn(rev(self$name_exogenous), 2),
                        2,
                        FUN = function(i) {
                          paste(i[1], i[2], sep = "<->")
                        }
                      ))
                    ),
                    y = self$specification$relation)
                } else {
                  relation_phi <-
                    setdiff(
                      x = paste0(self$name_eta,
                                 "<->",
                                 self$name_eta),
                      y = self$specification$relation
                    )
                }
                if (length(relation_phi) > 0) {
                  specification_phi <-
                    data.frame(
                      relation = relation_phi,
                      left = substr(
                        relation_phi,
                        start = 1,
                        stop = regexpr("<->", relation_phi) - 1
                      ),
                      right = substr(
                        relation_phi,
                        start = regexpr("<->", relation_phi) + 3,
                        stop = nchar(relation_phi)
                      ),
                      operator = "<=>",
                      prefix = NA_character_,
                      stringsAsFactors = FALSE
                    )
                } else {
                  specification_phi = data.frame()
                }
                self$specification <-
                  rbind(
                    self$specification,
                    specification_gamma,
                    specification_alpha,
                    specification_phi,
                    specification_psi,
                    make.row.names = FALSE,
                    stringsAsFactors = FALSE
                  )
              })



## \code{$expand_specification()} expands the specification table. ##
lslxModel$set("private",
              "expand_specification",
              function() {
                prefix_split <-
                  lapply(X = self$specification$prefix,
                         FUN = function(prefix_i) {
                           if (!is.na(prefix_i)) {
                             if ((substr(prefix_i, 
                                         start = 1, 
                                         stop = 2) == "c(") & 
                                 (substr(prefix_i, 
                                         start = nchar(prefix_i), 
                                         stop = nchar(prefix_i)) == ")")) {
                               prefix_i <- 
                                 substr(prefix_i, start = 3, stop = (nchar(prefix_i) - 1))
                               prefix_i_split <- character(0)
                               while (nchar(prefix_i) > 0) {
                                 idx_left_bracket <- regexpr("\\(", prefix_i)[1]
                                 idx_right_bracket <- regexpr("\\)", prefix_i)[1]
                                 idx_comma <- regexpr(",", prefix_i)[1]
                                 if (idx_comma == -1) {
                                   prefix_i_split <- c(prefix_i_split, prefix_i)
                                   prefix_i <- ""
                                 } else {
                                   if ((idx_left_bracket == -1) &
                                       (idx_right_bracket == -1)) {
                                     prefix_i_split <- c(prefix_i_split, strsplit(prefix_i, ",")[[1]])
                                     prefix_i <- ""
                                   } else {
                                     if ((idx_left_bracket < idx_comma) & 
                                         (idx_right_bracket > idx_comma)) {
                                       prefix_i_split <- c(prefix_i_split, 
                                                           substr(prefix_i, 
                                                                  start = 1, 
                                                                  stop = idx_right_bracket))
                                       if (idx_right_bracket < nchar(prefix_i)) {
                                         prefix_i <-
                                           substr(prefix_i, 
                                                  start = (idx_right_bracket + 2), 
                                                  stop = nchar(prefix_i))
                                       } else {
                                         prefix_i <- ""
                                       }
                                     } else {
                                       prefix_i_split <- c(prefix_i_split, 
                                                           substr(prefix_i, 
                                                                  start = 1, 
                                                                  stop = (idx_comma - 1)))
                                       prefix_i <-
                                         substr(prefix_i, 
                                                start = (idx_comma + 1), 
                                                stop = nchar(prefix_i))
                                     }
                                   }
                                 }
                               }
                             } else {
                               prefix_i_split <- prefix_i
                             }
                           } else {
                             prefix_i_split <- prefix_i
                           }
                           prefix_i_split <- 
                             ifelse(gsub(pattern = "\\(.*$", replacement = "", x = prefix_i_split) %in% 
                                      c("free", "fix", "pen", "start", "lab"),
                                    prefix_i_split,
                                    ifelse(is.na(prefix_i_split),
                                           prefix_i_split,
                                           ifelse(prefix_i_split == "NA",
                                                  NA,
                                                  ifelse(suppressWarnings(!is.na(as.numeric(prefix_i_split))),
                                                         paste0("fix", "(",  prefix_i_split, ")"),
                                                         paste0("lab", "(",  prefix_i_split, ")")))))
                           return(prefix_i_split)
                         })
                if (anyNA(self$reference_group)) {
                  if (any(sapply(X = prefix_split,
                                 FUN = function(prefix_split_i) {
                                   return(length(prefix_split_i) > 1)
                                 }))) {
                    stop("When 'reference_group = NA', vectorized prefix cannot be used.")
                  }
                  if (length(self$ordered_variable) > 0) {
                    stop("When 'reference_group = NA', response variable cannot be ordered.")
                  }
                  self$specification <- 
                    do.call(what = rbind,
                            args = lapply(
                              X = c("<NA>", self$level_group),
                              FUN = function(level_group_i) {
                                specification_i <- self$specification
                                specification_i$group <- level_group_i
                                specification_i$reference <-
                                  ifelse(level_group_i == "<NA>", TRUE, FALSE)
                                specification_i$operator <-
                                  ifelse(specification_i$group == "<NA>", 
                                         specification_i$operator, 
                                         gsub(pattern = "=",
                                              replacement = "~",
                                              x = specification_i$operator))
                                specification_i$prefix <-
                                  sapply(
                                    X = prefix_split,
                                    FUN = function(prefix_split_j) {
                                      prefix_split_j_verb <- 
                                        gsub(pattern = "\\(.*$", 
                                             replacement = "", 
                                             x = prefix_split_j)
                                      if (prefix_split_j_verb %in% "pen") {
                                        prefix_split_j <- 
                                          eval(parse(text = prefix_split_j))
                                      } 
                                      if (level_group_i != "<NA>") {
                                        if (prefix_split_j_verb %in% "fix") {
                                          prefix_split_j <- "fix(0)"
                                        } else {
                                          prefix_split_j <- NA
                                        }
                                      } 
                                      return(prefix_split_j)
                                    }
                                  )
                                rownames(specification_i) <-
                                  paste0(specification_i$relation,
                                         "/",
                                         specification_i$group)
                                return(specification_i)
                              }
                            ))
                } else {
                  if (any(sapply(X = prefix_split,
                                 FUN = function(prefix_split_i) {
                                   return((length(prefix_split_i) > 1) & 
                                          (length(prefix_split_i) != length(self$level_group)))
                                 }))) {
                    stop("The length of prefix vector should be 1 or equal to to the number of groups.")
                  }
                  self$specification <- 
                    do.call(what = rbind,
                            args = lapply(
                              X = self$level_group,
                              FUN = function(level_group_i) {
                                specification_i <- self$specification
                                specification_i$group <- level_group_i
                                specification_i$reference <-
                                  ifelse(
                                    is.null(self$reference_group),
                                    FALSE,
                                    ifelse(level_group_i == self$reference_group,
                                           TRUE,
                                           FALSE)
                                  )
                                specification_i$prefix <-
                                  sapply(
                                    X = prefix_split,
                                    FUN = function(prefix_split_j) {
                                      if (length(prefix_split_j) == 1) {
                                        prefix_split_j_verb <- 
                                          gsub(pattern = "\\(.*$", 
                                               replacement = "", 
                                               x = prefix_split_j)
                                        if (prefix_split_j_verb %in% "pen") {
                                          prefix_split_j <- 
                                            eval(parse(text = prefix_split_j))
                                        } 
                                        if (!is.null(self$reference_group)) {
                                          if (level_group_i != self$reference_group) {
                                            if (prefix_split_j_verb %in% 
                                                c("free", "fix", "start")) {
                                              prefix_split_j <- 
                                                paste0(prefix_split_j_verb, "(", 0, ")")
                                            } else if (prefix_split_j_verb %in% "pen") {
                                              prefix_split_j <- 
                                                paste0(substr(prefix_split_j, 
                                                              start = 1, 
                                                              stop = regexpr("=", prefix_split_j)[1]),
                                                       0,
                                                       substr(prefix_split_j, 
                                                              start = regexpr(",", prefix_split_j)[1], 
                                                              stop = nchar(prefix_split_j)))
                                            } else if (prefix_split_j_verb %in% c("lab")) {
                                              prefix_split_j <- "fix(0)"
                                            } else {
                                              prefix_split_j <- NA
                                            }
                                          } 
                                        }
                                      } else {
                                        prefix_split_j <-
                                          prefix_split_j[level_group_i == self$level_group]
                                        prefix_split_j_verb <- 
                                          gsub(pattern = "\\(.*$", 
                                               replacement = "", 
                                               x = prefix_split_j)
                                        if (prefix_split_j_verb %in% "pen") {
                                          prefix_split_j <- 
                                            eval(parse(text = prefix_split_j))
                                        } 
                                      }
                                      return(prefix_split_j)
                                    }
                                  )
                                rownames(specification_i) <-
                                  paste0(specification_i$relation,
                                         "/",
                                         specification_i$group)
                                return(specification_i)
                              }
                            ))
                }
                self$specification$matrix <-
                  factor(
                    x = ifelse(
                      self$specification$operator %in% c("|=", "|~"),
                      "gamma",
                      ifelse(self$specification$operator %in% c("*=*", "*~*"),
                             "psi",
                             ifelse(
                               self$specification$operator %in% c("<=>", "<~>"),
                               "phi",
                               ifelse(self$specification$right == "1",
                                      "alpha",
                                      "beta")
                             ))),
                    levels = c("gamma", "alpha", "beta", "phi", "psi")
                  )
                self$specification$block <-
                  with(self$specification, {
                    block_left <-
                      ifelse(left %in% self$name_response,
                             "y",
                             "f")
                    block_right <-
                      ifelse(matrix %in% "gamma",
                             "t",
                             ifelse(
                               right %in% self$name_response,
                               "y",
                               ifelse(
                                 right %in% self$name_factor,
                                 "f",
                                 "1"
                               )
                             ))
                    block_middle <-
                      ifelse(matrix %in% "gamma",
                             "|",
                             ifelse(matrix %in% "psi",
                                    "**",
                                    ifelse(matrix %in% c("alpha", "beta"),
                                           "<-",
                                           "<->")))
                    paste0(block_left, block_middle, block_right)
                  })
                
                self$specification$type <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    type <- 
                      ifelse(prefix_verb %in% "free",
                             "free",
                             ifelse(prefix_verb %in% "fix",
                                    "fixed",
                                    ifelse(prefix_verb %in% "pen",
                                           "pen",
                                           ifelse(
                                             operator %in% c("|=", "*=*", "<=:", "<=", "<=>"),
                                             "free",
                                             "pen"))))
                    return(type)
                  })
                self$specification$type <-
                  ifelse(self$specification$relation %in% 
                           c(paste0(self$ordered_variable, 
                                    "<->", 
                                    self$ordered_variable)),
                         "fixed",
                         self$specification$type)
                self$specification$start <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    prefix_value <- 
                      ifelse(prefix_verb %in% c("free", "fix", "start"),
                             substr(x = prefix, 
                                    start = (regexpr("\\(", prefix) + 1), 
                                    stop = (nchar(prefix) - 1)),
                             ifelse(prefix_verb %in% "pen",
                                    sapply(X = strsplit(prefix, split = "=|,|\\(|\\)"), 
                                           FUN = function(x) {
                                             x[3]},
                                           simplify = TRUE, USE.NAMES = FALSE),
                                    NA_real_))
                    start <- 
                      ifelse(prefix_verb %in% c("free", "fix", "pen", "start"),
                             suppressWarnings(as.numeric(prefix_value)),
                             NA_real_)
                    return(start)
                  })
                self$specification$start <-
                  ifelse(self$specification$relation %in% 
                           paste0(self$ordered_variable, 
                                  "<->", 
                                  self$ordered_variable),
                         NA_real_,
                         self$specification$start)
                self$specification$label <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    label <- 
                      ifelse(prefix_verb %in% "lab",
                             substr(x = prefix, 
                                    start = (regexpr("\\(", prefix) + 1), 
                                    stop = (nchar(prefix) - 1)),
                             NA_character_)
                    return(label)
                  })
                self$specification$penalty <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    penalty <- 
                      ifelse(prefix_verb %in% "pen",
                             sapply(X = strsplit(prefix, split = "=|,|\\(|\\)"), 
                                    FUN = function(x) {
                                      x[5]},
                                    simplify = TRUE, USE.NAMES = FALSE),
                             ifelse(type == "pen",
                                    "default",
                                    "none"))
                    return(penalty)
                  })
                self$specification$set <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    set <- 
                      ifelse(prefix_verb %in% "pen",
                             sapply(X = strsplit(prefix, split = "=|,|\\(|\\)"), 
                                    FUN = function(x) {
                                      as.numeric(x[7])},
                                    simplify = TRUE, USE.NAMES = FALSE),
                             ifelse(type == "pen", 1, 0))
                    return(set)
                  })
                self$specification$weight <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    weight <- 
                      ifelse(prefix_verb %in% "pen",
                             sapply(X = strsplit(prefix, split = "=|,|\\(|\\)"), 
                                    FUN = function(x) {
                                      x[9]},
                                    simplify = TRUE, USE.NAMES = FALSE),
                             ifelse(type == "pen", 1, 0))
                    return(as.numeric(weight))
                  })
                self$specification$operator <- NULL
                self$specification$prefix <- NULL
                if (length(self$ordered_variable) > 0) {
                  specification_gamma <- self$specification[self$specification$matrix == "gamma", ]
                  specification_non_gamma <- self$specification[self$specification$matrix != "gamma", ]
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
                      match(specification_non_gamma$right, c("1", self$name_threshold, self$name_eta)),
                      match(specification_non_gamma$left, self$name_eta),
                      method = "radix"
                    ),]
                  self$specification <- rbind(specification_gamma, 
                                              specification_non_gamma)
                  self$specification <-
                    self$specification[order(
                      self$specification$reference,
                      self$specification$group,
                      self$specification$matrix,
                      self$specification$block,
                      method = "radix"
                    ),]
                } else {
                  self$specification <-
                    self$specification[order(
                      self$specification$reference,
                      self$specification$group,
                      self$specification$matrix,
                      self$specification$block,
                      match(self$specification$right, c("1", self$name_threshold, self$name_eta)),
                      match(self$specification$left, self$name_eta),
                      method = "radix"
                    ),]
                }
                self$specification <-
                  self$specification[order(
                    self$specification$reference,
                    decreasing = TRUE,
                    method = "radix"
                  ), ]
              })
