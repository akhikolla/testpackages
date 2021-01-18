## \code{$free_directed()} / \code{$penalize_directed()} / \code{$fix_directed()} sets all the regression coefficients from \code{right} to \code{left} as FREE / PENALIZED / FIXED. ##
lslx$set("private",
         "set_directed",
         function(left,
                  right,
                  group,
                  penalty,
                  set,
                  action,
                  verbose = TRUE) {
           if (missing(left)) {
             stop("Argument 'left' must be given.")
           } else if (missing(right)) {
             stop("Argument 'right' must be given.")
           } else {
           }
           
           if (missing(group)) {
             group <-  private$model$level_group
           } else if (!all(group %in% private$model$level_group)) {
             stop(
               "Argument 'group' contains unknown group name.",
               "\n  Group name(s) currently recognized by 'lslx' is ",
               do.call(paste, as.list(private$model$level_group)),
               ".",
               "\n  Group name specified in 'group' is ",
               do.call(paste, as.list(group)),
               "."
             )
           } else {
           }
           name <- paste0(
             expand.grid(left, "<-", right)[, 1],
             expand.grid(left, "<-", right)[, 2],
             expand.grid(left, "<-", right)[, 3],
             "/",
             group
           )
           private$set_coefficient(name = name,
                                   action = action,
                                   penalty = penalty,
                                   set = set,
                                   verbose = verbose)
         })

lslx$set("public",
         "free_directed",
         function(left,
                  right,
                  group,
                  verbose = TRUE) {
           private$set_directed(
             left = left,
             right = right,
             group = group,
             action = "free",
             verbose = verbose
           )
         })

lslx$set("public",
         "fix_directed",
         function(left,
                  right,
                  group,
                  verbose = TRUE) {
           private$set_directed(
             left = left,
             right = right,
             group = group,
             action = "fix",
             verbose = verbose
           )
         })

lslx$set("public",
         "penalize_directed",
         function(left,
                  right,
                  group,
                  penalty,
                  set,
                  verbose = TRUE) {
           private$set_directed(
             left = left,
             right = right,
             group = group,
             penalty = penalty,
             set = set,
             action = "penalize",
             verbose = verbose
           )
         })
