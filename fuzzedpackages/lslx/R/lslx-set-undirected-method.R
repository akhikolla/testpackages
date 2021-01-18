## \code{$free_undirected()} / \code{$penalize_undirected()} / \code{$fix_undirected()} sets all the covariances among \code{both} as FREE / PENALIZED / FIXED. ##
lslx$set("private",
         "set_undirected",
         function(both,
                  group,
                  penalty,
                  set,
                  action,
                  verbose = TRUE) {
           if (length(both) == 1) {
             stop("Argument 'both' must contain more than one varaible.")
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
           name <- paste0(expand.grid(combn(both, 2, function(x)
             paste0(x[1], "<->", x[2])), group)[, 1],
             "/",
             expand.grid(combn(both, 2, function(x)
               paste0(x[1], "<->", x[2])), group)[, 2])
           private$set_coefficient(name = name,
                                   penalty = penalty,
                                   set = set,
                                   action = action,
                                   verbose = verbose)
         })

lslx$set("public",
         "free_undirected",
         function(both,
                  group,
                  verbose = TRUE) {
           private$set_undirected(
             both = both,
             group = group,
             action = "free",
             verbose = verbose
           )
         })

lslx$set("public",
         "fix_undirected",
         function(both,
                  group,
                  verbose = TRUE) {
           private$set_undirected(
             both = both,
             group = group,
             action = "fix",
             verbose = verbose
           )
         })

lslx$set("public",
         "penalize_undirected",
         function(both,
                  group,
                  penalty,
                  set,
                  verbose = TRUE) {
           private$set_undirected(
             both = both,
             group = group,
             penalty = penalty,
             set = set,
             action = "penalize",
             verbose = verbose
           )
         })
