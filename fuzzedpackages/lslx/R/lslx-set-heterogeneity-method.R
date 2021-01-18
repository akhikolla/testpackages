## \code{$free_heterogeneity()} / \code{$penalize_heterogeneity()} / \code{$fix_heterogeneity()} sets every target coefficient as FREE / PENALIZED / FIXED. ##
lslx$set("private",
         "set_heterogeneity",
         function(block,
                  group,
                  penalty,
                  set,
                  action,
                  hold_fixed,
                  verbose = TRUE) {
           if (!all(block %in% unique(private$model$specification$block))) {
             stop("Argument 'block' contains invalid relation. Please check.")
           }
           if (is.null(private$model$reference_group)) {
             stop(
               "Reference group is not specified.",
               "\n The 'set heterogeneity method' can be used only if reference group is specified"
             )
           } else {
             if (group == private$model$reference_group) {
               stop(
                 "Argument 'group' is set as the reference group, please select another group.",
                 "\n  Group name specified in 'group' is ",
                 do.call(paste, as.list(group)),
                 "."
               )
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
             }
           }
           if (hold_fixed) {
             relation <-
               private$model$specification[(private$model$specification$block %in% block) &
                                             (private$model$specification$type !=
                                                "fixed") &
                                             private$model$specification$reference, "relation"]
           } else {
             relation <-
               private$model$specification[(private$model$specification$block %in% block) &
                                             private$model$specification$reference, "relation"]
           }
           if (length(relation) == 0) {
             stop(
               "No valid relation ",
               do.call(paste, as.list(block[!(block %in% unique(private$model$specification$block))])),
               "under group ",
               do.call(paste, as.list(group)),
               "is found. Please check the settings."
             )
           } else {
             name <- paste0(relation, "/", group)
             private$set_coefficient(name = name,
                                     penalty = penalty,
                                     set = set,
                                     action = action,
                                     verbose = verbose)
           }
           
         })

lslx$set("public",
         "free_heterogeneity",
         function(block,
                  group,
                  hold_fixed = TRUE,
                  verbose = TRUE) {
           private$set_heterogeneity(
             block = block,
             group = group,
             action = "free",
             hold_fixed = hold_fixed, 
             verbose = verbose
           )
         })

lslx$set("public",
         "fix_heterogeneity",
         function(block,
                  group,
                  hold_fixed = TRUE,
                  verbose = TRUE) {
           private$set_heterogeneity(
             block = block,
             group = group,
             action = "fix",
             hold_fixed = hold_fixed, 
             verbose = verbose
           )
         })

lslx$set("public",
         "penalize_heterogeneity",
         function(block,
                  group,
                  penalty,
                  set,
                  hold_fixed = TRUE,
                  verbose = TRUE) {
           private$set_heterogeneity(
             block = block,
             group = group,
             penalty = penalty,
             set = set,
             action = "penalize",
             hold_fixed = hold_fixed, 
             verbose = verbose
           )
         })
