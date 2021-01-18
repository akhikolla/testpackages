## \code{$set_data()} reset the data field. ##
lslx$set("public",
         "set_data",
         function(data,
                  sample_cov,
                  sample_mean,
                  sample_size,
                  sample_moment_acov) {
           private$data <- 
             lslxData$new(data = data,
                          sample_cov = sample_cov,
                          sample_mean = sample_mean,
                          sample_size = sample_size,
                          sample_moment_acov = sample_moment_acov,
                          numeric_variable = private$model$numeric_variable,
                          ordered_variable = private$model$ordered_variable,
                          weight_variable = private$model$weight_variable,
                          auxiliary_variable = private$model$auxiliary_variable,
                          group_variable = private$model$group_variable,
                          reference_group = private$model$reference_group,
                          level_group = private$model$level_group,
                          nlevel_ordered = private$model$nlevel_ordered,
                          name_response = private$model$name_response)
           private$fitting <- NULL
         })
