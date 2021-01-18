## \code{$get_model()} returns a deep copy of \code{model} member in the current \code{lslx} object. ##
lslx$set("public",
         "get_model",
         function() {
           return(private$model$clone(deep = TRUE))
         })

## \code{$get_data()} returns a deep copy of \code{data} member in the current \code{lslx} object. ##
lslx$set("public",
         "get_data",
         function() {
           return(private$data$clone(deep = TRUE))
         })

## \code{$get_fitting()} returns a deep copy of \code{fitting} member in the current \code{lslx} object. ##
lslx$set("public",
         "get_fitting",
         function() {
           return(private$fitting$clone(deep = TRUE))
         })
