check.for.NA.in.trait <- function(trait=NULL)
{
     ## internal function for AM 
     ## to return the positions of NA in a trait
     ## ordered for largest to smallest (this is important for ReshapeM_rcpp code

       ## check for NA's in trait
        indxNA <- which(is.na(trait))
        if(length(indxNA)==0){
          indxNA <- NULL
        } else {
          ## place in reverse order
          indxNA <- sort(indxNA, decreasing = TRUE)
message(cat("\n\n The phenotypic data in rows ", indxNA, " have missing trait data but their other data are still included.  "))
          if(any(is.na(indxNA))){
            message("Error:  (internal).  indxNA contains NA values. ")
            message(" AM has terminated with errors. ")
            return(NULL)
          }
        }

      return(indxNA)
   }


