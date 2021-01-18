
doquiet <- function(dat, num_markers, lab){
     ## a diagnostic function for printing the contents of matrix or vector
     ## used for error checking

     if(dim(dat)[1] == 1 || dim(dat)[2] ==1 )
         dat <- as.numeric(dat)

     if(class(dat)=="matrix"){
          ### if dat is a matrix

         if(num_markers > 0){
           message(" Dimension of ", lab, " is ", dim(dat)[1], " by ", dim(dat)[2], " \n")
           message(" First few rows and ", num_markers, " columns of ", lab, " are: \n")
           if(nrow(dat) > 5 && ncol(dat) > num_markers){
             for(xx in 1:5)
               message(sprintf(" %f ", dat[xx, 1:num_markers]))

           }
           if(nrow(dat) <=5  &&  ncol(dat) > num_markers)
             for(xx in 1:nrow(dat))
               message(sprintf(" %f ", dat[xx, 1:num_markers]))
           if(nrow(dat) > 5  &&  ncol(dat) <=  num_markers)
             for(xx in 1:5)
               message(sprintf(" %f ", dat[xx, 1:ncol(dat)]))
           if(nrow(dat) <= 5  &&  ncol(dat) <=  num_markers)
             for(xx in 1:nrow(dat))
               message(sprintf(" %f ", dat[xx, 1:ncol(dat)]))
           message("\n\n")
         }
     } ## end if class(dat)

     if(class(dat)=="numeric" || class(dat)=="vector"){
       if(num_markers > 0){
          message(" Length of ", lab, " is ", length(dat), "\n")
          message(" The first ", num_markers, " elements of the vector are ", lab, "\n")
          if(length(dat) > num_markers)
             message(sprintf(" %f ", dat[1:num_markers]))
          if(length(dat) <= num_markers)
             message(sprintf(" %f ", dat[1:length(dat)]))
       message("\n\n")
       }
    }

    if(!(class(dat)=="matrix" || class(dat)=="vector" || class(dat)=="numeric"))
      message(" Internal error in function doquiet. dat not matrix or vector or numeric. \n")

}



