 .print_results <- function(itnum=NULL, selected_loci, map, extBIC)
 ## internal function: used by AM
 {  if(!is.null(itnum)){
       message(" Significant marker-trait association found. \n")
       message(" New results after iteration ", itnum, " are \n")
    }
    message(sprintf("%15s  %10s        %10s     %10s        %10s ",
                 "SNP", "Chrm", "Map Pos",  "Col Number",       "extBIC"))
    message(sprintf("%15s  %10s        %10s     %10s        %10s ",
                 "-----", "------", "---------",  "-----------",       "---------"))

    for(ii in 1:length(selected_loci)){
       if(is.na(selected_loci[ii])){
       message(sprintf("%15s  %10s        %10s        %8s           %-8.2f ",
        "Null Model", " ", " ", " ", extBIC[ii] ))
       }  else {
       message(sprintf("%15s  %10s        %10s       %8s            %-8.2f ",
        map[[1]][selected_loci[ii]], map[[2]][selected_loci[ii]], as.character(map[[3]][selected_loci[ii]]),
             selected_loci[ii], extBIC[ii] ))
     }  ## end if else 
   }
    message("\n\n")
 }




