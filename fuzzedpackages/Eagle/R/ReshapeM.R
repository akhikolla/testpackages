ReshapeM  <- function(fnameM, fnameMt, indxNA, dims){
   ## function to create a temp version of M.ascii and Mt.ascii where the rows and columns, 
   ## respectively have been removed for the elements in indxNA

   ## its indxNA-1 so that indexes start from 0 as in c++
   res <- ReshapeM_rcpp(fnameM=fnameM, fnameMt=fnameMt, indxNA=(indxNA-1), dims=dims)
   return(res)  ## returns integer vector with new dims of reshaped matrix M
}


