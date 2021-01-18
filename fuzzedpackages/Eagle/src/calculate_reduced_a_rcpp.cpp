// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "readblock.h"


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif





//--------------------------------------------
// Calculation of transformed blup a values
//--------------------------------------------
// [[Rcpp::export]]
Eigen::MatrixXd calculate_reduced_a_rcpp ( Rcpp::CharacterVector f_name_ascii, double varG,
                                           Eigen::Map<Eigen::MatrixXd> P,
                                           Eigen::Map<Eigen::MatrixXd>  y,
                                           double max_memory_in_Gbytes,
                                           std::vector <long> dims,
                                           Rcpp::NumericVector  selected_loci,
                                           bool quiet,
                                           Rcpp::Function message)
{
  // function to calculate the BLUPs for the dimension reduced model. 
  // It is being performed in Rcpp because it makes use of Mt. 
  // Args
  // f_name_ascii    path + file name of Mt.bin
  // varG          variance of polygenic component
  // P             calculate in R
  // y             response/trait  but read in as a row matrix
  // max_memory_in_Gbytes  working memory in gigabytes
  // dims          dimension (row, column), of M.

std::string
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

Eigen::MatrixXd
      ar(dims[1],1);  // column vector

std::ostringstream
      os;

Eigen::MatrixXd
      nullmat = Eigen::MatrixXd::Zero(1,1);  // 0 matrix of size 1,1 for null purposes 


// const size_t bits_in_double = std::numeric_limits<double>::digits;


   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt))
double mem_bytes_needed =   ( (double) dims[0]* (double) dims[1] + (double) dims[0]* (double) dims[0] +  (double) dims[0] ) *  ( sizeof(double)/( 1000000000.0));

if (!quiet){
    message("Inside internal function calculate_reduced_a_rcpp. Memory needed (gigabytes): ", mem_bytes_needed);
    message("Inside internal function calculate_reduced_a_rcpp. Memory available (gigabytes): ", max_memory_in_Gbytes);


}

if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory

   Eigen::MatrixXd
                   Mt;

   Mt = ReadBlockBin(fnamebin, 0, dims[0], dims[1]);
 //  Mt = ReadBlockFast(fnamebin, 0, dims[0], dims[1]);

  if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(long ii=0; ii < selected_loci.size() ; ii++)
          Mt.row(selected_loci(ii) ).setZero();
   }


   // ar  =    varG * Mt *  P   * y ;
    ar  =     P   * y ;
    ar  =    Mt * ar;
    ar  =    varG * ar;
} else {

      // calculation being processed in block form
      message(" Note:  Increasing availmemGb would improve performance... ");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algrebra but it is comes to the following
      long num_rows_in_block = (max_memory_in_Gbytes * 1000000000.0/sizeof(double) - (double) dims[0] * (double) dims[0] - (double) dims[0])/( (double) dims[0] ) ;

    if (num_rows_in_block < 0){
        message("\n");
        message("Error:  availmemGb is set to " , max_memory_in_Gbytes );
        message("        Cannot even read in a single row of data into memory." );
        message("        Please increase availmemGb for this data set." );
        message("\n");
        message( "AM has terminated with errors\n" );
        return nullmat;

      }


      // blockwise multiplication

      // find out the number of blocks needed
      long num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;

      if (!quiet ){
      message(" Block multiplication necessary. \n");
      message(" Number of blocks needing block multiplication is ... % d \n", num_blocks);
      }
      for(long i=0; i < num_blocks; i++){
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[1])
            num_rows_in_block1 = dims[1] - start_row1;

          Eigen::MatrixXd
                  Mt;
            Mt = ReadBlockBin(fnamebin, start_row1, dims[0], num_rows_in_block1) ;
        //   Mt = ReadBlockFast(fnamebin, start_row1, dims[0], num_rows_in_block1) ;

         Eigen::MatrixXd
             ar_tmp;

         if(!R_IsNA(selected_loci(0))){
         // setting columns (or row when Mt) to 0
            for(long ii=0; ii < selected_loci.size() ; ii++)
            {
            // since we are now dealing with Mt, and blocking on columns, 
            // because columns are rows in Mt, then we have to be careful
            // that we do not select loci outside the block bounds. Also 
            // the values have to be adjusted based on the block number
                if(selected_loci(ii) >= start_row1 && selected_loci(ii) < start_row1 + num_rows_in_block1 )
                {   // selected loci index is in block 
                long block_selected_loci = selected_loci(ii) - start_row1;
                Mt.row(block_selected_loci).setZero();
                }
             }
         }

         // ar_tmp  =  varG * Mt *  P  * y ;
          ar_tmp = P * y;
          ar_tmp = Mt * ar_tmp;
          ar_tmp = varG * ar_tmp;



         // assign block vector results to final vector (ar) of results
         long  counter = 0;
         for(long j=start_row1; j < start_row1 + num_rows_in_block1; j++){
              ar(j,0) = ar_tmp(counter,0);
              counter++;
         }

       if (!quiet  )  message( "block done ... ");
      } // end for long




}  // end if mem_bytes_needed

  return(ar);

} // end function 




