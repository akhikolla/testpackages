// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "readblock.h"


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif






// ------------------------------------------------------
//    Calculation of untransformed BLUP a values 
// ------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List   calculate_a_and_vara_rcpp(  Rcpp::CharacterVector f_name,
                                    Rcpp::NumericVector  selected_loci,
                                    Eigen::Map<Eigen::MatrixXd> inv_MMt_sqrt,
                                    Eigen::Map<Eigen::MatrixXd> dim_reduced_vara,
                                    double  max_memory_in_Gbytes,
                                    std::vector <long> dims,
                                    Eigen::VectorXd  a,
                                    bool  quiet,
                                    Rcpp::Function message)
{
// Purpose: to calculate the untransformed BLUP (a) and var(a) values from the 
//          dimension reduced BLUP and var value estimates. 
//          It is necessary to have a block multiplication form of this function. 
//          Also, since the matrix multiplications are reliant upon the BLAS library, only 
//          double precision matrix multiplication is possible. This means, the Mt matrix must 
//          be converted into a double precision matrix which has a large memory cost.  
// Note:
//      1. dims is the row, column dimension of the Mt matrix

//    dims[0] -----> number of SNP
//    dims[1] -----> number of genotypes




std::ostringstream
      os;

std::string
     fnamebin = Rcpp::as<std::string>(f_name);

 Eigen::MatrixXd
        ans(dims[0],1);

Eigen::MatrixXd
             ans_tmp,
             var_ans_tmp(dims[0] , dims[1]);



Eigen::MatrixXd
    var_ans = Eigen::MatrixXd::Zero(dims[0],1);  // initialize to 0
    // var_ans = Eigen::MatrixXd(dims[0],1);




   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt)%*%M)
 double mem_bytes_needed =   ( 2.5   * (double) dims[1]  *   (double) dims[0] * sizeof(double))/1000000000.0;

if (!quiet){
//   Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
   message("Inside internal function calculate_a_and_vara_rcpp: Need memory (gigabytes)  ", mem_bytes_needed);
}



if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory
     Eigen::MatrixXd Mt = ReadBlockBin(fnamebin, 0, dims[1], dims[0]);



   if(!R_IsNA(selected_loci(0))){
   // setting rows to 0
   for(long ii=0; ii < selected_loci.size() ; ii++){
           Mt.row(selected_loci(ii)).setZero();
    }
   }


  

    Eigen::MatrixXd  ans_part1;
    ans_part1.noalias() = inv_MMt_sqrt * a;


    ans.noalias() =   Mt  * ans_part1;

   Rcpp::IntegerVector f(dims[0]);

  #if defined(_OPENMP)
     #pragma omp for 
  #endif
 for(int i=0; i< dims[0] ; i++){
       f[i] = 0;
  }


  // changed this back to calculating all values so that we would get a 
  // test statistic across the entire genome and better for interpretation. 
   // f = (ans.array().abs()  > (ans.array().abs().maxCoeff() * 0.75 ) ) ; 
   f = (ans.array().abs()  > (ans.array().abs().maxCoeff() * 0.0 ) ) ; 
      
   long  NumOfaAboveThreshold  = Rcpp::sum(f);
  Rcpp::NumericVector indx(NumOfaAboveThreshold);
   long counter=0;   
   for( long ii=0; ii< dims[0]; ii++){
      if (f[ii]==1){
         indx[counter] = ii;
         counter++;
      }
   } 

 



  // calculate untransformed variances of BLUP values
    Eigen::MatrixXd var_ans_tmp_part1;
    var_ans_tmp_part1.noalias() =   dim_reduced_vara * inv_MMt_sqrt;
    var_ans_tmp_part1 = inv_MMt_sqrt * var_ans_tmp_part1;



//  Eigen::MatrixXd var_ans_tmp_part1 =  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;a

 Eigen::VectorXd ans1 ;


  #if defined(_OPENMP)
     #pragma omp for 
  #endif
 for(long i=0; i< NumOfaAboveThreshold ; i++){
       ans1 = (Mt.row(indx[i])) * var_ans_tmp_part1;
       var_ans(indx[i],0) =     ans1.dot(Mt.row(indx[i]) ) ;
  }






} else {
    //  -----------------------------------------
    //       BLOCK WISE UPDATE
    //  -----------------------------------------


      // calculation being processed in block form
      message(" Increasing maxmemGb would improve performance... \n");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algebra but it is comes to the following
      // Strictly, 2 should be used here but I want some extra memory to play with 
      long num_rows_in_block =  max_memory_in_Gbytes * ( 1000000000.0) /
                             ( 2.5  * (double) dims[1] *  sizeof(double) ) ;


      if (num_rows_in_block < 0){
        message("\n");
        message( "Error:  availmemGb is set to " , max_memory_in_Gbytes );
        message( "        Cannot even read in a single row of data into memory." );
        message( "        Please increase availmemGb for this data set." );
        message("\n");
        message(" multiple_locus_am has terminated with errors\n" );

       return Rcpp::List::create(Rcpp::Named("a")=0,
                            Rcpp::Named("vara") = 0);

      }


      // blockwise multiplication
      Eigen::MatrixXd
              vt1;


      // find out the number of blocks needed
      long num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;
      if (!quiet  ){
      message(" Block multiplication necessary. \n");
      message(" Number of blocks needing block multiplication is ... % d \n", num_blocks);
      }

      // This originally sat inside the block loop but it can sit outside it instead. 
      // variance calculation
      // vt.noalias() =  Mtd *  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;
      vt1.noalias() =  dim_reduced_vara * inv_MMt_sqrt;
      vt1 =  inv_MMt_sqrt * vt1;



      for(long i=0; i < num_blocks; i++){
         message("Performing block iteration ... " , i );
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[0])
            num_rows_in_block1 = dims[0] - start_row1;

           Eigen::MatrixXd Mt = ReadBlockBin(fnamebin, start_row1, dims[1], num_rows_in_block1) ;




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


           //  ans_tmp  =  Mtd *  inv_MMt_sqrt  * a ;
             ans_tmp.noalias()  =   inv_MMt_sqrt  * a ;
             ans_tmp = Mt * ans_tmp;


            // Dont need this sitting inside this loop. Moved it outside the block loop
            // variance calculation
            // vt.noalias() =  Mtd *  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;
            // vt1.noalias() =  dim_reduced_vara * inv_MMt_sqrt;
            // vt1 =  inv_MMt_sqrt * vt1;




   Rcpp::NumericVector f(num_rows_in_block1);
   f = (ans_tmp.array().abs()  > (ans_tmp.array().abs().maxCoeff() * 0.75 ) ) ;

   long  NumOfaAboveThreshold  = Rcpp::sum(f);
  Rcpp::NumericVector indx(NumOfaAboveThreshold);
   long counter=0;
   for( long ii=0; ii< num_rows_in_block1 ; ii++){
      if (f[ii]==1){
         indx[counter] = ii;
         counter++;
      }
   }

 Eigen::MatrixXd   var_ans_tmp;
 var_ans_tmp = Eigen::MatrixXd::Zero(num_rows_in_block1,1);

 Eigen::VectorXd ans1(num_rows_in_block1);

#if defined(_OPENMP)
#pragma omp  for  
#endif
for(long i=0; i< NumOfaAboveThreshold ; i++){
       ans1 = (Mt.row(indx[i])) * vt1;
       var_ans_tmp(indx[i],0) =     ans1.dot(Mt.row(indx[i]) ) ;

}






            // assign block vector results to final vector (ans) of results
            counter = 0;
            for(long j=start_row1; j < start_row1 + num_rows_in_block1; j++){
                 ans(j,0) = ans_tmp(counter,0);
                 var_ans(j,0) = var_ans_tmp(counter,0);
                 counter++;
            }

             if (!quiet  )  message( "block done ... " );


      } // end for long



}  //  end if block update


  return Rcpp::List::create(Rcpp::Named("a")=ans,
                            Rcpp::Named("vara") = var_ans);


}



