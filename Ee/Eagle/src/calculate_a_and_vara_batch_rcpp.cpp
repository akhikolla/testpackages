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
Rcpp::List   calculate_a_and_vara_batch_rcpp(  long numreps, 
                                    Rcpp::CharacterVector f_name,
                                    Rcpp::NumericVector  selected_loci,
                                    Eigen::Map<Eigen::MatrixXd> inv_MMt_sqrt,
                                    Eigen::Map<Eigen::MatrixXd> dim_reduced_vara,
                                    double  max_memory_in_Gbytes,
                                    std::vector <long> dims,
                                    Eigen::Map<Eigen::MatrixXd> a,
                                    bool  quiet,
                                    Rcpp::Function message)


{
// Purpose: to calculate the untransformed BLUP (a) values from the 
//          dimension reduced BLUP value estimates. 
//          It is necessary to have a block multiplication form of this function. 
//          Also, since the matrix multiplications are reliant upon the BLAS library, only 
//          double precision matrix multiplication is possible. This means, the Mt matrix must 
//          be converted into a double precision matrix which has a large memory cost.  
// Note:
//      1. dims is the row, column dimension of the Mt matrix
//    dims[0] -----> number of SNP
//    dims[1] -----> number of genotypes


//Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>  (a), 




std::ostringstream
      os;

std::string
     fnamebin = Rcpp::as<std::string>(f_name);

 Eigen::MatrixXd
       ans(dims[0],numreps);

Eigen::MatrixXd
             ans_tmp;



Eigen::MatrixXd
    var_ans = Eigen::MatrixXd(dims[0],1);
    var_ans.setZero();





   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt)%*%M)
 double mem_bytes_needed =   ( 3.0   * (double) dims[1]  *   (double) dims[0] * sizeof(double))/1000000000.0;




if (!quiet){
//   Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
   message("Inside internal function calculate_a_and_vara_rcpp: Need memory (gigabytes)  ", mem_bytes_needed);
}



if(mem_bytes_needed < max_memory_in_Gbytes){


 // calculation will fit into memory
     Eigen::MatrixXd Mt = ReadBlockBin(fnamebin, 0, dims[1], dims[0]);


   if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(long ii=0; ii < selected_loci.size() ; ii++){
           Mt.row(selected_loci(ii)).setZero();
    }
   }



   // ans contains the BLUPs across the entire genome 
   // for all the reps (numreps). It is a dims[0] x numreps matrix.
     Eigen::MatrixXd  ans_part1;
    ans_part1.noalias() = inv_MMt_sqrt * a;
    ans.noalias() =   Mt  * ans_part1;

   // 0,1 matrix for if the a value is above the threshold. 
   // Using f to limit the number of var calculations needed. 
   Rcpp::IntegerMatrix f(dims[0], numreps);


  // Initialisation of f to 0
  #if defined(_OPENMP)
     #pragma omp for 
  #endif
 for(int j=0; j < numreps; j++){
    for(int i=0; i< dims[0] ; i++){
          f(i,j) = 0;
     }
  }

  #if defined(_OPENMP)
     #pragma omp for 
  #endif
  for(int j=0; j < numreps; j++){
    double maxval = (ans.col(j).array().abs()).maxCoeff();
    for(int i=0; i < dims[0]; i++)
       f(i,j) = ( std::abs(ans(i,j))  > maxval * 0.75 )  ;
  }


  Rcpp::Function w("which");
   Rcpp::IntegerVector inds =  w(f==1) ;  
  inds = inds - 1; // - 1 since R returns indexes starting from 1 instead of O like in Cpp
  for(int i=0; i<inds.size(); i++)
     inds(i) =  inds(i) % dims[0];

  Rcpp::Function uniq("unique");
  Rcpp::IntegerVector indx = uniq(inds);





 


   long  NumOfaAboveThreshold  = indx.size();

  // calculate untransformed variances of BLUP values
    Eigen::MatrixXd var_ans_tmp_part1;
    var_ans_tmp_part1.noalias() =   dim_reduced_vara * inv_MMt_sqrt;
    var_ans_tmp_part1 = inv_MMt_sqrt * var_ans_tmp_part1;

 Eigen::VectorXd ans1 ;


  #if defined(_OPENMP)
     #pragma omp for 
  #endif
 for(long i=0; i< NumOfaAboveThreshold ; i++){
       ans1 = (Mt.row(indx[i])) * var_ans_tmp_part1;
       var_ans(indx[i],0) =     ans1.dot(Mt.row(indx[i]) ) ;
  }



  return Rcpp::List::create(Rcpp::Named("a")=ans,
                            Rcpp::Named("vara") = var_ans);





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
                             ( 3.0  * (double) dims[1] *  sizeof(double) ) ;
      



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

      // find out the number of blocks needed
      long num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;
      message(" Block multiplication necessary for calculation of BLUPS and their variances. ");
      message(" Number of blocks needing block multiplication is ", num_blocks);
      message(" Number of rows in block is ", num_rows_in_block);

      // calculate untransformed variances of BLUP values
      Eigen::MatrixXd var_ans_tmp_part1;
      var_ans_tmp_part1.noalias() =   dim_reduced_vara * inv_MMt_sqrt;
      var_ans_tmp_part1 = inv_MMt_sqrt * var_ans_tmp_part1;

     // holds all the unique snp indexes
     Rcpp::IntegerVector indxAcrossBlocks;


      for(long i=0; i < num_blocks; i++){
           message("Performing block iteration ... " , i );
           long start_row1 = i * num_rows_in_block;
           long num_rows_in_block1 = num_rows_in_block;
           if ((start_row1 + num_rows_in_block1) > dims[0])
                num_rows_in_block1 = dims[0] - start_row1;


           Eigen::MatrixXd Mt = ReadBlockBin(fnamebin, start_row1, dims[1], num_rows_in_block1) ;



           Eigen::MatrixXd
              vt1,
              ans_tmp1;

           Eigen::MatrixXd   var_ans_tmp(num_rows_in_block1,1);
           var_ans_tmp.setZero();

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


           // block multiplication to form untransposed a
           ans_tmp.noalias()  =   inv_MMt_sqrt  * a ;
           ans_tmp = Mt * ans_tmp;


           // 0,1 matrix for if the a value is above the threshold. 
           // Using f to limit the number of var calculations needed. 
           Rcpp::IntegerMatrix f(num_rows_in_block1 , numreps);


           // Initialisation of f to 0
           #if defined(_OPENMP)
               #pragma omp for 
           #endif
           for(int j=0; j < numreps; j++){
              for(int i=0; i< num_rows_in_block1  ; i++){
                    f(i,j) = 0;
               }   
           }

           #if defined(_OPENMP)
               #pragma omp for 
           #endif
           for(int j=0; j < numreps; j++){
              double maxval = (ans_tmp.col(j).array().abs()).maxCoeff();
              for(int i=0; i < num_rows_in_block1 ; i++)
                 f(i,j) = ( std::abs(ans_tmp(i,j))  > maxval * 0.75  )  ;
           }




           Rcpp::Function w("which");
           Rcpp::IntegerVector inds =  w(f==1);
           inds = inds - 1; // - 1 since R returns indexes starting from 1 instead of O like in Cpp
           for(int i=0; i<inds.size(); i++)
               inds(i) =  inds(i) % num_rows_in_block1 ;



           Rcpp::Function uniq("unique");
           Rcpp::IntegerVector indx = uniq(inds);

           long  NumOfaAboveThreshold  = indx.size();


           Eigen::VectorXd ans1 ;



           #if defined(_OPENMP)
               #pragma omp for 
           #endif
          for(long i=0; i< NumOfaAboveThreshold ; i++){
            ans1 = (Mt.row(indx[i])) * var_ans_tmp_part1;
            var_ans_tmp(indx[i],0) =     ans1.dot(Mt.row(indx[i]) ) ;
           }

          // assign block vector results to final vector (ans) of results
          long counter;
          counter = 0;
          for(long j=start_row1; j < start_row1 + num_rows_in_block1; j++){
              ans.row(j)  = ans_tmp.row(counter);
              var_ans(j,0) = var_ans_tmp(counter,0);
              counter++;
           }

           if (!quiet  )  message( "block done ... " );


  } // end for long

  return Rcpp::List::create(Rcpp::Named("a")=ans,
                            Rcpp::Named("vara") = var_ans) ;



}  //  end if block update


}




