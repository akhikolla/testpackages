// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "readblock.h"

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif






// [[Rcpp::export]]
Eigen::VectorXi  extract_geno_rcpp(Rcpp::CharacterVector f_name_ascii,
                                   double  max_memory_in_Gbytes,
                                    long  selected_locus,
                                    std::vector<long> dims)
{
  std::string
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

  long
     nind;

  nind = dims[0];



//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
double
  memory_needed_in_Gb =  ( (double) dims[0] *   (double) dims[1] *   sizeof(double) )/( (double) 1000000000.0) ;


Eigen::VectorXi
   column_of_genos(nind);



if(max_memory_in_Gbytes > memory_needed_in_Gb ){
   // reading entire data file into memory
     Eigen::MatrixXd genoMat =  ReadBlockBin(fnamebin,  0, dims[1], dims[0]);

   column_of_genos = genoMat.col(selected_locus).cast<int>() ;



}  else {
    long num_rows_in_block = ((double) max_memory_in_Gbytes  * (double) 1000000000.0 )/(sizeof(double) * (double) dims[1]);

         long num_blocks = dims[0]/num_rows_in_block;
          if (dims[0] % num_rows_in_block)
                 num_blocks++;


          for(long i=0; i < num_blocks; i++){
              long start_row1 = i * num_rows_in_block;
              long num_rows_in_block1 = num_rows_in_block;
              if ((start_row1 + num_rows_in_block1) > dims[0])
                     num_rows_in_block1 = dims[0] - start_row1  ;

              Eigen::MatrixXd
                genoMat_block1 ( ReadBlockBin(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;


              // dealing with assigning column_of_genos when some values 
              // may be missing due to having been removed. 
              long colindx = start_row1;
              for(long j=start_row1; j< start_row1+num_rows_in_block1 ; j++){
                  column_of_genos(colindx) = genoMat_block1.col(selected_locus)(j-start_row1);
                  colindx++;

              } // end for j

          } // end for  i


} // end if max_memory

return(column_of_genos);

}



