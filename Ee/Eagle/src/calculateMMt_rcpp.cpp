// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "readblock.h"

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


using Eigen::SelfAdjointEigenSolver;



// [[Rcpp::export]]
Eigen::MatrixXd  calculateMMt_rcpp(Rcpp::CharacterVector f_name,
                                   double  max_memory_in_Gbytes, int num_cores,
                                   Rcpp::NumericVector  selected_loci , std::vector<long> dims,
                                   bool  quiet, Rcpp::Function message)
{
// set multiple cores
Eigen::initParallel();

#ifdef _OPENMP
#include <omp.h>
   omp_set_num_threads(num_cores);
#endif




Eigen::setNbThreads(num_cores);
message(" Number of cores being used for calculation is .. ", num_cores);


std::string
   line;


std::ofstream
   fileOUT;

//int 
//   genoval;

std::string
     fnamebin = Rcpp::as<std::string>(f_name);

// gpu will only work with double precision matrices in Eigen. 
// Had to change code to be double precision. 
Eigen::MatrixXd
    MMt(Eigen::MatrixXd(dims[0], dims[0]).setZero());




//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
// Memory required for 
// MMt   dims[0] * dims[0] *  sizeof(double) 
// genoMat  dims[0] * dims[1] * sizeof(double) 
// genoMat transpose dims[0] * dims[1] * sizeof(double) 
// 
// Block update
//
// MMt   dims[0] * dims[0] *  sizeof(double) 
// genoMat  num_rows_in_block * dims[1] * sizeof(double) 
// genoMat transpose num_rows_in_block * dims[1] * sizeof(double) 

double
  memory_needed_in_Gb =   3  *(  dims[0] *    dims[1] *   sizeof(double) )/( (double) 1000000000.0) ;   




//-------------------------
// Perform MMt calculation
//-------------------------
if(max_memory_in_Gbytes > memory_needed_in_Gb ){
   // reading entire data file into memory
     Eigen::MatrixXd genoMat = ReadBlockBin(fnamebin,  0, dims[1], dims[0] );
  //   Eigen::MatrixXd genoMat = ReadBlockFast(fnamebin,  0, dims[1], dims[0] );
   if(!R_IsNA(selected_loci(0))){
     // setting columns to 0
     for(long ii=0; ii < selected_loci.size() ; ii++)
       genoMat.col(selected_loci(ii)).setZero();
   }




     MMt.noalias() = genoMat * genoMat.transpose();

     



} else {
    // based on user defined memory. Doing MMt via blockwise multiplication

     // max_memory_in_Gbytes = 2.5 *(  num_rows_in_block  *    dims[1] *   sizeof(double) )/( (double) 1000000000.0) ;
        long num_rows_in_block = max_memory_in_Gbytes * 1000000000.0 / ( 3  * dims[1] *  sizeof(double) );

    // double part1 = -2 *  dims[1];
    // double part2 = 4.0 * (double) dims[1] *  (double) dims[1] +  4.0 * (double) max_memory_in_Gbytes  * 1000000000.0/sizeof(double);
    // part2 = sqrt(part2);
    // long num_rows_in_block = (part1 + part2)/2.2;

           // blockwise multiplication

          // find out the number of blocks needed
          long num_blocks = dims[0]/num_rows_in_block;



          if (dims[0] % num_rows_in_block)
                 num_blocks++;
   if(!quiet){
    message( "   Number of rows in block for MMt calculation is " , num_rows_in_block);
    message( "   Number of blocks in MMt calculation is " , num_blocks);
  }

          for(long i=0; i < num_blocks; i++){
              long start_row1 = i * num_rows_in_block;
              long num_rows_in_block1 = num_rows_in_block;
              if ((start_row1 + num_rows_in_block1) > dims[0])
                     num_rows_in_block1 = dims[0] - start_row1;

                Eigen::MatrixXd
                     genoMat_block1 ( ReadBlockBin(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;



              Eigen::MatrixXd
                   MMtsub(Eigen::MatrixXd(num_rows_in_block1, num_rows_in_block1).setZero());

             if(!R_IsNA(selected_loci(0) )){
             // setting columns to 0
             for(long ii=0; ii < selected_loci.size() ; ii++)
                genoMat_block1.col(selected_loci(ii)).setZero();
             }
              MMtsub.noalias() = genoMat_block1 * genoMat_block1.transpose();
              //          i            j            num rows               num   cols
              MMt.block(start_row1, start_row1, num_rows_in_block1, num_rows_in_block1) = MMtsub;

              for(long j=i+1;j<num_blocks; j++){
                   long start_row2 = j * num_rows_in_block;
                   long num_rows_in_block2 = num_rows_in_block;
                   if ((start_row2 + num_rows_in_block2) > dims[0])
                          num_rows_in_block2 = dims[0] - start_row2;
                     Eigen::MatrixXd
                        genoMat_block2 ( ReadBlockBin(fnamebin,  start_row2, dims[1], num_rows_in_block2)) ;




                   Eigen::MatrixXd    MMtsub(Eigen::MatrixXd(num_rows_in_block1, num_rows_in_block2).setZero());

                  if(!R_IsNA(selected_loci(0) )){
                   // setting columns to 0
                   for(long jj=0; jj < selected_loci.size() ; jj++)
                      genoMat_block2.col(selected_loci(jj)).setZero();
                    
                   }
                   MMtsub.noalias() = genoMat_block1 * genoMat_block2.transpose();
                   //          i,        j,     num rows,              num cols
                   MMt.block(start_row1, start_row2, num_rows_in_block1, num_rows_in_block2) = MMtsub;
                   // and its symmetric block
                   MMt.block(start_row2, start_row1, num_rows_in_block2, num_rows_in_block1) = MMtsub.transpose();


            }  // end for int j





          } // end for int


 // }  // end inner if else

}  // end outer if else



  return MMt;

}


