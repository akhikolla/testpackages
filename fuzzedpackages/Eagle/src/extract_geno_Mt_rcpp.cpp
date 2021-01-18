// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "readblock.h"



/* Originally, I was extracting columns of data from M. This required for the entire matrix to be read in. 
   Here, I am using my brain. I am using Mt where I just need to extract a single column. Much smarter :-)
*/

// [[Rcpp::export]]
Eigen::VectorXi  extract_geno_Mt_rcpp(Rcpp::CharacterVector f_name,
                                    long  selected_locus,
                                    std::vector<long> dims)
{
  std::string
     fnamebin = Rcpp::as<std::string>(f_name);   // name of Mt bin file

  long
     nind;
  nind = dims[1];


// dims[0]   -- p SNPs
// dims[1]   -- n genotypes individuals

   Eigen::VectorXi
   row_of_genos(nind);

   // row_of_genos = ReadBlockBin(fnamebin, selected_locus, dims[1], 1);
   Eigen::MatrixXd genoMat =  ReadBlockBin(fnamebin,  selected_locus, dims[1], 1);

    

   row_of_genos = genoMat.row(0).cast<int>() ;


   return(row_of_genos);

}



