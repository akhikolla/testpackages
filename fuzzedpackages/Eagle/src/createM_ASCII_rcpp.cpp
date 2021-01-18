// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include  "createM_ASCII_rcpp.h"


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif







// [[Rcpp::export]]
bool  createM_ASCII_rcpp(Rcpp::CharacterVector f_name, Rcpp::CharacterVector f_name_ascii,
                  Rcpp::CharacterVector  type,
                  std::string AA,
                  std::string AB, 
                  std::string BB,
                  double  max_memory_in_Gbytes,  std::vector <long> dims,
                  bool quiet,
                  Rcpp::Function message, 
                  std::string missing)
{
  // Rcpp function to create space-removed ASCII file from ASCII and PLINK input files




// size_t found;


std::string
   line;


std::ofstream
   fileOUT;

//int 
//   genoval;

std::string
     ftype = Rcpp::as<std::string>(type),
     fname = Rcpp::as<std::string>(f_name),
     fnameascii = Rcpp::as<std::string>(f_name_ascii);



//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
double
  memory_needed_in_Gb;
  if (ftype == "PLINK"  ){
  // this is a PLINK ped file. Hence, we need to adjust the dims[1] to get the 
  // size of the genotype file in R land. 
    memory_needed_in_Gb =  ( (double) dims[0] *  ( (double) dims[1]-6.0)/2.0  *   sizeof(double) )/( (double) 1000000000.0) ;
  } else {
    // text file
    memory_needed_in_Gb =  ((double) dims[0] *  (double) dims[1] *   sizeof(double) )/( (double) 1000000000.0) ;
  }
  if ( ftype == "PLINK"  ){
     //------------------------------------
     // convert PLINK ped file into ASCII file with no spaces
     //----------------------------------------------
       bool it_worked = CreateASCIInospace_PLINK(fname, fnameascii, dims, quiet, message);
       if (!it_worked) // an error has occurred in forming ascii file
                 return false;

   }  else {
      //-------------------------------------------
      // convert text file into ASCII file
      //-----------------------------------------
      // Here, we do not need to worry about the amount of memory because 
      // we are processing a line of the file at a time. This is not the case when 
      // creating a ASCII Mt because we have to read in blocks before we can 
      // transpose. 
      if (!quiet  )
          message(" A text file is being assumed as the input data file type. ");

       if ( 1.0 * memory_needed_in_Gb   > max_memory_in_Gbytes){
           bool it_worked = CreateASCIInospace(fname, fnameascii, dims, AA, AB, BB, quiet, message, missing);
           if (!it_worked) // an error has occurred in forming ascii file
                 return false;
       } else {
           bool it_worked =  CreateASCIInospace(fname, fnameascii, dims, AA, AB, BB, quiet, message, missing);
          if (!it_worked) // an error has occurred in forming ascii file
                 return false;
       }


   }  // end if type == "PLINK" 


  return true;


}




