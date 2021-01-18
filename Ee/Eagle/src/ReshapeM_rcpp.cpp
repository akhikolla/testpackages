// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


// This code needs to be modified to be binary compliant. 
// Read a line buffer at a time, and write result to new temporary file.





// [[Rcpp::export]]
std::vector <long>    ReshapeM_rcpp( Rcpp::CharacterVector  fnameM,
               Rcpp::CharacterVector  fnameMt,
               std::vector <long> indxNA,
               std::vector <long> dims){


  // note indxNA starts from 0

  std::vector <long>
        newdims(2,0);

 std::ostringstream
      os;


   std::string
       line,
       FnameM = Rcpp::as<std::string>(fnameM),
       FnameMt = Rcpp::as<std::string>(fnameMt);



   //-------------------------------------------
   // converting M.ascii to reshaped M.asciitmp
   //-------------------------------------------
   // open file and check for its existence. 
   std::ifstream fileIN(FnameM.c_str());
   if(!fileIN.good()) {
        os << "\n\nERROR: Could not open  " << FnameM << "\n\n" << std::endl;
        Rcpp::stop(os.str() );
   }


   // change name for new no-space ASCII file with rows matching indxNA removed
   FnameM.append("tmp");


   // open ascii file that is to hold no-spaces genotype data
   std::ofstream fileOUT(FnameM.c_str(), std::ios::out );

   long rownum=0;
   bool writeline;
   while(fileIN.good()){
      while(getline(fileIN, line)){
          writeline = true;
          for(unsigned long ii=0; ii<indxNA.size(); ii++){
            if(indxNA[ii] == rownum)
                writeline = false;
          }
          if(writeline){
              fileOUT << line << std::endl;
              newdims[0]++;
          }
          rownum++;

      newdims[1] = line.length();

      }  // end inner while

   }  // end outer while(fileIN

 fileIN.close();
 fileOUT.close();

  //-------------------------------------------
  // converting Mt.ascii to reshaped Mt.asciitmp
  //-------------------------------------------

  // open file and check for its existence. 
  std::ifstream fileINt(FnameMt.c_str());
  if(!fileINt.good()) {
        os << "\n\nERROR: Could not open  " << FnameMt << "\n\n" << std::endl;
        Rcpp::stop(os.str() );
  }


  // change name for new no-space ASCII file with rows matching indxNA removed
  FnameMt.append("tmp");


   // open ascii file that is to hold no-spaces genotype data
   std::ofstream fileOUTt(FnameMt.c_str(), std::ios::out );

  while(fileINt.good()){
      while(getline(fileINt, line)){
          // removing columns
          // this is okay since indxNA is in decreasing size
          for(unsigned long ii=0; ii<indxNA.size(); ii++){
            line.erase ( indxNA[ii], 1 );
          }  // end for a
          fileOUTt << line << std::endl;
      }  // end inner while
 }  // end outer while(fileIN

 fileINt.close();
 fileOUTt.close();
  return newdims;

}  //end ReshapeM



