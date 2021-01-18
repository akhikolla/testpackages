// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif








//get number of rows and columns in marker file
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// [[Rcpp::export]]
std::vector<long>   getRowColumn(std::string fname)
{
  // Purpose:  to open the marker file where the marker data are kept.
  //           An error will be produced if the file cannot be found.
  //           I am assuming no row or column names

 std::string
   line;

 std::ostringstream
      os;


 std::vector<long> dimen(2,0)  ;  // dim[0] row number
                               // dim[1] col number 

 // open file and check for its existence. 
 std::ifstream fileIN(fname.c_str());
 if(!fileIN.good()) {
      os << "\n\n ERROR: Could not open  " << fname << "\n\n" << std::endl;
      Rcpp::stop(os.str() );
 }


 // Determine number of rows in file
 while(fileIN.good()){
      while(getline(fileIN, line )){
         dimen[0]++;
      }
 }


 // Determine number of columns in file
 fileIN.clear(); // returns to beginning of line
 fileIN.seekg(0, std::ios::beg);

 getline(fileIN, line );
 std::istringstream streamA(line);

  std::string
    token ;

 while(streamA >> token)
   dimen[1]++;
 fileIN.close();
 if(fileIN.bad())
 {
    os << "\n\nERROR:  There was a problem with reading the marker file - possibly strange ASCII characters.\n\n";
    Rcpp::stop(os.str() );
}
return dimen;

}

