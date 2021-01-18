// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif






// [[Rcpp::export]]
Eigen::MatrixXd  ReadBlockBin(std::string binfname,
                           long start_row,
                           long numcols,
                           long numrows_in_block)

{
 // reads in data from ASCII file 
 // to form M Eign double matrix 
std::ostringstream
      os;
std::string
   line;




Eigen::MatrixXd
      M(numrows_in_block, numcols) ;

long   numelem = numrows_in_block * numcols;

char * buffer ;
buffer = new char[numelem];


// Binary file
   std::ifstream fileIN(binfname.c_str(), std::ios::in | std::ios::binary);

   if(!fileIN.good()) {
      os << "ERROR: Could not open  " << binfname << std::endl;
      Rcpp::stop(os.str() );
   }

 // set pointer to beginning of block of binary data
 long begofblock = start_row *  numcols;  // number of characters to skip 
 fileIN.seekg(begofblock, std::ios::beg);

 // read block into buffer
 // writing vector to binary file
   fileIN.read( buffer, numelem );
 
  // check 
 #pragma omp for
 for(long ii=0; ii < numrows_in_block; ii++){
  for(long jj=0; jj < numcols ; jj++){


    // M(start_row + ii , jj) = (double) ((buffer[ ii*numcols + jj ] - '0') - 1);  // converting data into -1, 0, 1 
    M( ii , jj) = (double) ((buffer[ ii*numcols + jj ] - '0') - 1);  // converting data into -1, 0, 1 
  } 
 }



// Close the bin file
   fileIN.close();

delete [] buffer;
buffer = NULL;

 return M;

}



