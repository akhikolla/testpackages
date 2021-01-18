// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif





unsigned long getNumRows(std::string fname) {
        unsigned long numrows = 0;

        std::string line;

        std::ifstream fileIN(fname.c_str());

        std::istringstream streamA(line);

        // Determine number of columns in file
        fileIN.clear(); // returns to beginning of line
        fileIN.seekg(0, std::ios::beg);

        // Determine number of rows in file
                while(fileIN.good()){
                        while(getline(fileIN, line )){
                                numrows++;
                        }
                }

        return numrows;

}



