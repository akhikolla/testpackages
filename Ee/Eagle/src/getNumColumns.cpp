// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif







unsigned long getNumColumns(std::string fname) {
        unsigned long numcols = 0;

        std::string line;

        std::ifstream fileIN(fname.c_str());

        getline(fileIN, line);
        std::istringstream streamA(line);

        std::string token;

        numcols = line.length();

        return numcols;

}


