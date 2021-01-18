//createM_BIN_rcpp.h

#ifndef createM_BIN_rcpp_h_INCLUDED   // if x.h hasn't been included yet...
#define createM_BIN_rcpp_h_INCLUDED   //   #define this so the compiler knows it has been included



bool  CreateBIN(std::string fname, std::string binfname, std::vector<long> dims,
                         std::string  AA,
                         std::string AB,
                         std::string BB,
                         bool  quiet,
                         Rcpp::Function message,
                         std::string missing);



bool  CreateBIN_PLINK(std::string fname, std::string binfname, std::vector<long> dims,
                         bool quiet, Rcpp::Function message);


#endif 
