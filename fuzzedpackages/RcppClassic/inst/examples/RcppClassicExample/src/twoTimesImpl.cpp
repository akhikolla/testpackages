
#include <RcppClassic.h>

RcppExport SEXP twoTimesImpl(SEXP x){
    std::vector<int> iv = Rcpp::as<std::vector<int> >( x );
	for (size_t i=0; i<iv.size(); i++) {
        iv[i] = 2*iv[i];
    }
    RcppResultSet rs;
    rs.add("", iv);
    return(rs.getSEXP());
}
