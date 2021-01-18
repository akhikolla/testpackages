

# gets the function name from the code
get_function_name<-function(code)
{

preamble<-'// [[Rcpp::depends(RcppEigen)]]'
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'// [[Rcpp::depends(RcppEigenAD)]]',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'#include <Rcpp.h>')
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'#include <RcppEigen.h>')
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'typedef Eigen::MatrixXd ADmat;',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'// [[Rcpp::export]]',sep="")

    

eigencode<-paste(preamble,code)

return(sourceCpp(code=eigencode)[[1]])

}

