#include <Rcpp.h>
#include <fstream>

RcppExport SEXP exportPLINKSet(SEXP geneSets_, SEXP fname_){
    Rcpp::List geneSets(geneSets_);
    std::string fname = Rcpp::as<std::string>(fname_); 
    const int geneSetNum = geneSets.size();
    std::vector<std::string> geneNames = geneSets.names();
    std::ofstream myfile;
    myfile.open(fname.c_str());
    for(int i=0;i<geneSetNum;i++) {
        myfile<<geneNames[i]<<std::endl;
        std::vector<std::string> rsids = geneSets[i];
        for (std::vector<std::string>::iterator it = rsids.begin() ; it != rsids.end(); ++it) {
            myfile<<*it<<std::endl;
        }
        myfile<<"END"<<std::endl<<std::endl;
    }
    myfile.close();
    
    return Rcpp::wrap(true);
}

