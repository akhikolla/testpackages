#include <Rcpp.h>
#include "ETT_R.cpp"
#include <memory>
using namespace Rcpp;
using namespace std;

RCPP_EXPOSED_AS(ETT_R_Wrapper)
  
RCPP_MODULE(ETT) {
  Rcpp::class_<ETT_R_Wrapper>("ETT_R_Wrapper")
    .method("process",&ETT_R_Wrapper::processForR,"ETT processing entry")
    .method("getMachineIdentifiers",&ETT_R_Wrapper::getMachineIdentifiers,"Get the list of machine identifiers")
    .method("getCoincidenceMatrix",&ETT_R_Wrapper::getCoincidenceMatrix,"Get a coincidence matrix for the machine")
    .method("getCoincidenceValues",&ETT_R_Wrapper::getCoincidenceValues,"Get a coincidence list values for the machine")
    .method("printMachines",&ETT_R_Wrapper::printMachinesForR,"Print machine details")
    .method("mergeAllMachines",&ETT_R_Wrapper::mergeAllMachinesForR,"Merge all machines")
    .method("induceSubmachine",&ETT_R_Wrapper::induceSubmachine,"Induce submachines")
    .method("setStatePattern",&ETT_R_Wrapper::setStatePatternForR,"Set state pattern")
    .method("setTransitionPattern",&ETT_R_Wrapper::setTransitionPatternForR,"Set transition pattern")
    .method("cleanMachineKeys",&ETT_R_Wrapper::cleanMachineKeysForR,"Clean machine keys")
    .method("clone",&ETT_R_Wrapper::cloneForR,"Clone this wrapper")
    .method("compressMachines",&ETT_R_Wrapper::compressMachinesForR,"Compress machines")
    .method("serialize_ETT",&ETT_R_Wrapper::serialize_ETT,"Serialize");

  Rcpp::function("deserialize_ETT",&ETT_R_Wrapper::deserialize_ETT);
  Rcpp::function("create_ETT_wrapper",&ETT_R_Wrapper::create_ETTWrapper);
  Rcpp::function("echo_rcpp",&ETT_R_Wrapper::echo);
}
