//#include <Rcpp.h>
#include "distriFactory.h"

namespace lps {
  void DistriFactory::registerDistri(std::string distriID,
				     createLossFunc func) {
    creatorFunctions.insert(std::pair<std::string, createLossFunc>
			    (distriID, func));
  }

  Loss* DistriFactory::createLoss(std::string distriID,
				  const arma::mat& inputY,
				  const arma::mat& inputX) {
    // convert input string to be consistently in lower case
    for (unsigned i = 0; i < distriID.size(); i++) 
      distriID[i] = tolower(distriID[i]);
    std::map<std::string, createLossFunc>::const_iterator
      ptr = creatorFunctions.find(distriID);
    if (ptr == creatorFunctions.end()) {
      Rcpp::Rcout << distriID << " is an unknown distribution" << std::endl;
      return NULL;
    }
    return (ptr -> second)(inputY, inputX);
  }

  DistriFactory& DistriFactory::instance() {
    static DistriFactory theFactory;
    return theFactory;
  }	  
}

