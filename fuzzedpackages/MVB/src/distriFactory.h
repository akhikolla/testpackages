#ifndef DISTRIFACTORY_H_
#define DISTRIFACTORY_H_

#include <map>
#include <string>
#include "gaussian.h"

namespace lps {
  typedef Loss* (*createLossFunc) (const arma::mat&, const arma::mat&);
  // factory class
  class DistriFactory {
  public:
    static DistriFactory& instance();
    void registerDistri(std::string, createLossFunc);
    Loss* createLoss(std::string, const arma::mat&, const arma::mat&);
    ~DistriFactory() {};

  private:
    std::map<std::string, createLossFunc> creatorFunctions;
    DistriFactory() {};
    DistriFactory(const DistriFactory&) {};
    DistriFactory& operator=(const DistriFactory&) {return *this;};
  };

  // template factory helper
  template <class T>
    class DistriHelper {
  public:
    DistriHelper(std::string);
    static Loss* create(const arma::mat&, const arma::mat&);
  };

  template<class T>
    Loss* DistriHelper<T>::create(const arma::mat& inputY, 
				  const arma::mat& inputX) {
    return new T(inputY, inputX);
  }

  template <class T>
    DistriHelper<T>::DistriHelper(std::string id) {
    DistriFactory& theLossFactory = DistriFactory::instance();
    theLossFactory.registerDistri(id, DistriHelper<T>::create);
  }
}

#endif // DISTRIFACTORY_H_

