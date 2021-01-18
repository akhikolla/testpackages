#ifndef wpik_H
#define wpik_H


#include "distUnitk.h"
#include <RcppArmadillo.h>

arma::sp_mat wpik(arma::mat X,
                  arma::vec pik,
                  double bound,
                  bool tore,
                  bool shift,
                  double toreBound);



#endif
