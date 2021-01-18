#ifndef wpikInv_H
#define wpikInv_H

#include <RcppArmadillo.h>
#include "distUnitk.h"

arma::sp_mat wpikInv(arma::mat X,
                  arma::vec pik,
                  bool tore,
                  bool shift,
                  double toreBound);



#endif
