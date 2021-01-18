#ifndef distUnitk_H
#define distUnitk_H

#include <RcppArmadillo.h>

arma::vec distUnitk(arma::mat X,
                    int k,
                    bool tore,
                    double toreBound);

#endif
