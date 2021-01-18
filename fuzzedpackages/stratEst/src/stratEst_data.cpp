//#define ARMA_NO_DEBUG
#define ARMA_USE_CXX14
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec stratEst_data_cpp( arma::vec id, arma::vec game, arma::vec period, arma::vec input, arma::vec lagged_input, int lag, int num_ids ) {
    for(int i = 1; i <= num_ids; i++) {
      arma::vec games = unique( game( find( id == i ) ) );
      int num_games = games.n_elem;
      for(int j = 1; j <= num_games; j++) {
        int max_period = max( period( find( game == j && id == i ) ) );
        for( int k = 1; k <= max_period; k++){
          if( lag < k ){
            lagged_input( find( game ==j && id == i && period == k ) ) = input( find( game ==j && id == i && period == (k-lag) ) );
          }
        }
      }
    }
  return(lagged_input);
}
