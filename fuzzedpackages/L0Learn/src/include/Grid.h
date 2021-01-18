#ifndef GRID_H
#define GRID_H
#include <set>
#include "GridParams.h"
#include "FitResult.h"


class Grid
{
    private:
        arma::mat Xscaled;
        arma::vec yscaled;
        arma::vec BetaMultiplier;
        arma::vec meanX;
        double meany;

    public:

        GridParams PG;

        std::vector< std::vector<double> > Lambda0;
        std::vector<double> Lambda12;
        std::vector< std::vector<unsigned int> > NnzCount;
        std::vector< std::vector<arma::sp_mat> > Solutions;
        std::vector< std::vector<double> >Intercepts;
        std::vector< std::vector<bool> > Converged;

        Grid(const arma::mat& X, const arma::vec& y, const GridParams& PG);
        //~Grid();

        void Fit();

};

#endif
