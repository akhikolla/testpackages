#include "Grid2D.h"
#include "Grid1D.h"

Grid2D::Grid2D(const arma::mat& Xi, const arma::vec& yi, const GridParams& PGi)
{
    // automatically selects lambda_0 (but assumes other lambdas are given in PG.P.ModelParams)
    X = &Xi;
    y = &yi;
    p = Xi.n_cols;
    PG = PGi;
    G_nrows = PG.G_nrows;
    G_ncols = PG.G_ncols;
    G.reserve(G_nrows);
    Lambda2Max = PG.Lambda2Max;
    Lambda2Min = PG.Lambda2Min;
    LambdaMinFactor = PG.LambdaMinFactor;

    P = PG.P;
}


Grid2D::~Grid2D(){
    delete Xtr;
    if (PG.P.Specs.Logistic){delete PG.P.Xy;}
    if (PG.P.Specs.SquaredHinge){delete PG.P.Xy;}
}

std::vector< std::vector<std::unique_ptr<FitResult> > > Grid2D::Fit()
{

    arma::vec Xtrarma;
    if (PG.P.Specs.Logistic)
    {
      auto n = X->n_rows;
      double b0 = 0;
      arma::vec ExpyXB =  arma::ones<arma::vec>(n);
      if (PG.intercept){
        for (unsigned int t = 0; t < 50; ++t){
          double partial_b0 = - arma::sum( *y / (1 + ExpyXB) );
          b0 -= partial_b0 / (n * 0.25); // intercept is not regularized
          ExpyXB = arma::exp(b0 * *y);
        }
      }
      PG.P.b0 = b0;
      Xtrarma = arma::abs(- arma::trans(*y /(1+ExpyXB)) * *X).t(); // = gradient of logistic loss at zero
      //Xtrarma = 0.5 * arma::abs(y->t() * *X).t(); // = gradient of logistic loss at zero

      arma::mat Xy =  X->each_col() % *y;
      PG.P.Xy = new arma::mat;
      *PG.P.Xy = Xy;
    }

    else if (PG.P.Specs.SquaredHinge)
    {
      auto n = X->n_rows;
      double b0 = 0;
      arma::vec onemyxb =  arma::ones<arma::vec>(n);
      arma::uvec indices = arma::find(onemyxb > 0);
      if (PG.intercept){
        for (unsigned int t = 0; t < 50; ++t){
          double partial_b0 = arma::sum(2 * onemyxb.elem(indices) % (- y->elem(indices) ) );
          b0 -= partial_b0 / (n * 2); // intercept is not regularized
          onemyxb = 1 - (*y * b0);
          indices = arma::find(onemyxb > 0);
        }
      }
      PG.P.b0 = b0;
      Xtrarma = 2 * arma::abs(arma::trans(y->elem(indices) % onemyxb.elem(indices))* X->rows(indices)).t(); // = gradient of loss function at zero
      //Xtrarma = 2 * arma::abs(y->t() * *X).t(); // = gradient of loss function at zero
      arma::mat Xy =  X->each_col() % *y;
      PG.P.Xy = new arma::mat;
      *PG.P.Xy = Xy;
    }

    else
    {
        Xtrarma = arma::abs(y->t() * *X).t();
    }


    double ytXmax = arma::max(Xtrarma);

    unsigned int index;
    if (PG.P.Specs.L0L1)
    {
        index = 1;
        if(G_nrows != 1)
        {
            Lambda2Max = ytXmax;
            Lambda2Min = Lambda2Max * LambdaMinFactor;
        }

    }
    else if (PG.P.Specs.L0L2) {index = 2;}

    arma::vec Lambdas2 = arma::logspace(std::log10(Lambda2Min), std::log10(Lambda2Max), G_nrows);
    Lambdas2 = arma::flipud(Lambdas2);

    std::vector<double> Xtrvec = arma::conv_to< std::vector<double> >::from(Xtrarma);

    Xtr = new std::vector<double>(X->n_cols); // needed! careful


    PG.XtrAvailable = true;

    for(unsigned int i=0; i<Lambdas2.size();++i) //auto &l : Lambdas2
    {
        *Xtr = Xtrvec;

        PG.Xtr = Xtr;
        PG.ytXmax = ytXmax;

        PG.P.ModelParams[index] = Lambdas2[i];
        if (PG.LambdaU == true)
            PG.Lambdas = PG.LambdasGrid[i];

        //std::vector<std::unique_ptr<FitResult>> Gl();
        //auto Gl = Grid1D(*X, *y, PG).Fit();
        G.push_back(std::move(Grid1D(*X, *y, PG).Fit()));
    }

    return std::move(G);

}
