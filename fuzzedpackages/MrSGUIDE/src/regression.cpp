//
// Created by Peigen Zhou on 7/15/18.
//

#include "regression.hpp"

namespace SubGuide {
namespace RegSol {

    using arma::ones;

    void RegParm::update(const int& n_, const int& p_)
    {
        n = n_;
        p = p_;
        df = n - p;
        BIC = arma::datum::inf;
        loss = arma::datum::inf;
        beta.set_size(p);
        beta.zeros();
        betaSE.set_size(p);
        betaSE.zeros();
    }

    // Private functions

    RegParm LinReg::fit(const mat& X, const vec& Y)
    {
        // assert(arma::is_finite(Y));
        // assert(arma::is_finite(X));

        auto& n = X.n_rows;
        auto& p = X.n_cols;
        RegParm result(n, p);
        if (arma::cond(X) > 100) return result;
        result.beta = arma::solve(X, Y);
        const vec& Ypred = this->predict(X, result.beta);
        result.loss = this->getLoss(Y, Ypred);

        result.BIC = (double)n * log(result.loss / (double)n) + (double)p * log((double)n);
        return result;
    }

    vec LinReg::predict(const mat& X, const vec& beta)
    {
        if (X.n_cols != beta.n_elem) {
            cerr << "dimension mismatch in LinReg::predict \n";
            // abort();
        }
        return X * beta;
    }

    double LinReg::getLoss(const vec& Y, const vec& Ypred)
    {
        const vec& residuals = Y - Ypred;
        return std::inner_product(residuals.begin(), residuals.end(), residuals.begin(), 0.0);
    }

    vec LinReg::calStdErr(const mat& X, const vec& Y)
    {
        this->parm = fit(X, Y);
        return arma::sqrt(parm.loss / double(parm.df) * arma::diagvec(arma::inv(arma::trans(X) * X)));
    }

    RegParm stepWiseF(RegFun* fitMethod, const mat& X, const vec& Y,
        const uvec& holdIndex, const uvec& fitIndex, const int& K, arma::uvec& bestInd)
    {
        RegParm parm;

        const auto& p = X.n_cols;
        const auto& hp = holdIndex.n_elem;

        // assert(X.n_rows == Y.n_elem);
        // assert(p >= hp);

        bestInd.resize(hp + K);
        arma::uword minInd = p + 1;
        bestInd.fill(p + 1);

        arma::uword nhp = 0;
        for (auto k = 0; k < hp; k ++) {
            if (arma::is_finite(X.col(holdIndex[k]))) {
                bestInd[nhp] = holdIndex[k];
                nhp++;
            }
        }

        RegParm parmT = fitMethod->fit(X.cols(holdIndex), Y);
        double minBIC = parmT.BIC;

        arma::uword findp = 0;
        // Rcpp::Rcout << "p + 1 " << p + 1 <<  "BestInd:" << bestInd << '\n';
        // Rcpp::Rcout << "nhp:" << nhp << ", K: " << K << '\n';
        for (arma::uword i = nhp; i < nhp + K; i++) {
            // Rcpp::Rcout << "i: " << i << '\n';
            for (arma::uword j = 0; j < p; j++) {

                if (arma::any(bestInd == j))
                    continue;
                if (arma::all(fitIndex != j))
                    continue;
                if (!arma::is_finite(X.col(j)))
                    continue; // TODO: why add this?

                bestInd[i] = j;
                parmT = fitMethod->fit(X.cols(bestInd.head(i + 1)), Y);
                // Rcpp::Rcout << "(i, j) = (" << i << "," << j << "), BIC: " << parmT.BIC << '\n';
                if (parmT.BIC < minBIC) {
                    minInd = j;
                    minBIC = parmT.BIC;
                }
            }
            bestInd[i] = minInd;
            if (minInd == p + 1) {
                break;
            }
            // Rcpp::Rcout << i << ", BIC " <<minBIC << ", minInd " << minInd <<  '\n';
            minInd = p + 1;
            findp++;
            if (findp == K)
                break;
        }
        // Rcpp::Rcout << "p + 1 " << p + 1 <<  "BestInd:" << bestInd << '\n';
        bestInd = arma::sort(bestInd.rows(arma::find(bestInd < p + 1)));
        // Rcpp::Rcout << "BestInd:" << bestInd << '\n';

        parm = fitMethod->fit(X.cols(bestInd), Y);
        return parm;
    }
} // RegSol namespace
} //  SubGUIDE namespace
