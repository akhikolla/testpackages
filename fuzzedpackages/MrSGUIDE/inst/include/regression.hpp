//
// Created by Peigen Zhou on 7/15/18.
//

#ifndef SUBGUIDE_REGRESSION_HPP
#define SUBGUIDE_REGRESSION_HPP

#include "common.hpp"
#include "utils.hpp"

namespace SubGuide {
    namespace RegSol {
        using arma::colvec;
        using arma::ivec;
        using arma::mat;
        using arma::uvec;
        using arma::vec;

        struct RegParm {
            int n;
            int p;
            int df;
            vec beta;
            vec betaSE;
            
            double BIC;
            double loss; // sum of square error
            
            RegParm() = default;
            RegParm(const int &n_, const int &p_)
            : n(n_), p(p_), df(n - p), beta(p, arma::fill::zeros),
            betaSE(p, arma::fill::zeros), BIC(arma::datum::inf),
            loss(arma::datum::inf) {}
            void update(const int &n_, const int &p_);
            void print() {
                cout << "Parm: \n n: " << n << " p: " << p << " df: " << df
                << " loss: " << loss << "\n";
            }
        };
        
        class RegFun {
        protected: // changeable by inherited class but not outside
            RegParm parm;
            uvec offsetIndex;
            
        public:
            RegFun() : parm(0, 0){};
            
            virtual RegParm fit(const mat &X_, const vec &Y_) = 0;
            virtual vec predict(const mat &X) = 0;
            virtual vec calStdErr(const mat &X, const vec &Y) = 0;
            virtual vec
            predict(const mat &X,
                    const RegParm &parm) = 0; // TODO: check whether need this function
            
            virtual double getLoss(const vec &Y, const vec &Ypred) = 0;
            // Data access functions
            RegParm getParm() { return parm; }
        };
        
        class LinReg : public RegFun {
        public:
            LinReg() = default;
            
            RegParm fit(const mat &X, const vec &Y) override;
            vec predict(const mat &X) override { return predict(X, this->parm.beta); };
            
            vec predict(const mat &X, const RegParm &parm) override {
                return predict(X, parm.beta);
            }
            double getLoss(const vec& Y, const vec& Ypred) override;
            vec calStdErr(const mat &X, const vec &Y) override;
            // Data access functions
        private:
            vec predict(const mat &X, const vec &beta);
        };
        
        // TODO: Finish fit function. Use MLE to estimate beta. I think need to use
        // greadient descent.
        class LogReg : public RegFun {
        public:
            LogReg() = default;
            
            RegParm fit(const mat &X, const vec &Y) override;
            vec predict(const mat &X) override { return predict(X, this->parm.beta); };
            vec predict(const mat &X, const RegParm &parm) override {
                return predict(X, parm.beta);
            }
            double getLoss(const vec& Y, const vec& Ypred) override;
            vec calStdErr(const mat &X, const vec &Y) override;
            
        private:
            vec predict(const mat &X, const vec &beta);
        };
        
        RegParm stepWiseF(RegFun *fitMethod, const mat &X, const vec &Y,
                          const uvec &holdIndex, const uvec &fitIndex, const int &K,
                          arma::uvec &bestInd);
        inline double getLoss(const mat& y, const mat& yhat, RegFun *fitMethod)   {
            double testLoss = 0.0;
            for (auto i = 0; i < yhat.n_cols; i++)
                testLoss += fitMethod->getLoss(y.col(i), yhat.col(i));
            return testLoss;
        }
    } // namespace RegSol
} // namespace SubGuide

#endif // SUBGUIDE_REGRESSION_HPP
