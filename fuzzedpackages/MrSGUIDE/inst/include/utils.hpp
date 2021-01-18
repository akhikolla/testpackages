//
// Created by Peigen Zhou on 7/14/18.
//

#ifndef SUBGUIDE_UTILS_HPP
#define SUBGUIDE_UTILS_HPP

#include "common.hpp"

namespace SubGuide {
    
    arma::vec quantile(const arma::vec& X, const int & part);
    /**
     * QuartileX used for numeric feature, default we discritzie the feature into 4
     * category based on
     * @param feature
     * @param parts
     * @return vector, each integer refer to the group indicator
     */
    arma::ivec quartileX(const arma::vec &x, int parts = 4);
    
    /**
     * Assume feature only start with 1
     * @param feature categorical feature.
     * @param intercept intercept
     * @return Design matrix for categorical variable
     */
    arma::mat hotCoding(const arma::ivec &cx, const arma::ivec &levels,
                        bool reference = true);
    
    arma::mat hotCoding(const arma::ivec &cx, bool reference = true);
    
    /**
     * Design matrix for interaction.
     * @param x1 design matrix for x1
     * @param x2 design matrix for x2
     * @return design matrix for interaction term of x1 and x2
     */
    arma::mat designInt(const arma::mat &x1, const arma::mat &x2);
    
    /**
     * Get all levels combination of a categorical variable
     * x = {1, 2, 3}, return a matrix, {{1, 0, 0} {0, 1, 0} {0, 0, 1}}
     * @param cx categorical variable
     * @return matrix where each row is 0, 1 indicate whether the current level
     * should be include.
     */
    arma::umat getLevels(const arma::ivec &cx);
    
    /**
     * match, find cx %in% x as R, return 1 means in x, 0 mean not in x
     * @param cx
     * @param x
     * @return
     */
    arma::uvec match(const arma::ivec &cx, const arma::ivec &x);
    
    arma::mat imputeMean(const arma::mat &X);
    arma::mat imputeValue(const arma::mat &X, const arma::vec &Xmean);
    arma::vec colMean(const arma::mat &X);
    arma::mat transVec(std::vector<arma::vec> tmp);

} // namespace SubGuide

#endif // SUBGUIDE_UTILS_HPP
