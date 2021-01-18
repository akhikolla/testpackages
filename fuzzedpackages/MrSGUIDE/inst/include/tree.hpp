//
// Created by Peigen Zhou on 8/13/18.
//

#ifndef SUBGUIDE_TREE_HPP
#define SUBGUIDE_TREE_HPP

#include <algorithm>
#include "common.hpp"
#include "node.hpp"
#include "regression.hpp"
#include "split.hpp"

namespace SubGuide {
namespace Tree {

using arma::imat;
using arma::irowvec;
using arma::ivec;
using arma::mat;
using arma::rowvec;
using arma::umat;
using arma::uvec;
using arma::uword;
using arma::vec;
using Node::node;

/*
 * Tree should only contains operations for tree,
 * Building operation can be implement inside each derived class
 */
class tree {
   protected:
    // std::shared_ptr<spdlog::logger> logger;
    // Data info
    uword N, yp, np, cp;

    // Settings
    int maxDepth = 3;
    int minData = 10;
    int batchNum = 1;
    int CVFold = 10;
    double CVSE = 0.0;

    uvec splitIndex = {};    //  Numerical Variable used in Split
    uvec fitIndex = {};      // Numerical Variable used in fit pool (step wise)
    uvec holdIndex = {};     // Numerical Variable keep in fit model

    uvec fixIndex = {};    //  For provided comX, which are fixed.
    int bestK;

    // Results
    node* root = nullptr;
    std::vector<node*> mainTreeMap;


    // Methods
    node* predict(node* leaf, const rowvec& numX, const irowvec& catX);
    void evaluateTrain(node* leaf);
    vec getAlphaVector(const std::vector<node*>& nodesMap);
    void destroy(node* leaf);
    uvec predictNode(node* Root, const arma::mat& numX, const arma::imat& catX);

   public:
    tree() = default;
    tree(const int& maxDepth_, const int& minData_, const int& batchNum_,
         const int& CVFold_, const double& CVSE_);

    void display() { Node::display(root, 0); }
    
    std::string writeTree(const std::vector<std::string>& varNames) {
        return Node::writeNode(root, 0, varNames);
    };
    
    node* predict(const rowvec& numX, const irowvec& catX) {
        return predict(root, numX, catX);
    };
    void getImportance() {
        Rcpp::Rcout << importanceScoreN << '\n'
                    << importanceScoreC << std::endl;
    }

    // Get function
    node* getRoot() { return root; }
    uvec predictNode(const mat& numX, const imat& catX) {
        return predictNode(this->root, numX, catX);
    };
    ~tree() = default;

    vec importanceScoreN = {};
    vec importanceScoreC = {};
   private:
    void pruneAlpha(node* leaf, const double& alpha);
};

class RegTree : public tree {
   public:
    RegTree(const int& maxDepth_, const int& minData, const int& minTrt,
            const int& batchNum_, RegSol::RegFun* fitMethod);

    void fit(const mat& numX, const imat& catX, const mat& Y);
    void crossValidation(const mat& numX, const imat& catX, const mat& Y,
                         const int& CVFold, const double& CVSE);

   private:
    // Methods
    RegSol::RegFun* fitMethod;         // Node model
    SplitSol::GiSplit* splitMethod;    // Split method

    // Settings
    node* buildTree(const mat& numX, const imat& catX, const mat& Y,
                    const size_t& id, const int& depth_);
    // void evaluateTest(node* leaf, const mat& numX, const imat& catX, const mat& Y);
    
};

class ClaTree : public tree {
   public:
    ClaTree(const int& maxDepth_, const int& minData, const int& minTrt,
            const int& batchNum_, RegSol::RegFun* fitMethod);

    void fit(const mat& numX, const imat& catX, const mat& Y);
    void crossValidation(const mat& numX, const imat& catX, const mat& Y,
                         const int& CVFold, const double& CVSE);

   private:
    // Methods
    RegSol::RegFun* fitMethod;       // Node model
    SplitSol::Split* splitMethod;    // Split method
    // Settings
    node* buildTree(const mat& numX, const imat& catX, const mat& Y,
                    const size_t& id, const int& depth_);
    // void evaluateTest(node* leaf, const mat& numX, const imat& catX, const mat& Y);
};

}    // namespace Tree
void GiStepWise(const arma::mat& numX, const arma::imat& catX,
                const arma::mat& Y, const arma::ivec trt,
                const arma::uvec& splitIndex, const arma::uvec& fitIndex,
                const arma::uvec& holdIndex, const int& bestK,
                const int& maxDepth, const int& minData, const int& minTrt,
                const int& batchNum, const int& CVFold, const double& CVSE,
                const std::vector<std::string>& varName,
                const std::string& filename);
}    // namespace SubGuide

#endif    // SUBGUIDE_TREE_HPP
