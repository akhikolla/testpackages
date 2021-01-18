// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "regression.hpp"
#include "split.hpp"
#include "gitree.hpp"
#include <vector>
#include <fstream>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]



//' Character vector to Integer vector use Rcpp
//'
//' @param x original character vector
//' @param levels the unique value form x, where the index is the output integer vector value
//' @return integer vector
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector characterToInteger(const Rcpp::CharacterVector &x, const Rcpp::CharacterVector &levels) {

    const auto &N = x.length();
    const auto &p = levels.length();

    Rcpp::IntegerVector result(N);

    for (int i = 0; i < p; i++) {
        for(int j = 0; j < N; j ++) {
            if (Rcpp::CharacterVector::is_na(x(j))) result(j) = p + 1;
            if (x(j) == levels(i)) result(j) = i + 1;
        }
    }
    return result;
}

//' R dataframe to numeric vector
//'
//' Function will change the original numerical dataframe to a numerical matrix, where missing value will change to Inf.
//'
//' @param numX numerical dataframe
//' @return numerical matrix
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix dataFramToNumeric(const Rcpp::DataFrame &numX) {
    const auto &p = numX.ncol();
    const auto &N = numX.nrow();

    Rcpp::NumericMatrix result(N, p);

    for (int i = 0; i < p; i++) {
        const Rcpp::NumericVector &x = numX[i];
        for (int j = 0; j < N; j++) {
            result(j, i) = Rcpp::NumericVector::is_na(x(j)) ? R_PosInf : x(j);
        }
    }
    return result;
}

//' Character or factor dataframe change to integer matrix.
//'
//' Character dataframe, for each column, use \code{characterToInteger()} with unique value from character vector
//'
//' @param charX character dataframe or factor dataframe
//' @param levels, the desired levels for each character vector
//' @return list intX the integer matrix x, with corresponding levels.
//' @noRd
// [[Rcpp::export]]
Rcpp::List characterDict(const Rcpp::DataFrame &charX, const Rcpp::List &levels) {
    const auto &p = charX.ncol();
    const auto &N = charX.nrow();

    Rcpp::IntegerMatrix cX(N, p);

    for (int i = 0; i < p; i ++) {
        const Rcpp::CharacterVector &x = charX[i];
        cX.column(i) = characterToInteger(x, levels[i]);
    }

    return Rcpp::List::create(Rcpp::Named("intX") = cX,
                        Rcpp::Named("Levels") = levels);
}



/*
//' Initial logger use spdlog library.
//'
// [[Rcpp::export]]
void initLog(const int level) {
     auto logger = spdlog::stdout_logger_st("tree");
     auto loggerS = spdlog::stdout_logger_st("split");
     auto loggerD = spdlog::stdout_logger_st("debug");

     if (level == 1) {
         loggerS->set_level(spdlog::level::off);
         logger->set_level(spdlog::level::off);
         loggerD->set_level(spdlog::level::off);
     } else if (level == 2) {
         loggerS->set_level(spdlog::level::info);
         logger->set_level(spdlog::level::info);
         loggerD->set_level(spdlog::level::info);
     } else {
         loggerS->set_level(spdlog::level::debug);
         logger->set_level(spdlog::level::debug);
         loggerD->set_level(spdlog::level::debug);
     }
}
*/

//' MrSGUIDE C++ function
//'
//' @title MrSGUIDE stepwise
//'
//' @author Peigen Zhou
//'
//'
//' @param numX numerical X matrix
//' @param catX categorical X matrix
//' @param Y outcome Y matrix
//' @param trt treatment vector
//' @param splitIndex variable used for split
//' @param fitIndex variables can be used for fit
//' @param holdIndex variable must include in the fitting model
//' @param bestK maximal number of variables used in the outcome model for prognostic control
//' @param maxDepth maximal depth
//' @param minData minimum sample in each node
//' @param minTrt minimum treatment and placebo sample in each node
//' @param batchNum related with exhaustive search for numerical split variable
//' @param CVFold cross validation times
//' @param CVSE cross validation SE
//' @param bootNum bootstrap number
//' @param alpha desire alpha levels for confidence interval with respect to treatment parameters
//' @param faster related with tree split searching
//' @param display Whether display tree in the end
//' @param varName variable names
//' @param treeName yaml file for save the tree
//' @param nodeName file same for each node
//' @param bootName file save bootstrap calibrate alpha
//' @param impName important variable file name
//' @noRd
//'
// [[Rcpp::export]]
void GiStepWisePure(const arma::mat &numX, const arma::imat &catX,
                const arma::mat &Y, const arma::ivec trt,
                const arma::uvec &splitIndex, const arma::uvec &fitIndex,
                const arma::uvec &holdIndex, const int &bestK,
                const int &maxDepth, const int& minData,
                const int &minTrt, const int &batchNum,
                const int &CVFold, const double &CVSE,
                const int &bootNum, const double &alpha,
                const bool &faster, const bool &display,
                const std::vector<std::string>&varName,
                const std::string &treeName, const std::string &nodeName,
                const std::string &bootName, const std::string &impName) {
    SubGuide::RegSol::RegFun *linearNode;
    linearNode = new SubGuide::RegSol::LinReg;

    SubGuide::Tree::GiTree resTree(linearNode, maxDepth, minData, minTrt, batchNum, CVFold, CVSE, bootNum, alpha);

    resTree.fit(numX, catX, Y, trt, splitIndex, fitIndex, holdIndex, bestK, faster);

    if (display) resTree.display();


    std::ofstream myfile;

    myfile.open(impName);
    if (resTree.importanceScoreN.n_rows > 0)    myfile << resTree.importanceScoreN;
    if (resTree.importanceScoreC.n_rows > 0)    myfile << resTree.importanceScoreC;
    myfile.close();

    myfile.open(treeName);
    myfile << resTree.writeTree(varName);
    myfile.close();

    myfile.open(nodeName);
    myfile << "node\n";
    myfile << resTree.predictNode(numX, catX);
    myfile.close();

    if (bootNum > 10) {
        myfile.open(bootName);
        myfile << resTree.getBootAlpha();
        myfile.close();
    }
}

/**
 * Export Variable selection unit
 */

/*
// [[Rcpp::export]]
Rcpp::List GiNumThresh(const arma::vec &x, const arma::mat &numX, const arma::mat &Y, const arma::ivec trt,
                       const std::vector<arma::uvec> &bestInd, const int &batchNum, const int& minData,
const int &minTrt) {

SubGuide::RegSol::RegFun *linearNode;
linearNode = new SubGuide::RegSol::LinReg;

SubGuide::SplitSol::GiSplit splitRes(linearNode, batchNum, minData, minTrt);
splitRes.findNumThresh(x, numX, trt, Y, bestInd);

Rcpp::List result = Rcpp::List::create(Rcpp::Named("Threshold") = splitRes.getThreshold(),
                                       Rcpp::Named("MisDirection") = splitRes.getMissDir());
return result;
}

// [[Rcpp::export]]
Rcpp::List GiCateThresh(const arma::ivec &x, const arma::mat &numX, const arma::mat &Y, const arma::ivec trt,
                        const std::vector<arma::uvec> &bestInd, const int &batchNum, const int& minData,
const int &minTrt) {

SubGuide::RegSol::RegFun *linearNode;
linearNode = new SubGuide::RegSol::LinReg;

SubGuide::SplitSol::GiSplit splitRes(linearNode, batchNum, minData, minTrt);
splitRes.findCateThresh(x, numX, trt, Y, bestInd);

Rcpp::List result = Rcpp::List::create(Rcpp::Named("ThreshSet") = splitRes.getThreshSet());
return result;
}
*/
