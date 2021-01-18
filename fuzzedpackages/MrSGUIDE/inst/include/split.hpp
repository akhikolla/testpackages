//
// Created by Peigen Zhou on 7/21/18.
//

#ifndef SUBGUIDE_SPLIT_HPP
#define SUBGUIDE_SPLIT_HPP

#include "common.hpp"
#include "regression.hpp"
#include "utils.hpp"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>

namespace SubGuide {
namespace SplitSol {

    using arma::icolvec;
    using arma::imat;
    using arma::ivec;
    using arma::join_rows;
    using arma::mat;
    using arma::ones;
    using arma::umat;
    using arma::uvec;
    using arma::uword;
    using arma::vec;

    /**
 * Tree split function class
 * Split is a base class which use to store all the data info
 * It should be a general class, not design for treatment
 * Todo delete numX, cateX, Y, Trt inside class.
 */
    class Split {
    protected:
        // std::shared_ptr<spdlog::logger> logger;
        // Setting
        int batchNum = 1;
        int minData = 10;

        // Result
        uword varID = {};
        char role = 'n';
        char misDirection = 'n';
        double threshold = 0.0;
        ivec threshSet = {};
        uvec optLeft,
            optRight; // TODO: optLeft, optRight should change in buildTree
            // procedure.

        // Data Info
        uword N, yp, np, cp;
        // std::vector<uvec> bestInd;

        double loss = 0;

        void dataCheck(const mat& nx_, const imat& cx_, const mat& ys_);

    public:
        Split();
        Split(const int& batchNum_, const int& minData_);

        void clear(); // clean previous result

        uword getVarID() const { return this->varID; }
        char getVarRole() const { return this->role; }
        double getThreshold() const { return this->threshold; }
        char getMissDir() const { return this->misDirection; }
        ivec getThreshSet() const { return this->threshSet; }
        double getLoss() const { return this->loss; }

        // dir in {'L', 'R'}
        uvec getSplitVec(const char& dir) const
        {
            return dir == 'L' ? this->optLeft : this->optRight;
        }
        vec chiN, chiC;
    };

    /**
 * Gi Split. Assume nx potentially has inf value. Gi split is design for
 * treatment ys must not has inf value. Must be imputed before pass to
 * GiSplit
 * TODO: It looks wired. Since how to use GiSplit directly?
 * Write a new function to expose?
 */
    class GiSplit : public Split {
    public:
        GiSplit(RegSol::RegFun* fitMethod, const int& batchNum_, const int& minData_, 
            const int& minTrt_, const ivec& trtLevel_);

        void findSplit(const mat& numX, const imat& catX, const mat& Y,
            const icolvec& trt, const std::vector<uvec>& bestInd,
            const uvec& SplitIndex);

        void findThresh(const mat& numX, const imat& catX, const mat& Y, const ivec trt,
            const mat& comX, const uvec& fixIndex,
            const uvec& fitIndex, const int& bestK,
            const std::vector<uvec>& bestInd_, const bool& faster);

        ivec findCateThresh(const ivec& x, const ivec trt, const mat& ys,
            const mat& comX, const uvec& fixIndex,
            const uvec& fitIndex, const int& bestK,
            const std::vector<uvec>& bestInd_, const bool& faster);

        double findNumThresh(const vec& x, const ivec trt, const mat& ys,
            const mat& comX, const uvec& fixIndex,
            const uvec& fitIndex, const int& bestK,
            const std::vector<uvec>& bestInd_, const bool& faster);

    private:
        int minTrt = 5; // TODO: change minTrt Only in GiSplit
        RegSol::RegFun* nodeFitMethod;

        ivec trtLevel; // Trt Level
        // uword tp;      // trt level length

        // Derived data
        mat trtDes; // Trt Design matrix should I hold this?

        // methods
        double lackOfFit(const mat& Main, const mat& cx, const vec& Y_,
            const mat& trtDes);
        void dataChecking(const mat& nx_, const imat& cx_, const mat& ys_,
            const ivec trt_);

        inline bool checkNodeData(const uvec& index, const ivec& trt,
            const ivec& trtLevel, const int& minData,
            const int& minTrt);
        // uvec numStepIndex(const mat &comX, const uvec &offInd, int K, const vec
        // &Y_);
    };

    class RegSplit : Split {
    };
    class ClaSplit : Split {
    };
} // namespace SplitSol
} // namespace SubGuide

#endif // SUBGUIDE_SPLIT_HPP
