#ifndef SUBGUIDE_GITREE_HPP
#define SUBGUIDE_GITREE_HPP

#include "tree.hpp"

namespace SubGuide {
namespace Tree {
    class GiTree : public tree {
    public:
        GiTree(RegSol::RegFun* fitMethod, const int& maxDepth_, const int& minData,
            const int& minTrt, const int& batchNum_, const int& CVFold_,
            const double& CVSE_, const int& bootNum, const double& alpha);

        void fit(const mat& numX, const imat& catX, const mat& Y, const ivec& Trt,
            const uvec& splitIndex_, const uvec& fitIndex_,
            const uvec holdIndex_, const int& K_, const bool faster_);
        mat getBootAlpha() const {return this->BootAlpha;}

    private:
        // Methods
        RegSol::RegFun* fitMethod; // Node model
        SplitSol::GiSplit* splitMethod; // Split method
        // Treatment
        ivec trtLevel;
        uword tp;
        int minTrt;

        // Settings
        bool faster = false;

        // inference
        int bootNum;
        double alphaLevel;
        mat BootAlpha;

        node* buildTree(const mat& numX, const imat& catX, const mat& Y,
            const ivec& Trt, const size_t& id, const int& depth_);

        void crossValidation(node* cvroot, const mat& numX, const imat& catX,
            const mat& Y, const ivec& Trt, const int& CVFold,
            const double& CVSE);

        mat boostrapCI(const mat& numX, const imat& catX, const mat& Y,
            const ivec& Trt, const int& bootNum, const double& alpha);

        mat predictY(node* root, const mat& numX, const imat& catX, const ivec& trt);
    };

    inline std::vector<uvec> refineFit(const std::vector<uvec>& x,
        const int& np);

    inline std::vector<vec> refineTrt(const std::vector<RegSol::RegParm>& parms,
        const int& tp, const bool& beta);

    std::vector<mat> GetCoverMat(const std::vector<vec>& TrueTrt,
        const std::vector<vec>& center,
        const std::vector<vec>& se, const vec& Alpha);

    mat GetBootAlpha(const vec& gamma, const vec& theta, const double& alphaL,
        const vec& alphaK);
    
    inline arma::mat createFitMatrix(const arma::mat& numX, const arma::ivec& Trt, const arma::ivec &trtlevel);
    
}
}

#endif // SUBGUIDE_TREE_HPP
