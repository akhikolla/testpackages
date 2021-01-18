//
//  node.hpp
//  SubGUIDE
//
//  Created by Peigen Zhou on 8/18/18.
//

#ifndef node_hpp
#define node_hpp

#include "common.hpp"
#include "regression.hpp"
#include <boost/math/distributions/normal.hpp>
#include <iomanip>
#include <numeric>
#include <sstream>

namespace SubGuide {

    namespace Node {
        using arma::ivec;
        using arma::mat;
        using arma::umat;
        using arma::uvec;
        using arma::uword;
        using arma::vec;

        // TODO: design logic is not clear
        /*
         struct NodeModel {
         double trainLoss;
         double testLoss;
         RegSol::RegFun *fitMethod;

         // results
         std::vector<uvec> bestInds;
         std::vector<RegSol::RegParm> parms;

         // helper Data
         vec XMean; // use to impute X columns.

         NodeModel(RegSol::RegFun *fitMethod_) : fitMethod(fitMethod_){};
         // NodeModel(RegSol::RegFun* fitMethod_, const mat& X, const mat& Y, const
         // uvec& fixIndex, const uvec& fitIndex, const int& bestK, const vec&
         XMean_);

         void fit(const mat &X, const mat &Y, const uvec &fixIndex,
         const uvec &fitIndex, const int &bestK, const vec &XMean_);
         double evaluateLoss(const mat &X, const mat &Y);
         mat predict(const mat &X);
         };
         */

        class NodeModel {
        protected:
            RegSol::RegFun *fitMethod;
            vec XMean;

        public:
            double trainLoss;
            double testLoss;
            std::vector<RegSol::RegParm> parms;
            NodeModel(RegSol::RegFun *fitMethod_)
            : fitMethod(fitMethod_), trainLoss(-0.1), testLoss(-0.1){};
            virtual mat predict(const mat &X) = 0;
            double evaluateLoss(const mat &X, const mat &Y);
        };

        class stepNodeModel : public NodeModel {
        public:
            std::vector<uvec> bestInds;
            stepNodeModel(RegSol::RegFun *fitMethod_) : NodeModel(fitMethod_){};

            void fit(const mat &X, const mat &Y, const uvec &holdIndex,
                     const uvec &fitIndex, const int &betsK, const vec &XMean_);
            std::vector<arma::vec> fitBeta(const mat &X, const mat &Y,
                                           const uvec &holdIndex, const uvec &fitIndex,
                                           const int &bestK, const int &tailn);
            mat predict(const mat &X);

            // results
        };

        // TODO: Add GiNode
        struct node {
            size_t NodeID;
            bool terminal = true;

            uword SplitID;
            char SplitRole;

            double threshold;
            ivec threshSet;
            char misDirection;

            node *left;
            node *right;
            stepNodeModel *nodeModel;

            std::vector<uvec> fitInds; // Only shows numx fitted variable
            std::vector<vec> trtBeta;
            std::vector<vec> trtSE;

            node(const size_t &id)
            : NodeID(id), threshold(arma::datum::inf), threshSet({}), left(nullptr),
            right(nullptr), nodeModel(nullptr), trainLoss(-1.0), testLoss(-1.0),
            termNodeSize(0.0), termTrainLoss(-1.0), CriticalAlpha(0.0) {}

            // cvInfo
            int N;            // Node Size
            double trainLoss; // mean sse
            double testLoss = 0.0;  // mean sse
            double termNodeSize;
            double termTrainLoss;
            double CriticalAlpha = 0.0;
            double SplitChi = 0.0;
        };

        // Methods
        double getTermSize(node *leaf, const double &alpha);
        double getInterSize(node *leaf);
        inline double getTermSize(node *leaf) { return getTermSize(leaf, -1); };

        double getTermLoss(node *leaf);
        double getTermTestLoss(node *leaf);

        void updateNodeMap(node *leaf, std::vector<node *> &result);
        void destroy(node *leaf);
        void pruneAlpha(node *leaf, const double &alpha);

        // given single observation numX and catX, within current node, predict left or
        // Right
        uword predLR(node *leaf, const arma::rowvec &numX, const arma::irowvec &catX);

        // given data matrix, predict left or Right
        uvec predictLR(node *leaf, const arma::mat &numX, const arma::imat &catX);

        void display(node *leaf, const int &indent);
        std::string writeNode(node *leaf, const int &indent,
                              const std::vector<std::string> &varName);

        void evaluateTrain(node* leaf);

    } // namespace Node
} // namespace SubGuide

#endif /* node_hpp */
