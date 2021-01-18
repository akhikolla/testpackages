//
// Created by Peigen Zhou on 8/13/18.
//

#include "tree.hpp"

namespace SubGuide {
    namespace Tree {
        using arma::is_finite;
        using Node::NodeModel;
        using Node::stepNodeModel;
        
        tree::tree(const int &maxDepth_, const int &minData_, const int &batchNum_,
                   const int &CVFold_, const double &CVSE_)
        : maxDepth(maxDepth_), minData(minData_), batchNum(batchNum_),
        CVFold(CVFold_), CVSE(CVSE_) {
            // logger = spdlog::get("tree");
        };
        
        /**
         * predict observations nodeid
         * Missing Data (only numerical)
         */
        node *tree::predict(node *leaf, const rowvec &numX, const irowvec &catX) {
            // assert(numX.n_elem == np);
            // assert(catX.n_elem == cp);
            node *current = leaf;
            
            while (!current->terminal) {
                current =
                Node::predLR(current, numX, catX) == 1 ? current->left : current->right;
            }
            return current;
        }
        
        uvec tree::predictNode(node *Root, const arma::mat &numX,
                               const arma::imat &catX) {
            const uword &N = numX.n_rows;
            uvec result(N, arma::fill::zeros);
            
            for (uword i = 0; i < N; i++) {
                result(i) = this->predict(Root, numX.row(i), catX.row(i))->NodeID;
            }
            return result;
        }
        
        vec tree::getAlphaVector(const std::vector<node *> &nodesMap) {
            std::vector<double> alpha;
            alpha.push_back(0.0);
            for (auto &item : nodesMap) {
                if (!item->terminal) {
                    alpha.push_back(item->CriticalAlpha - 0.01);
                }
            }
            vec result(alpha);
            return arma::sort(result);
        }
        

        /*
         void GiTree::evaluateTest(node* leaf, const mat& numX, const imat& catX,
         const mat& Y, const ivec& Trt)
         {
         if (Y.n_rows == 0) {
         leaf->testLoss = 0.0;
         return;
         }
         // assert(this->trtLevel.n_elem > 1);
         
         const mat& tmp = hotCoding(Trt, this->trtLevel, true);
         mat comX(tmp.n_rows, this->np + this->tp, arma::fill::ones);
         
         for (uword i = 0; i < this->np + this->tp; i++) {
         if (i < this->np) {
         comX.col(i) = numX.col(i);
         } else {
         if (i == this->np) continue; // jump one;
         comX.col(i) = tmp.col(i - (this->np + 1));
         }
         }
         
         leaf->testLoss = leaf->nodeModel->evaluateLoss(comX, Y) / (double)Y.n_elem;
         
         // assert(leaf->testLoss > 0.0);
         
         if (leaf->terminal)
         return;
         
         uvec directions = Node::predictLR(leaf, numX, catX);
         
         const uvec& left = arma::find(directions == 0);
         const uvec& right = arma::find(directions == 1);
         
         this->evaluateTest(leaf->left, numX.rows(left), catX.rows(left),
         Y.rows(left), Trt.rows(left)); this->evaluateTest(leaf->right, numX.rows(right),
         catX.rows(right), Y.rows(right), Trt.rows(right));
         }
         */
        
    } // namespace Tree
    /*
     void GiStepWise(const arma::mat &numX, const arma::imat &catX,
     const arma::mat &Y, const arma::ivec trt,
     const arma::uvec &splitIndex, const arma::uvec &fitIndex,
     const arma::uvec &holdIndex, const int &bestK,
     const int &maxDepth, const int &minData, const int &minTrt,
     const int &batchNum, const int &CVFold, const double &CVSE,
     const std::vector<std::string> &varName,
     const std::string &filename) {
     RegSol::RegFun *linearNode;
     linearNode = new RegSol::LinReg;
     Tree::GiTree resTree(maxDepth, minData, minTrt, batchNum, linearNode);
     
     resTree.fit(numX, catX, Y, trt, splitIndex, fitIndex, holdIndex, bestK);
     
     if (CVFold > 0) {
     resTree.crossValidation(numX, catX, Y, trt, CVFold, CVSE);
     }
     
     resTree.display();
     
     std::ofstream myfile;
     myfile.open(filename);
     myfile << resTree.writeTree(varName) << "\n";
     myfile.close();
     }
     */
} // namespace SubGuide
