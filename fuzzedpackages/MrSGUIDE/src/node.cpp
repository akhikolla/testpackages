//
//  node.cpp
//  SubGUIDE
//
//  Created by Peigen Zhou on 8/18/18.
//

#include "node.hpp"

namespace SubGuide {

namespace Node {
    using arma::is_finite;

    double NodeModel::evaluateLoss(const mat& X, const mat& Y)
    {
        mat yhat = predict(X);
        this->testLoss = RegSol::getLoss(Y, yhat, this->fitMethod);
        return this->testLoss;
    }

    void stepNodeModel::fit(const mat& X, const mat& Y, const uvec& holdIndex,
        const uvec& fitIndex, const int& bestK,
        const vec& XMean_)
    {
        this->XMean = XMean_;
        this->bestInds.resize(0);
        for (auto i = 0; i < Y.n_cols; i++) {
            RegSol::RegParm parm;
            arma::uvec bestInd;

            parm = stepWiseF(this->fitMethod, X, Y.col(i), holdIndex, fitIndex, bestK, bestInd);
            this->bestInds.push_back(bestInd);

            const mat& tmpX = X.cols(this->bestInds.at(i));

            parm.betaSE = this->fitMethod->calStdErr(tmpX, Y.col(i));

            this->parms.push_back(parm);
            this->trainLoss += parm.loss;
        }
    }

    std::vector<arma::vec> stepNodeModel::fitBeta(const mat& X, const mat& Y,
        const uvec& holdIndex,
        const uvec& fitIndex,
        const int& bestK,
        const int& tailn)
    {
        std::vector<arma::vec> result;

        for (auto i = 0; i < Y.n_cols; i++) {
            RegSol::RegParm parm;
            mat tmp = X.cols(this->bestInds.at(i));
            tmp = imputeValue(tmp, this->XMean(this->bestInds.at(i)));

            parm = this->fitMethod->fit(tmp, Y.col(i)); // original index to fit.

            // parm = stepWiseF(this->fitMethod, X, Y.col(i), holdIndex, fitIndex, bestK, bestInd); // step wise to fit
            result.push_back(parm.beta.tail(tailn));
        }
        return result;
    }

    mat stepNodeModel::predict(const mat& X)
    {
        const uword& N = X.n_rows;
        auto yp = this->bestInds.size();

        mat result(N, yp);
        for (auto i = 0; i < yp; i++) {
            const mat& tmpX = X.cols(this->bestInds.at(i));
            if (is_finite(tmpX)) {
                result.col(i) = this->fitMethod->predict(tmpX, this->parms.at(i));
            } else {
                result.col(i) = this->fitMethod->predict(
                    imputeValue(tmpX, this->XMean(this->bestInds.at(i))),
                    this->parms.at(i));
            }
        }
        return result;
    }

    double getTermSize(node* leaf, const double& alpha)
    {
        if (leaf->terminal || leaf->CriticalAlpha < alpha) {
            return 1.0;
        } else {
            return getTermSize(leaf->left, alpha) + getTermSize(leaf->right, alpha);
        }
    }

    double getInterSize(node* leaf)
    {
        if (leaf->terminal) {
            return 0.0;
        } else {
            return 1.0 + getInterSize(leaf->left) + getInterSize(leaf->right);
        }
    }

    double getTermLoss(node* leaf)
    {
        if (leaf->terminal) {
            return leaf->trainLoss;
        } else {
            return getTermLoss(leaf->left) + getTermLoss(leaf->right);
        }
    }

    double getTermTestLoss(node* leaf)
    {
        if (leaf->terminal) {
            // assert(leaf->testLoss > -1.0);
            return leaf->testLoss;
        } else {
            return getTermTestLoss(leaf->left) + getTermTestLoss(leaf->right);
        }
    }

    void updateNodeMap(node* leaf, std::vector<node*>& result)
    {
        if (leaf) {
            if (leaf->terminal) {
                result.push_back(leaf);
                return;
            } else {
                result.push_back(leaf);
                updateNodeMap(leaf->left, result);
                updateNodeMap(leaf->right, result);
            }
        }
    }

    void destroy(node* leaf)
    {
        if (leaf != nullptr) {
            destroy(leaf->left);
            destroy(leaf->right);
            delete leaf;
        }
    }

    void pruneAlpha(node* leaf, const double& alpha)
    {
        if (leaf->terminal)
            return;
        if (leaf->CriticalAlpha < alpha) {
            leaf->terminal = true;
            destroy(leaf->left);
            destroy(leaf->right);
            leaf->left = nullptr;
            leaf->right = nullptr;
        } else {
            pruneAlpha(leaf->left, alpha);
            pruneAlpha(leaf->right, alpha);
        }
    }

    void evaluateTrain(node* leaf)
    {
        if (leaf->terminal) return;

        leaf->termNodeSize = getTermSize(leaf);
        leaf->termTrainLoss = getTermLoss(leaf);
        leaf->CriticalAlpha = (leaf->trainLoss - leaf->termTrainLoss) / (leaf->termNodeSize - 1.0);

        evaluateTrain(leaf->left);
        evaluateTrain(leaf->right);
    }

    uword predLR(node* leaf, const arma::rowvec& numX, const arma::irowvec& catX)
    {
        uword direction; // 1: left, 2: right

        if (leaf->SplitRole == 'n') {
            const uword& ind = leaf->SplitID;
            const double& val = numX(ind);
            if (leaf->misDirection == 'A') {
                direction = arma::is_finite(val) ? 2 : 1;
            } else {
                if (arma::is_finite(val)) {
                    direction = (val <= leaf->threshold) ? 1 : 2;
                } else {
                    direction = leaf->misDirection == 'L' ? 1 : 2;
                }
            }
        } else {
            const arma::uword& ind = leaf->SplitID - numX.n_cols;
            direction = arma::any(leaf->threshSet == catX(ind)) ? 1 : 2;
        }
        return direction;
    }

    uvec predictLR(SubGuide::Node::node* leaf, const arma::mat& numX,
        const arma::imat& catX)
    {
        // assert(catX.n_rows == numX.n_rows);
        const uword& N = catX.n_rows;
        uvec directions(N);
        for (auto i = 0; i < N; i++)
            directions(i) = predLR(leaf, numX.row(i), catX.row(i));
        return directions;
    }

    void display(node* leaf, const int& indent)
    {
        if (leaf) {
            if (indent) {
                cout << std::setw(indent) << ' ';
            }
            cout << "Node ID: " << leaf->NodeID << ", Critical Alpha: " << leaf->CriticalAlpha << ", ";
            if (leaf->terminal) {
                cout << "[Terminal] n = " << leaf->N << "\n";
                if (indent) {
                    cout << std::setw(indent) << ' ';
                }
                cout << "Parm: \n";
                for (auto& item : leaf->nodeModel->parms) {
                    if (indent) {
                        cout << std::setw(indent) << ' ';
                    }
                    cout << item.beta.t();
                }
            } else {
                // cout << "Critical Alpha: " << leaf->CriticalAlpha;
                cout << "Split Var: " << leaf->SplitID
                     << ", Role: " << leaf->SplitRole << ", ";
                if (leaf->SplitRole == 'n') {
                    cout << "Threshold: " << leaf->threshold
                         << ", Miss Dir: " << leaf->misDirection;
                } else {
                    cout << "ThreshSet: { ";
                    for (auto& item : leaf->threshSet)
                        cout << item << " ";
                    cout << "}";
                }
                cout << std::endl;
            }
            if (leaf->left) {
                display(leaf->left, indent + 4);
            }

            if (leaf->right) {
                display(leaf->right, indent + 4);
            }
        }
    }

    template <class ArmaT>
    void writeArma(std::ostringstream& out, const ArmaT& vec)
    {
        out << "[";
        for (auto i = 0; i < vec.n_elem; i++) {
            out << std::fixed << std::setprecision(6) << vec[i];
            if (i < vec.n_elem - 1) {
                out << ", ";
            }
        }
        out << "]";
    }

    void writeArma(std::ostringstream& out, const uvec& vec,
        const std::vector<std::string>& varName)
    {
        out << "[";
        for (auto i = 0; i < vec.n_elem; i++) {
            out << varName[vec[i]];
            if (i < vec.n_elem - 1) {
                out << ", ";
            }
        }
        out << "]";
    }

    void writeVec(std::ostringstream& out, const std::vector<uvec>& bestInd,
        const int& indent)
    {
        const int& N = bestInd.size();
        for (auto i = 0; i < N; i++) {
            out << "  - ";
            writeArma(out, bestInd[i]);
            if (i < N - 1) {
                out << "\n";
                if (indent) {
                    out << std::setw(indent) << ' ';
                }
            }
        }
    }

    void writeVec(std::ostringstream& out, const std::vector<uvec>& bestInd,
        const int& indent, const std::vector<std::string>& varName)
    {
        const int& N = bestInd.size();
        for (auto i = 0; i < N; i++) {
            out << "  - ";
            if (!varName.empty()) {
                writeArma(out, bestInd[i], varName);
            } else {
                writeArma(out, bestInd[i]);
            }

            if (i < N - 1) {
                out << "\n";
                if (indent) {
                    out << std::setw(indent) << ' ';
                }
            }
        }
    }

    std::string writeNode(node* leaf, const int& indent,
        const std::vector<std::string>& varName)
    {
        std::ostringstream treeOut;
        if (leaf) {
            if (indent) {
                treeOut << std::setw(indent) << ' ';
            }
            treeOut << "ID: " << leaf->NodeID << "\n";

            if (indent) {
                treeOut << std::setw(indent) << ' ';
            }

            treeOut << "Type: ";

            if (leaf->terminal) {
                treeOut << "Terminal\n";

                if (indent) {
                    treeOut << std::setw(indent) << ' ';
                }

                treeOut << "FitIndex: \n";

                if (indent) {
                    treeOut << std::setw(indent) << ' ';
                }

                writeVec(treeOut, leaf->fitInds, indent, varName);
                treeOut << "\n";

                if (indent) {
                    treeOut << std::setw(indent) << ' ';
                }

                treeOut << "Parms: \n";

                for (auto& item : leaf->nodeModel->parms) {
                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }
                    treeOut << "  - ";
                    writeArma(treeOut, item.beta);
                    treeOut << "\n";
                }

                if (indent) {
                    treeOut << std::setw(indent) << ' ';
                }

                treeOut << "Trts: \n";

                for (auto& item : leaf->trtBeta) {
                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }
                    treeOut << "  - ";
                    writeArma(treeOut, item);
                    treeOut << "\n";
                }

                if (indent) {
                    treeOut << std::setw(indent) << ' ';
                }

                treeOut << "SEs: \n";

                for (auto& item : leaf->trtSE) {
                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }
                    treeOut << "  - ";
                    writeArma(treeOut, item);
                    treeOut << "\n";
                }

                return treeOut.str();

            } else {
                treeOut << "Intermediate\n";

                if (indent) {
                    treeOut << std::setw(indent) << ' ';
                }

                if (!varName.empty()) {
                    treeOut << "SplitVar: " << varName[leaf->SplitID] << "\n";
                } else {
                    treeOut << "SplitVar: " << leaf->SplitID << "\n";
                }

                if (leaf->SplitRole == 'n') {

                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }

                    treeOut << "Role: 'num' \n";

                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }

                    treeOut << "Threshold: " << leaf->threshold << "\n";

                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }

                    treeOut << "MisDirection: " << leaf->misDirection << "\n";
                } else {
                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }

                    treeOut << "Role: 'char' \n";


                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }

                    treeOut << "ThreshSet: ";

                    writeArma(treeOut, leaf->threshSet);

                    treeOut << "\n";
                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }

                    treeOut << "MisDirection: " << leaf->misDirection << "\n";
                }

                if (leaf->left) {
                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }
                    treeOut << "Left: \n"
                            << writeNode(leaf->left, indent + 2, varName);
                }
                if (leaf->right) {
                    if (indent) {
                        treeOut << std::setw(indent) << ' ';
                    }
                    treeOut << "Right: \n"
                            << writeNode(leaf->right, indent + 2, varName);
                }
            }
        }
        return treeOut.str();
    }

} // namespace Node
} // namespace SubGuide
