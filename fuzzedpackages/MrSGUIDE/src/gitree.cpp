#include "gitree.hpp"

namespace SubGuide
{
namespace Tree
{

using arma::is_finite;
using Node::NodeModel;
using Node::stepNodeModel;

GiTree::GiTree(RegSol::RegFun *fitMethod_, const int &maxDepth_,
               const int &minData_, const int &minTrt_, const int &batchNum_,
               const int &CVFold_, const double &CVSE_, const int &bootNum_,
               const double &alpha_)
    : tree(maxDepth_, minData_, batchNum_, CVFold_, CVSE_),
      fitMethod(fitMethod_), minTrt(minTrt_), bootNum(bootNum_),
      alphaLevel(alpha_) {}

void GiTree::fit(const mat &numX, const imat &catX, const mat &Y,
                 const ivec &Trt, const uvec &splitIndex_,
                 const uvec &fitIndex_, const uvec holdIndex_, const int &K_,
                 const bool faster_)
{
    this->N = Y.n_rows;
    this->np = numX.n_cols;
    this->cp = catX.n_cols;
    this->yp = Y.n_cols;

    this->splitIndex = splitIndex_;
    this->fitIndex = fitIndex_;
    this->holdIndex = holdIndex_;
    this->bestK = K_;
    this->faster = faster_;

    this->trtLevel = arma::unique(Trt);
    this->splitMethod = new SplitSol::GiSplit(
        fitMethod, this->batchNum, this->minData, this->minTrt, this->trtLevel);

    this->tp = trtLevel.n_elem;
    this->np = numX.n_cols;

    importanceScoreN.resize(np);
    importanceScoreC.resize(cp);
    importanceScoreN.zeros();
    importanceScoreC.zeros();

    const uword &hp = holdIndex.n_elem;

    this->fixIndex.resize(hp + this->tp);
    this->fixIndex.head(hp) = holdIndex;
    for (uword i = 0; i < this->tp; i++)
        this->fixIndex(i + hp) = this->np + i;

    // logger->debug("trtLevel {} \nfixIndex: {}", trtLevel.t(), fixIndex.t());
    // logger->info("Start building tree");

    this->root = buildTree(numX, catX, Y, Trt, 1, 1);
    // logger->info("Finish building");

    // logger->info("Start Cross validation: {} folds, {} SE", CVFold, CVSE);
    this->crossValidation(this->root, numX, catX, Y, Trt, this->CVFold,
                          this->CVSE);

    this->mainTreeMap.clear();
    Node::updateNodeMap(this->root, this->mainTreeMap);

    this->BootAlpha.resize(3, this->yp);
    this->BootAlpha.fill(this->alphaLevel);

    if (this->bootNum > 10)
    {
        //logger->info("Start bootstrap confidence interval\nbootstrap sample: {}, alpha: {}", this->bootNum, this->alphaLevel);
        if (this->root->terminal)
        {
            cout << "Tree only has one node. No need for bootstrap!\n";
            //logger->info("Tree only has one node. No need bootstrap!");
        }
        else
        {
            //logger->info("Use Bootstarp for confidence interval! SLOW!");
            cout << "Use bootstrap for confidence interval! SLOW!\n";
            this->BootAlpha = this->boostrapCI(numX, catX, Y, Trt, this->bootNum, this->alphaLevel);
        }
    }
}

/*
* Missing Data should considerate
* 1. Fit a node model first
* 2. Split variable based on best node model
* 3. Find best threshold or set based on split variables
* TODO: Not safe if use multicore in fitMethod and splitMethod
*/
node *GiTree::buildTree(const mat &numX, const imat &catX, const mat &Y,
                        const ivec &Trt, const size_t &id, const int &depth)
{
    Rcpp::checkUserInterrupt();
    // combine X: numX + 1 + trt2 + .. + trt_tp
    mat comX = createFitMatrix(numX, Trt, this->trtLevel);

    vec Xmean = colMean(comX);

    node *result = new node(id);
    stepNodeModel *nodeFit = new stepNodeModel(fitMethod);
    //Rcpp::Rcout << "NodeId: " << id << '\n';
    //logger->debug("NodeID: {}, bestK: {}", id, this->bestK);
    // Rcpp::Rcout << "fitIndex" << this->fitIndex << "\n";
    nodeFit->fit(imputeValue(comX, Xmean), Y, this->fixIndex, this->fitIndex,
                 this->bestK, Xmean);
    result->nodeModel = nodeFit;
    result->fitInds = refineFit(nodeFit->bestInds, this->np);
    result->trtBeta = refineTrt(nodeFit->parms, this->tp, true); // for beta

    // for (int vi = 0; vi < Y.n_cols; vi++) {
    //   Rcpp::Rcout <<  "vi: " << vi << ", fiinde: " << result->fitInds[vi].t() << '\n'
    //   << result->trtBeta[vi] << '\n';
    // }

    result->trtSE = refineTrt(nodeFit->parms, this->tp, false);  // for SE

    // To make alpha sample size irrelevant

    result->trainLoss = nodeFit->trainLoss;
    result->N = Y.n_rows;

    //logger->debug("Tree Id: {}, depth: {}, loss: {}", id, depth,
    //              result->trainLoss);

    if (depth >= this->maxDepth || Trt.n_elem < 2 * this->minData)
    {
        //    logger->debug("Terminal due to maxDepth and minData");
        return result;
    }
    // cout << "Node ID: " << id << endl;

    //Rcpp::Rcout << "Fitinds: " << result->fitInds[0] << '\n';

    splitMethod->findSplit(numX, catX, Y, Trt, result->fitInds, this->splitIndex);
    // Rcpp::Rcout << "ChiN: " << splitMethod->chiN.t() << '\n';
    // Rcpp::Rcout << "ChiC: " << splitMethod->chiC.t() << '\n';
    // Rcpp::Rcout << "Split at: " << splitMethod->getVarID() << '\n';
    this->importanceScoreN += splitMethod->chiN * Trt.n_elem;
    this->importanceScoreC += splitMethod->chiC * Trt.n_elem;

    splitMethod->findThresh(numX, catX, Y, Trt, comX, this->fixIndex,
                            this->fitIndex, bestK, nodeFit->bestInds,
                            this->faster);

    if (splitMethod->getLoss() >= nodeFit->trainLoss)
    {
        //logger->debug(
        //             "Terminal due to loss not smaller. Original loss {}, split loss {}",
        //             nodeFit->trainLoss, splitMethod->getLoss());
        return result;
    }

    result->SplitID = splitMethod->getVarID();
    result->SplitRole = splitMethod->getVarRole();
    result->misDirection = splitMethod->getMissDir();

    if (result->SplitRole == 'n')
    {
        result->threshold = splitMethod->getThreshold();
        if (!is_finite(result->threshold) && result->misDirection != 'A')
        {
            //logger->debug("Terminal due to can not find threshold");
            return result;
        }
        //logger->debug("Split Var: {}, Split Threshold: {}", result->SplitID,
        //             result->threshold);
    }
    else
    {
        result->threshSet = splitMethod->getThreshSet();
    }

    uvec left = splitMethod->getSplitVec('L');
    uvec right = splitMethod->getSplitVec('R');
    // cout << "left :" << left.n_elem << "right: " << right.n_elem << endl;
    splitMethod->clear();
    result->terminal = false;

    result->left = buildTree(numX.rows(left), catX.rows(left), Y.rows(left),
                             Trt.rows(left), 2 * id, depth + 1);
    result->right = buildTree(numX.rows(right), catX.rows(right), Y.rows(right),
                              Trt.rows(right), 2 * id + 1, depth + 1);

    return result;
}

/**
         * Cross validation for tree
         * 1. First split data based on cvFold
         * 2. Generate tree based on each part of data
         * 3. evaluate tree based on test sample
         */
void GiTree::crossValidation(node *cvroot, const arma::mat &numX,
                             const arma::imat &catX, const arma::mat &Y,
                             const arma::ivec &Trt, const int &CVFold,
                             const double &CVSE)
{
    if (CVFold <= 5)
        return;

    if (!cvroot)
    {
        cerr << "Please run .fit function first";
        throw;
    }

    Node::evaluateTrain(cvroot);
    std::vector<node *> cvTreeMap;
    Node::updateNodeMap(cvroot, cvTreeMap);

    vec alphaVec = getAlphaVector(cvTreeMap);
    const uword &alphap = alphaVec.n_elem;
    vec nodeSize(alphap, arma::fill::zeros);

    for (int i = 0; i < alphap; i++)
    {
        alphaVec(i) = (i < alphap - 1) ? std::sqrt(alphaVec(i) * alphaVec(i + 1))
                                       : alphaVec(i) + 1.0;
        nodeSize(i) = Node::getTermSize(cvroot, alphaVec(i));
    }

    ivec SampleIndex =
        arma::randi<ivec>(Y.n_rows, arma::distr_param(0, CVFold - 1));
    mat CVRthat(CVFold, alphap, arma::fill::zeros);

    for (int i = 0; i < CVFold; i++)
    {
        Rcpp::checkUserInterrupt();
        //logger->debug("Start CV: {}", i);

        uvec trainIndex = arma::find(SampleIndex != i);
        uvec testIndex = arma::find(SampleIndex == i);

        node *CVRootTrain =
            this->buildTree(numX.rows(trainIndex), catX.rows(trainIndex),
                            Y.rows(trainIndex), Trt.rows(trainIndex), 1, 1);

        Node::evaluateTrain(CVRootTrain);

        for (auto j = 0; j < alphaVec.n_elem; j++)
        {
            Node::pruneAlpha(CVRootTrain, alphaVec(j));
            mat yhat = this->predictY(CVRootTrain, numX.rows(testIndex),
                                      catX.rows(testIndex), Trt.rows(testIndex));
            CVRthat(i, j) = RegSol::getLoss(Y.rows(testIndex), yhat, this->fitMethod);
        }
        Node::destroy(CVRootTrain);
    }

    vec aveLoss = arma::mean(CVRthat, 0).t();
    vec stdLoss = arma::stddev(CVRthat, 0).t() / std::sqrt(double(CVFold - 1));

    //logger->debug("Node Size: {}", nodeSize.t());
    //logger->debug("Average Loss: {}", aveLoss.t());
    //logger->debug("Std Loss: {}", stdLoss.t());

    auto minInd = arma::index_min(aveLoss);
    // Rcpp::Rcout << "aveLoss: " << aveLoss.t() << '\n';
    // Rcpp::Rcout << "CVSE: " << CVSE << ", aveLoss(minInd): " << aveLoss(minInd) << '\n';
    double minMSE = aveLoss(minInd) + CVSE * stdLoss(minInd);
    // Rcpp::Rcout << "minMSE: " << minMSE << '\n';

    const uvec &lossInd = arma::find(aveLoss <= minMSE + 1e-4);
    // Rcpp::Rcout << "aveLoss: " << aveLoss.t() << '\n';
    // Rcpp::Rcout << "nodeSize: " << nodeSize.t() << '\n';
    // Rcpp::Rcout << "AlphaVec: " << alphaVec.t() << '\n';
    // Rcpp::Rcout << "lossInd: " << lossInd.t() << '\n';
    auto mid = lossInd(arma::index_min(nodeSize(lossInd)));
    // if (mid == 0) mid++;
    // Rcpp::Rcout << "Mid: " << mid << '\n';
    // Rcpp::Rcout << "nodeSize(lossInd): "<< nodeSize(lossInd) << '\n';
    double minAlpha = alphaVec(mid);
    // Rcpp::Rcout << "minAlpha: " << minAlpha << '\n';
    Node::pruneAlpha(cvroot, minAlpha);

    //logger->debug("Minimal Alpha: {}", minAlpha);
    //logger->debug("Finish Cross-validation.");
}


mat GiTree::boostrapCI(const arma::mat &numX, const arma::imat &catX,
                       const arma::mat &Y, const arma::ivec &Trt,
                       const int &bootNum, const double &alpha)
{
    int K = 1000;
    const uword &yp = Y.n_cols;

    arma::vec alphaK = arma::linspace<arma::vec>(1e-9, alpha / 2.0, K);
    // arma::mat gamma(K, this->yp, arma::fill::zeros),
    //    theta(K, this->yp, arma::fill::zeros);
    arma::vec gamma(K, arma::fill::zeros), theta(K, arma::fill::zeros);

    // omp_set_num_threads(2);
    // # pragma omp parallel for schedule(static)
    for (int bi = 0; bi < bootNum; bi++)
    {
        Rcpp::checkUserInterrupt();
        // arma::mat gammaU(K, this->yp, arma::fill::zeros), thetaU(K, this->yp, arma::fill::ones);
        arma::vec gammaU(K, arma::fill::zeros), thetaU(K, arma::fill::ones);

        arma::uvec samInd =
            arma::randi<arma::uvec>(this->N, arma::distr_param(0, N - 1));
        const mat &numXBoot = numX.rows(samInd);
        const imat &catXBoot = catX.rows(samInd);
        const mat &YBoot = Y.rows(samInd);
        const ivec &TrtBoot = Trt.rows(samInd);

        node *bootRoot = GiTree::buildTree(numXBoot, catXBoot, YBoot, TrtBoot, 1, 1);

        GiTree::crossValidation(bootRoot, numXBoot, catXBoot, YBoot, TrtBoot,
                              this->CVFold, this->CVSE);

        std::vector<node *> bootTreeMap;
        Node::updateNodeMap(bootRoot, bootTreeMap);

        uvec bootNodeId = this->predictNode(bootRoot, numXBoot, catXBoot);
        uvec oriNodeId = this->predictNode(bootRoot, numX, catX);

        int TNodeSize = 0;
        for (auto &bootTerm : bootTreeMap)
        {
            if (!bootTerm->terminal)
                continue;

            uvec curInd = arma::find(oriNodeId == bootTerm->NodeID);
            const mat &numXcur = numX.rows(curInd);
            const mat &Ycur = Y.rows(curInd);
            const ivec &Trtcur = Trt.rows(curInd);

            mat comX = createFitMatrix(numXcur, Trtcur, this->trtLevel);

            std::vector<arma::vec> trueBeta = bootTerm->nodeModel->fitBeta(
                comX, Ycur, this->fixIndex, this->fitIndex, this->bestK,
                this->tp);
            std::vector<arma::mat> CoverMat =
                GetCoverMat(trueBeta, bootTerm->trtBeta, bootTerm->trtSE, alphaK);

            for (int k = 0; k < K; k++)
            {
                // gammaU.row(k) += arma::sum(CoverMat.at(k), 0);  // Count coverage
                // thetaU.row(k) %= arma::prod(CoverMat.at(k), 0); // Count all coverage
                gammaU(k) += arma::accu(CoverMat.at(k));
                thetaU(k) *= CoverMat.at(k).min();
            }
            TNodeSize++;
        }

        gamma += gammaU / (double)(TNodeSize * tp * yp);
        theta += thetaU;
        // if ((bi > 0) && (bi % 10 == 0)) logger->info("Finish {} bootstrap samples", bi);
        if ((bi > 0) && (bi % 10 == 0))
            cout << "finish " << bi << " bootstrap samples." << endl;
    }
    gamma /= (double)bootNum;
    theta /= (double)bootNum;
    // logger->debug("Gamma: \n{}\n Theta: {}\n", gamma, theta);
    return arma::join_rows(alphaK, arma::join_rows(gamma, theta));
}

mat GiTree::predictY(node *leaf, const mat &numX, const imat &catX,
                     const ivec &trt)
{
    const uword &N = trt.n_elem;
    mat result(N, this->yp, arma::fill::zeros);
    for (uword i = 0; i < N; i++)
    {
        node *pred = tree::predict(leaf, numX.row(i), catX.row(i));
        mat X = createFitMatrix(numX.row(i), trt.row(i), this->trtLevel);
        result.row(i) = pred->nodeModel->predict(X);
    }
    return result;
}

/**
 get Coverage matrix based on true value, center and SE

@param TrueTrt TRUE beta
@param center point estimate
@param sev standard error
@param Alpha alpha vector
@return coverage matrix 0 / 1
*/
std::vector<arma::mat> GetCoverMat(const std::vector<arma::vec> &TrueTrt,
                                   const std::vector<arma::vec> &center,
                                   const std::vector<arma::vec> &sev,
                                   const vec &Alpha)
{
    const int &N = TrueTrt.size();         // number of Y COL
    const uword &P = TrueTrt.at(0).n_elem; // number of beta ROW
    const uword &K = Alpha.n_elem;
    std::vector<arma::mat> CoverMat;

    boost::math::normal norm;
    for (uword k = 0; k < K; k++)
    {
        arma::mat covMat(P, N, arma::fill::zeros);
        double normLen = quantile(complement(norm, Alpha(k)));
        for (int col = 0; col < N; col++)
        {
            for (uword row = 0; row < P; row++)
            {
                const double &tru = TrueTrt.at(col)(row);
                const double &cen = center.at(col)(row);
                const double &se = sev.at(col)(row);
                double len = normLen * se;
                double upper = cen + len;
                double lower = cen - len;

                if ((tru <= upper) & (tru >= lower))
                    covMat(row, col) = 1.0;
            }
        }
        CoverMat.push_back(covMat);
    }
    return CoverMat;
}

mat GetBootAlpha(const vec &gamma, const vec &theta, const double &alphaL,
                 const vec &alphaK)
{
    const uword &yp = gamma.n_cols; // gamma each col is outcome, each row is
    // gamma k in original paper
    mat res(3, yp);
    res.fill(alphaL);

    for (int yi = 0; yi < yp; yi++)
    {

        uvec kpv = arma::find(gamma.col(yi) < 1 - alphaL, 1);
        uvec kqv = arma::find(theta.col(yi) < 1 - alphaL, 1);

        if (kpv.n_elem > 0)
        {
            uword kp = kpv(0);
            if (kp == 0)
            {
                cerr << "Consider increase bootstrap number!!!\n";
                gamma.print("Individual coverage rate!");
                return res;
            }
            double f = (gamma(kp - 1, yi) - 1.0 + alphaL) /
                       (gamma(kp - 1, yi) - gamma(kp, yi));
            res(1, yi) = (1.0 - f) * alphaK(kp - 1) + f * alphaK(kp);
        }

        if (kqv.n_elem > 0)
        {
            uword kq = kqv(0);
            if (kq == 0)
            {
                cerr << "Consider increase bootstrap number!!!\n";
                theta.print("Overall coverage rate!");
                return res;
            }
            double g = (theta(kq - 1, yi) - 1.0 + alphaL) /
                       (theta(kq - 1, yi) - theta(kq, yi));
            res(2, yi) = (1.0 - g) * alphaK(kq - 1) + g * alphaK(kq);
        }
    }
    return res;
}

inline std::vector<uvec> refineFit(const std::vector<uvec> &x, const int &np)
{
    std::vector<uvec> result;
    for (auto &item : x)
    {
        result.push_back(item.elem(arma::find(item < np)));
    }
    return result;
}

inline std::vector<vec> refineTrt(const std::vector<RegSol::RegParm> &parms,
                                  const int &tp, const bool &beta)
{
    std::vector<vec> result;
    for (auto &item : parms)
    {
        if (beta)
        {
            result.push_back(item.beta.tail(tp));
        }
        else
        {
            result.push_back(item.betaSE.tail(tp));
        }
    }
    return result;
}

inline arma::mat createFitMatrix(const arma::mat &numX, const arma::ivec &trt, const arma::ivec &trtLevel)
{
    const arma::uword &np = numX.n_cols;
    const arma::uword &tp = trtLevel.n_elem;
    const arma::mat &tmp = hotCoding(trt, trtLevel, false);

    // combine X: numX + 1 + trt2 + .. + trt_tp
    arma::mat comX(tmp.n_rows, np + tp);

    for (arma::uword i = 0; i < np + tp; i++)
    {
        if (i < np)
        {
            comX.col(i) = numX.col(i);
        }
        else if (i == np)
        {
            comX.col(i).fill(1);
        }
        else
        {
            comX.col(i) = tmp.col(i - np);
        }
    }
    return comX;
}
} // namespace Tree
} // namespace SubGuide
