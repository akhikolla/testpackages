//
// Created by Peigen Zhou on 7/21/18.
//

#include "split.hpp"

namespace SubGuide {

    namespace SplitSol {
        using arma::is_finite;

        double chiApproximate(const double &ss_t, const double &ss_pe, const int &df_t,
                              const int &df_pe);

        inline void updateLRres(RegSol::RegParm &LRes, RegSol::RegParm &RRes, const uvec &left, const uvec &right, RegSol::RegFun *fitMethod, const mat &conX, const mat &comX, const vec &yi,  const uvec &fixIndex, const uvec &fitIndex, const int& bestK,  const bool &faster);

        // Split::Split() : logger(spdlog::stdout_color_st("Split")) {}
        Split::Split() {}

        Split::Split(const int &batchNum_, const int &minData_)
        : batchNum(batchNum_), minData(minData_) {
            //logger = spdlog::get("split");
        }

        void Split::clear() {
            varID = {};
            threshold = 0.0;
            threshSet.clear();
            loss = -1.0;
            this->optLeft.clear();
            this->optRight.clear();
        }

        void Split::dataCheck(const mat &nx_, const imat &cx_, const mat &ys_) {
            N = ys_.n_rows;
            np = nx_.n_cols;
            cp = cx_.n_cols;
            yp = ys_.n_cols;

            chiN.resize(np);
            chiC.resize(cp);
            chiN.zeros();
            chiC.zeros();

            // assert(nx_.n_rows == ys_.n_rows);
            // assert(cx_.n_rows == ys_.n_rows);
        }

        GiSplit::GiSplit(RegSol::RegFun *fitMethod, const int &batchNum_,
                         const int &minData_, const int &minTrt_, const ivec &trtLevel_)
        : Split(batchNum_, minData_), minTrt(minTrt_), nodeFitMethod(fitMethod), trtLevel(trtLevel_) {}

        void GiSplit::dataChecking(const mat &nx_, const imat &cx_, const mat &ys_,
                                   const ivec trt_) {
            Split::dataCheck(nx_, cx_, ys_);
            // assert(ys_.n_rows == trt_.n_rows);
            // logger->info("Data checking GiSplit passed;");

            // if (!arma::approx_equal(trtLevel, arma::unique(trt_), "absdiff", 0)) {
            //     Rcpp::Rcerr << "Error in GiSplit TrtLevel not match \n" << "Original: " << trtLevel.t() << "Now: " << arma::unique(trt_).t() << "\n" << std::endl;
                 // abort();
            // }
            this->trtDes = hotCoding(trt_, false); // TODO: trtDesign whether need this
        }


        /**
         find Split variable based on current data set.

         @param numX numerical X // TODO: Should numX first column is intercept?
         @param catX categorical X
         @param Y outcome matrix, n * yp
         @param trt treatment indicator
         @param bestInd_ step wise result, best fitted x index in numX
         @param SplitIndex splitIndex, variable can used to split data
         */
        void GiSplit::findSplit(const mat &numX, const imat &catX, const mat &Y,
                                const icolvec &trt, const std::vector<uvec> &bestInd_,
                                const uvec &SplitIndex) {

            this->dataChecking(numX, catX, Y, trt);
            // assert(trt.n_rows == Y.n_rows);

            // double maxChi = 0;

            for (auto i = 0; i < yp; i++) {
                //logger->info("Start {}th outcomes", i);

                const vec& Yi = Y.col(i);
                uvec bInd = bestInd_.at(i);
                // logger->info("bInd: {}", bInd.t());

                for (auto jn = 0; jn < SplitIndex.n_elem; jn++) {

                    // logger->info("Start {}th numerical variable", jn);

                    const vec &nX = numX.col(SplitIndex(jn));

                    if (!arma::is_finite(nX)) {
                        const arma::uvec &ninf_v = arma::find_finite(nX);
                        if (!checkNodeData(ninf_v, trt(ninf_v), this->trtLevel, this->minData, this->minTrt)) {
                            // logger->debug("Skip {}th numerical variable, due to missing data is too much", jn);
                            continue;
                        }
                    }

                    ivec cNx = quartileX(nX, 4);
                    mat DcNx = hotCoding(cNx, true);
                    double lchi = lackOfFit(numX.cols(bInd), DcNx, Yi, this->trtDes);
                    this->chiN(SplitIndex(jn)) += lchi;
                    // Rcpp::Rcout << "jn: " << jn << ", cp: " << cp << ", lackoffit" << lchi << '\n';
                    /*
                    if (this->chiN(jn) > maxChi) {
                        maxChi = this->chiN(jn);
                        this->varID = SplitIndex(jn);
                        this->role = 'n';
                    }
                     */
                }
                // logger->debug("ChiN: {}", chiN.t());
                // logger->info("Finish numerical variables in {}th outcome", i);

                for (auto jc = 0; jc < cp; jc++) {
                    // logger->info("Start {}th categorical variable", jc);

                    const ivec& Cx = catX.col(jc);
                    ivec levels = arma::unique(Cx);

                    if (levels.n_elem == 1) {
                        // logger->info("Skip due to unique levels");
                        this->chiC(jc) += 0.0;
                    } else {
                        mat DcCx = hotCoding(Cx, true);
                        double lchi = lackOfFit(numX.cols(bInd), DcCx, Yi, this->trtDes);
                        this->chiC(jc) += lchi;
                        // Rcpp::Rcout << "bind: " << bInd.t() << '\n';
                        // Rcpp::Rcout << "jc: " << jc << ", cp: " << cp << ", lackoffit" << lchi << '\n';
                    }
                    /*
                    if (this->chiC(jc) > maxChi) {
                        maxChi = this->chiC(jc);
                        this->varID = np + jc;
                        this->role = 'c';
                    }
                     */
                }
                // logger->debug("ChiC: {}", chiC.t());
            }

            if ((this->np > 0) && (this->cp > 0)) {
                bool wNC = chiN.max() >= chiC.max();
                // logger->info("wNC {}", wNC);
                this->varID = wNC ? arma::index_max(chiN) : arma::index_max(chiC) + np;
                this->role = wNC ? 'n' : 'c';

            } else {
                this->varID = this->np != 0 ? arma::index_max(chiN) : arma::index_max(chiC) + np;
                this->role = this->np != 0 ? 'n' : 'c';
            }

            // logger->info("The selected variable ID: {}, role: {}", varID, role);
        }

        void GiSplit::findThresh(const mat &numX, const imat &catX, const mat &Y, const ivec trt,
                       const mat &comX, const uvec &fixIndex,
                       const uvec &fitIndex, const int &bestK,
                       const std::vector<uvec> &bestInd_, const bool& faster) {
            this->dataChecking(numX, catX, Y, trt);
            // assert(trt.n_rows == Y.n_rows);
            this->optLeft.clear();
            this->optRight.clear();

            if (role == 'n') {
                // logger->debug("Start Finding best threshold: {}", varID);
                this->threshold = findNumThresh(numX.col(varID), trt, Y, comX, fixIndex, fitIndex, bestK, bestInd_, faster);
                // logger->info("Best threshold: {}", threshold);
            }

            if (role == 'c') {
                // logger->debug("Start Finding best threshSet: {}", varID - np);
                this->threshSet = findCateThresh(catX.col(varID - np), trt, Y, comX, fixIndex, fitIndex, bestK, bestInd_, faster);
            }
        }

        /**
         Find categorical threshset, missing data is recoded as 123456789

         @param x categorical x variable
         @param numX numerical X used in the fitted function
         @param trt treatment variable
         @param ys ydesing
         @param bestInd_ fitted variable for numX
         @return return threshSet ivce
         */
        ivec GiSplit::findCateThresh(const ivec &x, const ivec trt,
                                     const mat &ys, const mat &comX,
                                     const uvec &fixIndex, const uvec &fitIndex, const int &bestK,
                                     const std::vector<uvec> &bestInd_, const bool& faster) {
            ivec levels = arma::unique(x);
            umat xL = getLevels(x); // all combination of categorical set

            mat threshLoss(xL.n_rows, 2, arma::fill::zeros);

            for (auto j = 0; j < yp; j++) {

                const vec &yi = ys.col(j);

                uvec bInd = bestInd_.at(j);
                const mat& conX = comX.cols(bInd);

                for (auto i = 0; i < xL.n_rows; i++) {
                    ivec threshSetTmp = levels.elem(arma::find(xL.row(i) == 1));
                    uvec ind = match(x, threshSetTmp);

                    uvec left = arma::find(ind == 1);
                    if (!checkNodeData(left, trt(left), trtLevel, minData, minTrt))
                        continue;

                    uvec right = arma::find(ind == 2);
                    if (!checkNodeData(right, trt(right), trtLevel, minData, minTrt))
                        continue;

                    // assert(left.n_elem >= minData);
                    // assert(right.n_elem >= minData);

                    RegSol::RegParm LRes, RRes;

                    updateLRres(LRes, RRes, left, right, this->nodeFitMethod, conX, comX, yi, fixIndex, fitIndex, bestK, faster);

                    threshLoss(i, 0) += LRes.loss + RRes.loss;
                    threshLoss(i, 1) = i;
                }
            }
            threshLoss = threshLoss.rows(arma::find(threshLoss.col(0) > 0));

            if (threshLoss.n_rows == 0) {
                loss = arma::datum::inf;
                this->threshSet = {};
                return this->threshSet;
            }

            loss = arma::as_scalar(arma::min(threshLoss.col(0)));
            const uword &ind = arma::index_min(threshLoss.col(0));

            this->threshSet = levels.elem(arma::find(xL.row(threshLoss(ind, 1)) == 1));
            // Rcpp::Rcout << "Levels: " << levels << '\n';
            // Rcpp::Rcout << "Threshset" << this->threshSet << '\n';
            uvec Tind = SubGuide::match(x, this->threshSet);

            this->optLeft = arma::find(Tind == 1);
            this->optRight = arma::find(Tind == 2);
            // Rcpp::Rcout << "Left size: " << this->optLeft.n_elem << '\n';
            // Rcpp::Rcout << "Right size: " << this->optRight.n_elem << '\n';

            if (arma::any(arma::find(levels == miss))) {
                this->misDirection = arma::any(arma::find(this->threshSet == miss)) ? 'L' : 'R';
            } else {
                this->misDirection = this->optLeft.n_elem >= this->optRight.n_elem ? 'L' : 'R';
            }

            // assert(this->optLeft.n_elem >= minData);
            // assert(this->optRight.n_elem >= minData);

            return this->threshSet;
        }

        /**
         Find threshold to cut the data based on numerical variable X.
         Missing Data need a special treatment.

         1. Missing Data in one side. Also need to pass the checkNodeData
         function.
         2. Missing Data left and loop
         3. Missing Data right and loop

         - parameters:

         1. x: numerical variable
         2. numX: control variable
         3. trt: treatment assignment
         4. ys: outcome variables
         5. bestInd_: prefound best indicators

         - return:
         threshold for numerical variable.
         */
        double GiSplit::findNumThresh(const vec &x, const ivec trt, const mat &ys,
                       const mat &comX, const uvec &fixIndex,
                       const uvec &fitIndex, const int &bestK,
                       const std::vector<uvec> &bestInd_, const bool& faster) {
            uvec splitSort = arma::sort_index(x);
            const uvec missInd = arma::find_nonfinite(x);
            const uvec nonMiss = arma::find_finite(x);
            const uword &nonN = nonMiss.n_elem;
            bool missOnSide =
            checkNodeData(missInd, trt(missInd), trtLevel, minData, minTrt);

            bool missStatue = missInd.n_elem > 0;
            double missLoss = missStatue ? -1.0 : 0.0;

            // threshLossL:
            // 1. no missing data, treshloss holding
            // 2. Missing all left holding
            // threshLossR: Missing all right holding
            mat threshLossL(N, 2), threshLossR(N, 2);
            threshLossL.fill(-1.0);
            threshLossR.fill(-1.0);

            for (auto j = 0; j < yp; j++) {
                const vec &yi = ys.col(j);
                uvec bInd = bestInd_[j];
                const mat &conX = comX.cols(bInd);

                if (missOnSide) {
                    RegSol::RegParm LRes, RRes;

                    updateLRres(LRes, RRes, missInd, nonMiss, this->nodeFitMethod, conX, comX, yi, fixIndex, fitIndex, bestK, faster);

                    missLoss += LRes.loss + RRes.loss;
                    // logger->debug("missing onside Loss: {}", missLoss);
                }

                for (int i = 1; i < N - missInd.n_elem; i += batchNum) {
                    uvec left = splitSort.head(i);
                    uvec right = splitSort(arma::span(i, nonN - 1));
                    // uvec right = splitSort.tail(N - missInd.n_elem - i);
                    if (arma::approx_equal(x.row(left[i - 1]), x.row(right[0]), "absdiff", 1e-5)) continue;
                    if (!missStatue) {
                        if (!checkNodeData(left, trt(left), trtLevel, minData, minTrt))
                            continue;
                        if (!checkNodeData(right, trt(right), trtLevel, minData, minTrt))
                            continue;

                        RegSol::RegParm LRes, RRes;
                        updateLRres(LRes, RRes, left, right, this->nodeFitMethod, conX, comX, yi, fixIndex, fitIndex, bestK, faster);

                        threshLossL(i, 0) += LRes.loss + RRes.loss;
                        threshLossL(i, 1) = i;

                        // logger->debug("Current i: {}, x ind: {}, loss: {}, Left: {}", i,
                        //              threshLossL(i, 1), threshLossL(i, 0), left.n_elem);
                    } else {
                        const uvec &leftNew = join_cols(left, missInd);
                        const uvec &rightNew = join_cols(right, missInd);

                        RegSol::RegParm LRes, RRes;

                        if (checkNodeData(leftNew, trt(leftNew), trtLevel, minData, minTrt) &&
                            checkNodeData(right, trt(right), trtLevel, minData, minTrt)) {

                            updateLRres(LRes, RRes, leftNew, right, this->nodeFitMethod, conX, comX, yi, fixIndex, fitIndex, bestK, faster);

                            threshLossL(i, 0) += LRes.loss + RRes.loss;
                            threshLossL(i, 1) = i;
                            // logger->debug("Current i: {}, x ind: {}, loss: {}, Left: {}", i,
                            //              threshLossL(i, 1), threshLossL(i, 0), left.n_elem);
                        }

                        if (checkNodeData(left, trt(left), trtLevel, minData, minTrt) &&
                            checkNodeData(rightNew, trt(rightNew), trtLevel, minData, minTrt)) {

                            updateLRres(LRes, RRes, left, rightNew, this->nodeFitMethod, conX, comX, yi, fixIndex, fitIndex, bestK, faster);

                            threshLossR(i, 0) += LRes.loss + RRes.loss;
                            threshLossR(i, 1) = i;
                            // logger->debug("Current i: {}, x ind: {}, loss: {}, right: {}", i,
                            //               threshLossR(i, 1), threshLossR(i, 0), right.n_elem);
                        }
                    }
                }
            }

            threshLossL = threshLossL.rows(arma::find(threshLossL.col(0) > 0));
            // logger->debug("missStatus: {}", missStatue);

            /**
             If has missing value:
             1. check one side loss whether smaller than threshLossL and threshLossR
             Yes: return misDIrection = 'A';
             No: if threshLossL give smaller loss: find the splitSort index i; append left
             and missind if no missing value:
             */
            if (missStatue) {
                threshLossR = threshLossR.rows(arma::find(threshLossR.col(0) > 0));

                double lossL, lossR;
                lossL = threshLossL.n_rows > 0
                ? arma::as_scalar(arma::min(threshLossL.col(0)))
                : arma::datum::inf;
                lossR = threshLossR.n_rows > 0
                ? arma::as_scalar(arma::min(threshLossR.col(0)))
                : arma::datum::inf;

                this->loss = lossL < lossR ? lossL : lossR;

                if ((missLoss > 0.0) && (missLoss < this->loss)) {
                    this->misDirection = 'A';
                    this->loss = missLoss;
                    this->threshold = arma::datum::inf;
                    this->optLeft = missInd;
                    this->optRight = nonMiss;
                    return this->threshold;
                } else {
                    if (!is_finite(lossL) || !is_finite(lossR)) {
                        this->loss = arma::datum::inf;
                        this->threshold = arma::datum::inf;
                        return this->threshold;
                    }

                    int minIndT;
                    if (lossL < lossR) {
                        this->misDirection = 'L';
                        minIndT = threshLossL(arma::index_min(threshLossL.col(0)), 1);
                        this->optLeft = join_cols(splitSort.head(minIndT), missInd);
                        this->optRight = splitSort(arma::span(minIndT, nonN - 1));
                        // this->optRight = splitSort.tail(N - missInd.n_elem - minIndT);

                    } else {
                        this->misDirection = 'R';
                        minIndT = threshLossR(arma::index_min(threshLossR.col(0)), 1);
                        this->optLeft = splitSort.head(minIndT);
                        this->optRight =
                        join_cols(splitSort(arma::span(minIndT, nonN - 1)), missInd);
                    }
                    this->threshold = (x(splitSort(minIndT - 1)) + x(splitSort(minIndT))) / 2.0;
                }
            } else {
                this->misDirection = 'L';
                if (threshLossL.n_rows == 0) {
                    this->loss = arma::datum::inf;
                    this->threshold = arma::datum::inf;
                    return this->threshold;
                } else {
                    this->loss = arma::as_scalar(arma::min(threshLossL.col(0)));
                    int minIndT = threshLossL(arma::index_min(threshLossL.col(0)), 1);
                    this->threshold = (x(splitSort(minIndT - 1)) + x(splitSort(minIndT))) / 2.0;

                    this->optLeft = arma::find(x <= this->threshold);
                    this->optRight = arma::find(x > this->threshold);

                    misDirection = this->optLeft.n_elem > this->optRight.n_elem ? 'L' : 'R';

                    if (!checkNodeData(this->optLeft, trt(this->optLeft), trtLevel, minData, minTrt) ||
                        !checkNodeData(this->optRight, trt(this->optRight), trtLevel, minData, minTrt)) {
                        this->loss = arma::datum::inf;
                        this->threshold = arma::datum::inf;
                        return this->threshold;
                    }
                }
            }

            // logger->debug("loss: {}. this->optLeft n: {}, this->optRight n: {}", loss, this->optLeft.n_elem,
            //               this->optRight.n_elem);

            // assert(checkNodeData(this->optLeft, trt(this->optLeft), trtLevel, minData, minTrt));
            // assert(checkNodeData(this->optRight, trt(this->optRight), trtLevel, minData, minTrt));
            // assert(this->optLeft.n_elem >= minData);
            // assert(this->optRight.n_elem >= minData);

            return threshold;
        }

        /**
         Lack of fit test of split variable selection
         missing value will be imputed in Main

         @param Main First part of data
         @param cx X design
         @param Y_ outcome vector
         @param trtDes treatment design matrix
         @return chi-square value
         */
        double GiSplit::lackOfFit(const mat &Main, const mat &cx, const vec &Y_, const mat &trtDes) {
            // assert(trtDes.n_cols > 0);

            mat xTotal = is_finite(Main)
            ? join_rows(Main, join_rows(cx, trtDes))
            : join_rows(imputeMean(Main), join_rows(cx, trtDes));

            const mat &tmp = designInt(cx, trtDes.tail_cols(trtDes.n_cols - 1));
            mat xPe = join_rows(xTotal, tmp);

            if (xPe.n_cols == xTotal.n_cols) {
                return 0.0;
            }

            if (xPe.n_rows <= xPe.n_cols) {
                return 0.0;
            }

            RegSol::RegParm parmTotal = nodeFitMethod->fit(xTotal, Y_);
            RegSol::RegParm parmPe = nodeFitMethod->fit(xPe, Y_);

            if (parmPe.loss > parmTotal.loss)
                return 0.0;

            if (parmPe.loss == arma::datum::inf || parmTotal.loss == arma::datum::inf) 
                return 0;
            
            // Rcpp::Rcout << "Total loss: " << parmTotal.loss << ", Pe loss: " << parmPe.loss << ", Total df: " << parmTotal.df << ", Pe df: " << parmPe.df << '\n';
            double result =
            chiApproximate(parmTotal.loss, parmPe.loss, parmTotal.df + 2, parmPe.df + 1);
            // logger->info("Chi-value: {}", result);
            return result;
        }

        double chiApproximate(const double &ss_t, const double &ss_pe, const int &df_t,
                              const int &df_pe) {
            double chi_value;
            double ss_lof = ss_t - ss_pe;
            int df_lof = df_t - df_pe;
            double u = df_lof;
            double v = df_pe;
            double phi = v / (v - 2.0);
            double tau = 2.0 * std::pow(v, 2) +
            (u + v - 2.0) / (u * std::pow(v - 2.0, 2) * (v - 4.0));
            double F_value = (ss_lof / u) / (ss_pe / v);
            bool con1 = (v < 10.0) && (F_value < 3000.0 * tau + phi);
            bool con2 = (v >= 10.0) && (F_value < 150.0 * tau + phi);

            boost::math::fisher_f f_dist(df_lof, df_pe);
            boost::math::chi_squared chi_dist(1);

            if (con1 || con2) {
                double p_value = boost::math::cdf(f_dist, F_value);
                if (std::fabs(p_value - 1.0) < 1e-8)
                    return 1e10;
                chi_value = boost::math::quantile(chi_dist, p_value);
                // logger->info("\ndf_lof: {} df_pe: {}\n F_value: {}, p_value:
                // {}\n chi_value: {}", df_lof, df_pe, F_value, p_value,
                // chi_value);
            } else {
                double a = u * F_value / 3.0;
                double b = (2.0 * v + a + u - 2.0) / (2.0 * (v + 2.0 * a));
                chi_value = b * v * F_value;
                if (std::fabs(v - 1.0) < 1e-8)
                    return chi_value;

                double w1;
                w1 = std::pow(std::sqrt(2.0 * chi_value) - std::sqrt(2.0 * v - 1.0) + 1.0,
                              2) /
                2.0;
                double w2;
                w2 = std::max(0.0,
                              std::pow(7.0 / 9.0 + std::sqrt(v) *
                                       (std::pow(chi_value / v, 1.0 / 3.0) -
                                        1.0 + 2.0 / (9.0 * v)),
                                       3));
                chi_value = chi_value < v + 10.0 * std::sqrt(2.0 * v)
                ? w2
                : chi_value > w2 ? (w1 + w2) / 2.0 : w1;
            }
            return chi_value;
        }

        inline bool GiSplit::checkNodeData(const uvec &index, const ivec &trt,
                                           const ivec &level, const int &minData,
                                           const int &minTrt) {
            // assert(index.n_elem == trt.n_elem);

            if (index.n_elem <= minData)
                return false;

            if (trt.n_elem < level.n_elem * minTrt)
                return false;

            for (const auto &item : level) {
                uvec tmpInd = arma::find(trt == item);
                if (tmpInd.n_elem < minTrt)
                    return false;
            }
            return true;
        }

        inline void updateLRres(RegSol::RegParm &LRes, RegSol::RegParm &RRes, const uvec &left, const uvec &right, RegSol::RegFun *fitMethod, const mat &conX, const mat &comX, const vec &yi,  const uvec &fixIndex, const uvec &fitIndex, const int& bestK, const bool &faster) {

            if (faster) {
            LRes = fitMethod->fit(imputeMean(conX.rows(left)), yi(left));
            RRes = fitMethod->fit(imputeMean(conX.rows(right)), yi(right));
            } else {
                arma::uvec Btmp;
                arma::mat completeX;

                completeX = imputeMean(comX.rows(left));
                LRes = RegSol::stepWiseF(fitMethod, completeX, yi(left), fixIndex, fitIndex, bestK, Btmp);

                completeX = imputeMean(comX.rows(right));
                RRes = RegSol::stepWiseF(fitMethod, completeX, yi(right), fixIndex, fitIndex, bestK, Btmp);
            }
        }


        /*
         uvec GiSplit::numStepIndex(const mat &comX, const uvec &offInd, int K,
         const vec &Y_) {

         RegSol::stepWise stepRes(nodeFitMethod, comX, Y_, offInd, K);
         uvec bestInd = stepRes.bestInd;
         bestInd      = bestInd.elem(arma::find(bestInd < np));

         // bestInd.print("Best fitted varID: \n");
         return bestInd;
         }
         */

    } // namespace SplitSol
} // namespace SubGuide
