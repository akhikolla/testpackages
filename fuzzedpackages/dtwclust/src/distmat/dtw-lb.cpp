#include "R-gateways.h"

#include <cstddef> // std::size_t
#include <memory> // shared_ptr

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/calculators.h"
#include "../utils/ParallelWorker.h"
#include "../utils/utils.h" // get_grain, id_t

namespace dtwclust {

// =================================================================================================
/* worker to update DTW distance in parallel */
// =================================================================================================

class DtwDistanceUpdater : public ParallelWorker {
public:
    // constructor
    DtwDistanceUpdater(const SurrogateMatrix<bool>& id_changed,
                       const SurrogateMatrix<int>& id_nn,
                       Rcpp::NumericMatrix& distmat,
                       const std::shared_ptr<DistanceCalculator>& dist_calculator,
                       int margin,
                       const int grain)
        : ParallelWorker(grain, 1000, 10000)
        , id_changed_(id_changed)
        , id_nn_(id_nn)
        , distmat_(distmat)
        , dist_calculator_(dist_calculator)
        , margin_(margin)
    { }

    // parallel loop across specified range
    void work_it(std::size_t begin, std::size_t end) override {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // update distances
        if (margin_ == 1) {
            for (std::size_t i = begin; i < end; i++) {
                if (is_interrupted(i)) break; // nocov

                if (id_changed_[i]) {
                    int j = id_nn_[i];
                    distmat_(i,j) = dist_calculator->calculate(i,j);
                }
            }
        }
        else {
            for (std::size_t j = begin; j < end; j++) {
                if (is_interrupted(j)) break; // nocov

                if (id_changed_[j]) {
                    int i = id_nn_[j];
                    distmat_(i,j) = dist_calculator->calculate(i,j);
                }
            }
        }

        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    // input vectors
    const SurrogateMatrix<bool>& id_changed_;
    const SurrogateMatrix<int>& id_nn_;
    // output matrix
    RcppParallel::RMatrix<double> distmat_;
    // distance calculator
    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    // margin for update
    const int margin_;
};

// =================================================================================================
/* find nearest neighbors */
// =================================================================================================

void set_nn(const Rcpp::NumericMatrix& distmat, SurrogateMatrix<int>& nn, const int margin)
{
    if (margin == 1) {
        for (int i = 0; i < distmat.nrow(); i++) {
            double d = distmat(i,0);
            nn[i] = 0;
            for (int j = 1; j < distmat.ncol(); j++) {
                double temp = distmat(i,j);
                if (temp < d) {
                    d = temp;
                    nn[i] = j;
                }
            }
        }
    }
    else {
        for (int j = 0; j < distmat.ncol(); j++) {
            double d = distmat(0,j);
            nn[j] = 0;
            for (int i = 1; i < distmat.nrow(); i++) {
                double temp = distmat(i,j);
                if (temp < d) {
                    d = temp;
                    nn[j] = i;
                }
            }
        }
    }
}

// =================================================================================================
/* check if updates are finished based on indices */
// =================================================================================================

bool check_finished(const SurrogateMatrix<int>& nn,
                    const SurrogateMatrix<int>& nn_prev,
                    SurrogateMatrix<bool>& changed)
{
    bool finished = true;
    for (id_t i = 0; i < nn.nrow(); i++) {
        if (nn[i] != nn_prev[i]) {
            changed[i] = true;
            finished = false;
        }
        else {
            changed[i] = false;
        }
    }
    return finished;
}

// =================================================================================================
/* main C++ function */
// =================================================================================================

void dtw_lb_cpp(const Rcpp::List& X,
                const Rcpp::List& Y,
                Rcpp::NumericMatrix& distmat,
                const SEXP& DOTS,
                const int margin,
                const int num_threads)
{
    auto dist_calculator = DistanceCalculatorFactory().create("DTW_BASIC", DOTS, X, Y);

    int len = margin == 1 ? distmat.nrow() : distmat.ncol();
    SurrogateMatrix<int> id_nn(len, 1);
    SurrogateMatrix<int> id_nn_prev(len, 1);
    SurrogateMatrix<bool> id_changed(len, 1);

    int grain = get_grain(len, num_threads);
    DtwDistanceUpdater dist_updater(id_changed, id_nn, distmat, dist_calculator, margin, grain);

    set_nn(distmat, id_nn, margin);
    for (id_t i = 0; i < id_nn.nrow(); i++) id_nn_prev[i] = id_nn[i] + 1; // initialize different

    while (!check_finished(id_nn, id_nn_prev, id_changed)) {
        // update nn_prev
        for (id_t i = 0; i < id_nn.nrow(); i++) id_nn_prev[i] = id_nn[i];
        // calculate dtw distance if necessary
        parallel_for(0, len, dist_updater, grain);
        // update nearest neighbors
        set_nn(distmat, id_nn, margin);
    }
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

extern "C" SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    Rcpp::NumericMatrix distmat(D);
    dtw_lb_cpp(X, Y, distmat, DOTS, Rcpp::as<int>(MARGIN), Rcpp::as<int>(NUM_THREADS));
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
