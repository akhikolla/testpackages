//
// Created by Peigen Zhou on 7/22/18.
//

#include "utils.hpp"

namespace SubGuide {
    
    arma::vec quantile(const arma::vec& X, const int & part) {
        // assert(arma::is_finite(X));
        const arma::uword N = X.n_elem;
        // assert(N > part);
        
        arma::vec result(part);
        arma::vec quartile = arma::linspace(0.0, 1.0, part + 1);
        arma::vec sortX = arma::sort(X);
        result[0] = arma::min(X);
        
        for (auto i = 1; i < part; i++) {
            const double & p = quartile[i];
            double h = (N + 1.0) * p;
            int hd = static_cast<int>(std::floor(h));
            result[i] = sortX[hd - 1] + (h - hd)*(sortX[hd] - sortX[hd - 1]);
        }
        
        return arma::unique(result);
    }

    /**
     Divided X into parts based on sample quartile, missing data will add to part + 1

     @param x numerical variable used to discritzied
     @param parts number of part
     @return integer vector
     */
    arma::ivec quartileX(const arma::vec &x, int parts) {
        int n = x.n_elem;
        arma::ivec label_index(n);
        arma::uvec inf_v = arma::find_nonfinite(x);
        arma::uvec fin_v = arma::find_finite(x);
        
        label_index.fill(parts);
        label_index.rows(inf_v).fill(parts + 1);
        
        if (fin_v.n_elem <= parts) {
            return label_index;
        }
        
        arma::vec quar = quantile(x.rows(fin_v), parts);
        for (int i = 1; i < quar.n_elem; i++) {
            arma::uvec index;
            if (i == 1) {
                index = arma::find((x <= quar[i]) && (x >= quar[i - 1]));
            } else {
                index = arma::find((x <= quar[i]) && (x > quar[i - 1]));
            }
            label_index.rows(index).fill(i);
        }
        return label_index;
    }
    
    arma::mat hotCoding(const arma::ivec &cx, const arma::ivec &levels,
                        bool reference) {
        const auto &n = cx.n_elem;
        const auto &p = reference ? levels.n_elem - 1 : levels.n_elem;
        arma::uword start = reference ? 1 : 0;
        arma::mat loading;
        arma::uvec rowIndex;
        arma::uvec colIndex;
        
        loading = arma::zeros<arma::mat>(n, p);
        
        for (arma::uword i = 0; i < p; start++, i++) {
            rowIndex = arma::find(cx == levels[start]);
            colIndex = {i};
            loading.elem(rowIndex, colIndex).ones();
        }
        return loading;
    }
    
    arma::mat hotCoding(const arma::ivec &cx, bool reference) {
        const arma::ivec &levels = arma::unique(cx);
        return hotCoding(cx, levels, reference);
    }
    
    arma::mat designInt(const arma::mat &x1, const arma::mat &x2) {
        // assert(x1.n_rows == x2.n_rows);
        
        const auto &n = x1.n_rows;
        const auto &p1 = x1.n_cols;
        const auto &p2 = x2.n_cols;
        
        arma::mat loading(n, p1 * p2);
        loading.zeros();
        arma::uword index = 0;
        
        for (auto i = 0; i < p1; i++) {
            for (auto j = 0; j < p2; j++) {
                loading.col(index) = x1.col(i) % x2.col(j);
                index++;
            }
        }
        const arma::uvec &nonZero = arma::find(arma::sum(loading, 0) != 0);
        return loading.cols(nonZero);
    }
    
    arma::umat getLevels(const arma::ivec &cx) {
        arma::ivec levels = arma::unique(cx);
        int nl = levels.n_elem;
        int mi = std::pow(2, nl - 1) - 1;
        arma::umat result(mi, nl);
        for (int i = 0; i < mi; i++) {
            int ii = i + 1;
            for (int j = 0; j < nl; j++) {
                result(i, j) = ii % 2L == 0 ? 0 : 1;
                ii /= 2L;
            }
        }
        return result;
    }
    
    arma::uvec match(const arma::ivec &cx, const arma::ivec &Xset) {
        const arma::uword &n = cx.n_elem;
        arma::uvec result(n);
        for (auto i = 0; i < n; i++)
            result(i) = arma::any(Xset == cx(i)) ? 1 : 2;
        return result;
    }
    
    arma::vec colMean(const arma::mat &X) {
        const arma::uword &p = X.n_cols;
        arma::vec result(p);
        
        if (arma::is_finite(X)) {
            return arma::mean(X, 0).t();
        } else {
            for (auto i = 0; i < p; i++) {
                if (arma::is_finite(X.col(i))) {
                    result[i] = arma::mean(X.col(i));
                    continue;
                } else {
                    const arma::vec &a = X.col(i);
                    const arma::uvec &find = arma::find_finite(a);
                    if (find.n_elem > 0)
                        result[i] = arma::mean(a(find));
                }
            }
        }
        return result;
    }
    
    arma::mat imputeValue(const arma::mat &X, const arma::vec &Xmean) {
        // assert(X.n_cols == Xmean.n_elem);
        arma::mat result(X);
        int i = 0;
        result.each_col([&i, &Xmean](arma::vec &a) {
            if (arma::is_finite(a)) {
                i++;
                return;
            }
            if (arma::is_finite(Xmean[i])) {
                a.elem(arma::find_nonfinite(a)).fill(Xmean[i]);
                i++;
            } else {
                return;
            }
        });
        return result;
    }
    
    arma::mat imputeMean(const arma::mat &X) {
        if (arma::is_finite(X)) return X;
        
        arma::vec Xmean = colMean(X);
        return imputeValue(X, Xmean);
    }
    
    arma::mat transVec(std::vector<arma::vec> tmp) {
        const int &N = tmp.size();
        const arma::uword &p = tmp.at(0).n_elem;
        arma::mat result(p, N);
        for(int i = 0; i < N; i++) {
            result.col(i) = tmp.at(i);
        }
        return result;
    }

    
} // namespace SubGuide
