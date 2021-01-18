#ifndef RCPPDL_H_
#define RCPPDL_H_

#include <vector>

#include "util.h"
#include "deeplearning/SdA.h"
#include "deeplearning/DBN.h"

class RcppDA {
  public:
    RcppDA();
    void init(SEXP x);
    void train();
    Rcpp::NumericMatrix reconstruct(SEXP test);
    Rcpp::List show();
    void setlr(double lr) {
        learning_rate = lr;
    };
    void setcl(double cl) {
        corruption_level = cl;
    };
    void setTE(int te) {
        training_epochs = te;
    };
    void setHidden(int h) {
        n_hidden = h;
    };

  private:
    dA * da;
    int ** train_X;
    double learning_rate;
    double corruption_level;
    int training_epochs;
    int train_N;
    int test_N;
    int n_visible;
    int n_hidden;

};


class RcppSDA {

  public:
    RcppSDA();
    void init(SEXP x, SEXP y, SEXP hidden);
    void pretrain();
    void finetune();
    Rcpp::NumericMatrix predict(SEXP test);
    Rcpp::List show();

    void setPlr(double plr) {
        pretrain_lr = plr;
    };

    void setFlr(double flr) {
        finetune_lr = flr;
    };

    void setcl(double cl) {
        corruption_level = cl;
    };

    void setPE(int pe) {
        pretraining_epochs = pe;
    };

    void setFE(int fe) {
        finetune_epochs = fe;
    };

  private:
    SdA * sda;
    int ** train_X;
    int ** train_Y;
    double pretrain_lr;
    double corruption_level;
    int pretraining_epochs;
    double finetune_lr;
    int finetune_epochs;
    int train_N;
    int n_ins;
    int n_outs;
    std::vector<int> hidden_layer_sizes;

};

class RcppRBM {

  public:
    RcppRBM();
    void init(SEXP x);
    void train();
    Rcpp::NumericMatrix reconstruct(SEXP test);
    Rcpp::List show();
    void setlr(double lr) {
        learning_rate = lr;
    };
    void setTE(int t) {
        training_epochs = t;
    };
    void setStep(int t) {
        step = t;
    };
    void setHidden(int h) {
        n_hidden = h;
    };

  private:
    RBM * rbm;
    double learning_rate;
    int training_epochs;
    int step;
    int n_hidden;
    int n_visible;
    int ** train_X;
    int train_N;
};

class RcppDBN {

  public:
    RcppDBN();
    void init(SEXP x, SEXP y, SEXP hidden);
    void pretrain();
    void finetune();
    Rcpp::NumericMatrix predict(SEXP test);
    Rcpp::List show();
    void setPlr(double plr) {
        pretrain_lr = plr;
    };
    void setPE(int p) {
        pretraining_epochs = p;
    };
    void setFlr(double flr) {
        flr = finetune_lr;
    };
    void setFE(int t) {
        finetune_epochs = t;
    };
    void setStep(int t) {
        step = t;
    };

  private:
    DBN * dbn;
    int ** train_X;
    int ** train_Y;

    double pretrain_lr;
    int pretraining_epochs;
    double finetune_lr;
    int finetune_epochs;

    int step;

    int train_N;
    int n_ins;
    int n_outs;
    std::vector<int> hidden_layer_sizes;
};

#endif
