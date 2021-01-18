#ifndef DBN_H_
#define DBN_H_

#include <iostream>
#include <cstdlib>

#include "utils.h"
#include "HiddenLayer.h"
#include "RBM.h"
#include "LogisticRegression.h"

class DBN {

  public:
    int N;
    int n_ins;
    int *hidden_layer_sizes;
    int n_outs;
    int n_layers;
    HiddenLayer **sigmoid_layers;
    RBM **rbm_layers;
    LogisticRegression *log_layer;
    DBN(int, int, int*, int, int);
    ~DBN();
    void pretrain(int**, double, int, int);
    void finetune(int**, int**, double, int);
    void predict(int*, double*);
};

#endif
