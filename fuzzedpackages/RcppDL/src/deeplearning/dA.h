#ifndef DA_H_
#define DA_H_

#include "utils.h"

class dA {

  public:
    dA(int, int, int , double**, double*, double*);
    dA(int, int, int);
    ~dA();
    void get_corrupted_input(int*, int*, double);
    void get_hidden_values(int*, double*);
    void get_reconstructed_input(double*, double*);
    void train(int*, double, double);
    void reconstruct(int*, double*);

  private:
    int N;
    int n_visible;
    int n_hidden;
    double **W;
    double *hbias;
    double *vbias;

};

#endif
