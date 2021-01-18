#ifndef LR_H_
#define LR_H_


class LogisticRegression {

  public:
    int N;
    int n_in;
    int n_out;
    double **W;
    double *b;
    LogisticRegression(int, int, int);
    ~LogisticRegression();
    void train(int*, int*, double);
    void softmax(double*);
    void predict(int*, double*);
};

#endif
