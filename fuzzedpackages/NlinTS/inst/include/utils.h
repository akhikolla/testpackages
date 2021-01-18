#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <algorithm>    //
#include <ctime>        // std::time
#include <random>       // std::default_random_engine, random_shuffle
using namespace std;

typedef vector<double> VectD;
typedef vector <VectD> MatD;

/*****************************/
template<class A, class B>
void shuffle_X_y (vector<A> & X, vector<B> & Y, unsigned seed = 0)
{
    if (seed == 0)
        seed = unsigned (std::time (0));

      // using the same seed for both vectors
      shuffle (X.begin(), X.end(), std::default_random_engine(seed));
      shuffle (Y.begin(), Y.end(), std::default_random_engine(seed));
}


/*****************************/
template<typename T>
void copy_vector (vector<T> & V1, const vector<T> & V2)
{
    V1. clear();
    V1. reserve (V2. size ());
    for (const T & val : V2)
        V1.push_back (val);
}

double get_random();
vector <double> random_vector (unsigned long size);

std::vector<double> get_col (const std::vector<std::vector<double>> & A, unsigned long j);
std::vector<std::vector<double>> Transpose (const std::vector<std::vector<double>> & M);


MatD reshape (const VectD & A, unsigned axis);
void matrix_dot (MatD & A, double a);
// hammart product
VectD matrix_dot (const VectD & A, const VectD & B);
VectD matrix_dot (const MatD & A, const VectD & B);
MatD matrix_dot (const MatD & A, const MatD & B);
VectD matrix_sum (const VectD & A, const VectD & B);
MatD matrix_sum (const MatD & A, const MatD & B);
VectD matrix_sum (const VectD & A, const VectD & B, const VectD & C, const VectD & D);
MatD matrix_sum (const MatD & A, const MatD & B, const MatD & C, const MatD & D);

double sum_vect (const vector<double> & Vect);
double min_vect (const std::vector<double> & Vect);
double max_vect (const std::vector<double> & Vect);
std::vector<std::vector<double>> Normalise (std::vector<std::vector<double>> & mat);
std::vector <double> r_score (const std::vector<std::vector<double> > &pred, const std::vector<std::vector<double> > &real);
//double sum_vect (const std::vector<double> & A);

// activation functions and their derivatives
double sigmoid (double x);
double deriveSigmoid (double x);
double relu (double x);
double deriveRELU (double x);
double tanh_ac (double x);
double deriveTanh (double x);

// activation functions and their derivatives for vectors
VectD sigmoid_v (const VectD & A);
VectD relu_v (const VectD & A);
VectD tanh_v (const VectD & A);
VectD sigmoid_diff (const VectD & A);
VectD relu_diff (const VectD & A);
VectD tanh_diff (const VectD & A);
VectD vect_activation (const VectD & A, const string & f);
VectD diff_activation (const VectD & A, const string & f);

#endif
