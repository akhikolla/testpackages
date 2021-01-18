/**
 ** Feed forward layer
 ** Jule 1, 2020
 ** Author: Youssef Hmamouche

   Copyright (c) 2020 Youssef Hmamouche.
   This file is part of NlinTS. NlinTS is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
   NlinTS is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
**/



#ifndef DENSE_H
#define DENSE_H

#include<vector>
#include<string>

#include"layer.h"

using namespace std;

//typedef vector<double> VectD;
//typedef vector <VectD> MatD;

class Dense: public Layer
{
private:

    unsigned long n_neurons;
    string activation;
    double learning_rate_init; // initial learning rate
    unsigned input_dim;
    unsigned bias;
    bool output_layer;
    string algo;

    // last sum of input (without activation function)
     VectD net;

    // last input
    VectD input;
    // output
    VectD O;
    // errors
    VectD E;

    // learning rate of each weight
    MatD alpha;

    // Parameters for the Adam algorithm
    double beta_1; // decay rates for M
    double beta_2; // decay rates for V
    MatD M; //gradient
    MatD V; //squared gradient

    // Weights of the neurons related to the previous layer
    MatD W;

    // Gradients of Weights
    MatD DeltaW;

public:

    Dense (unsigned long n_neurons, string activation = "sigmoid", double learning_rate_init_ = 0.01, bool bias_ = 1,  const string & alg = "sgd");

    void set_input_dim (vector<unsigned long> in_dim);

    bool contains_bias ();
    bool is_output ();
    void set_output_layer (bool last);

    vector<unsigned long> get_output_dim ();
    vector<unsigned long> get_input_dim ();

    MatD simulate (const MatD & input, bool store);


    void computeErrors(const MatD & nextErrors);
    void updateWeights (unsigned numb_iter);
    MatD get_output (){return {O};}
    MatD get_errors ();

    MatD get_weights ();
    string getType ();
};

#endif // LAYER_H
