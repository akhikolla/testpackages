/*
  This file is part of NlinTS. NlinTS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
  NlinTS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#include<Rcpp.h>
#include"../inst/include/matrixOperations.h"
#include"../inst/include/varnn.h"
#include"../inst/include/dense.h"

using namespace MatrixOperations;

VARNN::VARNN(const vector<unsigned long> &sizeOfLayers, unsigned p, bool bias, double learning_rate_init /*= 0.1*/, const  std::vector<string> & activations /*= {}*/, const string & algo /*= "sgd"*/):
    bias (bias),
    lag (p),
    learning_rate_init (learning_rate_init),
    algo (algo),
    sizeOfLayers (sizeOfLayers),
    activations (activations)
{
    numLayers = unsigned (sizeOfLayers.size ()) + 1;
    //this->mlp = shared_ptr <Network> (new Network ());
    mlp = Network ();
    //this->mlp = new Network ();
}

void VARNN::fit(const CMatDouble &M, unsigned iterations)
{
    CMatDouble A, Real;
    CVDouble result, res_out;
    CMatDouble minMax;

    MatD X, Y, Predictions;

    // Input matrix
    inputMat = M;

    // Input dimensions
    Nb_Ln = M[0].size ();
    Nb_Cl = M.size ();

    // Normalisation
    minMax = inputMat.Normalise ();

    // activations function for eac layer
    if (activations. size () != numLayers)
    {
        activations. clear ();
        // relu for hidden layers
        for (unsigned i = 0; i < sizeOfLayers. size (); ++i)
            activations. push_back("relu");

        // sigmoid for output layer
         activations. push_back("sigmoid");
    }

    // Construct the matrix of lagged variables (VAR (p) representation)
    for (auto & vec:inputMat)
        P_Part (vec, Real, A, lag);

    X = A.to_Mat();
    Y = Real.to_Mat();

    // Define the network structure
    vector<unsigned long> input_dim;
    input_dim. push_back (X[0].size ());

    // input dimension of the network
    this->mlp.set_input_size (input_dim);

    // hidden layers
    for (unsigned i = 0; i < sizeOfLayers. size (); ++i)
    {
        Dense * hidenLayer = new Dense (sizeOfLayers[i], activations[i], learning_rate_init, bias, algo);
        mlp.addLayer (hidenLayer);
    }

    // output layer
    Dense * outputLayer = new Dense (Real.size (), activations. back (), learning_rate_init, bias, algo);
    mlp.addLayer (outputLayer);

    //Summary
    //mlp.summary ();

    // train the network
    mlp.fit (X, Y, iterations, true);

    // compute the prediction and transforming them inti CMAtDouble class type
    Predictions = mlp.predict (X);

    // Evaluate the model on training data
    SSR.clear ();
    SSR.resize (Y[0].size (), 0);

    for (unsigned long i = 0 ; i < Y[0].size () ; ++i)
    {
        for (unsigned j = 0 ; j < Y.size() ; j++)
            SSR [i] += pow (Predictions[j][i] - Y[j][i], 2);
    }
}

CMatDouble VARNN::forecast(const CMatDouble &M)
{
    CMatDouble A, Predictions;

    // we can use the output matrix temporarely
    Predictions = M;

    CMatDouble minMax;

    minMax = Predictions.Normalise ();

    if (Predictions.size () == Nb_Cl and Predictions[0].size () >= lag)
    {
        for (auto vec:Predictions)
            Pr_Part (vec, A, lag);

        Predictions.clear ();
        Predictions.Init_Mat (mlp.predict (A.to_Mat()));
    }
    else
    {
        Rcpp::Rcout << "Error in input dimensions.\n";
        Rcpp::Rcout << "The model expects a matrix of " << Nb_Cl << " columns. \n";
        Rcpp::Rcout << "While the input matrix conatins " << Predictions.size () << " columns. \n";
        Rcpp::stop ("\n.");
    }

    Predictions.Denormalising (minMax);

    return Predictions;
}

void VARNN::train(const CMatDouble &M)
{
    CMatDouble B, present, Res, minMax;

    Res = M;

    minMax = Res.Normalise();

    if (Res.size () == Nb_Cl and Res[0].size () >= lag)
    {
        for (auto vec:Res)
            P_Part (vec, present, B, lag);

        // Train the model
        mlp.train (B.to_Mat(), present.to_Mat(), 0);
    }
    else
    {
        Rcpp::Rcout << "Error in input dimensions.\n";
        Rcpp::Rcout << "The model expects a matrix of " << Nb_Cl << " columns. \n";
        Rcpp::Rcout << "While the input matrix conatins " << Res.size () << " columns.";
        Rcpp::stop ("\n.");
    }
}

std::vector<double> VARNN::getSSR() {return SSR;}
