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

#ifndef VARNN_H
#define VARNN_H

#include <cstdlib>
#include <memory>

#include"struct.h"
#include"network.h"


using namespace std;
using namespace Struct;

class VARNN
{
  private:
      unsigned  numLayers;
      bool bias;
      unsigned  lag;  // lag parameter
      double learning_rate_init; // initial learning rate
      unsigned long  Nb_Ln ;    // numbre of observations
      unsigned long  Nb_Cl ;    // nombre of variables

      string activation; //activation function
      string algo; //backpropagation algorithm

      std::vector<unsigned long> sizeOfLayers;
      std::vector<string> activations;
      //std::shared_ptr <Network> mlp;
      Network mlp;
      std::vector<double> SSR;
      Struct::CMatDouble inputMat ;    // input data

  public:
      VARNN (const vector<unsigned long> & sizeOfLayers, unsigned p, bool bias, double learning_rate_init = 0.1, const  std::vector<string> & activations = {}, const string & algo = "sgd");
      VARNN(){}
     ~VARNN(){}

      /*********************************************************/
      void fit (const CMatDouble & M, unsigned iterations);

      /*********************************************************/
      Struct::CMatDouble forecast (const Struct::CMatDouble & M);

      /*********************************************************/
      void train (const Struct::CMatDouble & M);

      Struct::CVDouble getMSE ();
      Struct::CVDouble getMAE ();

      // Sum of squared errors
      std::vector<double> getSSR ();
};

#endif // VARNN_H
