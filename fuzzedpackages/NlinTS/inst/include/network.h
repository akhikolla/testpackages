/**
 **  class for a neural network
 ** Jule 5, 2020
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
  */


#ifndef NETWORK_H
#define NETWORK_H

#include<vector>
#include"layer.h"

using namespace std;


class Network
{
	private:
		vector<unsigned long> input_dim;
		vector<Layer*> layers;
		VectD input;
		unsigned nb_layers;

	public:
	  Network (vector<unsigned long> _input_dim);
	  Network();

	  void set_input_size (vector<unsigned long> _input_dim);
	  ~Network();
	  void addLayer (Layer * layer);
	  MatD simulate (const MatD & input, bool store);
	  void backpropagation (const vector <double> & output);
	  void updateWeight (unsigned numb_iter);
	  double loss (const VectD & preds, const VectD & real);
	  double average_loss (const MatD &preds, const MatD &real);
	  void train (const MatD & X, const MatD & y, unsigned numb_iter);
	  //void train (const vector<MatD> &X, const MatD &y);

	  void fit (const MatD & X, const MatD & y, unsigned n_iters, bool shuffle = true);
	  //void fit (const vector<MatD> & X, const MatD & y, int n_iters, bool shuffle = true);

	  MatD predict (const MatD & X);
	  //MatD predict (const vector<MatD> & X);
	  VectD score(const MatD & X, const MatD & y);
	  //VectD score(const vector<MatD> &X, const MatD &y);

	  void summary ();

};

#endif // NETWORK_H
