/**
 ** Abstract class for layers of a network
 ** Jule 5, 2020
 ** Author: Youssef Hmamouche


  This file is part of NlinTS. NlinTS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
  NlinTS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

 **/

#ifndef LAYER_H
#define LAYER_H

#include<vector>
#include<string>

#include "utils.h"

using namespace std;

class Layer
{
 public:
    //Layer (){}
    virtual ~Layer()=0;
    virtual void set_input_dim (vector<unsigned long> in_dim)=0;
    virtual bool contains_bias ()=0;
    virtual bool is_output ()=0;
    virtual void set_output_layer (bool last)=0;
    virtual vector<vector<double>> get_output ()=0;

    virtual vector<unsigned long> get_output_dim ()=0;
    virtual vector<unsigned long> get_input_dim ()=0;
    //virtual vector<double> simulate (const vector<double> & input, bool store)=0;
    virtual vector<vector<double>> simulate (const vector<vector<double>> & input, bool store)=0;
    virtual void computeErrors(const vector<vector<double>>  & nextErrors)=0;
    virtual void updateWeights (unsigned numb_iter)=0;

    //virtual  void get_errors (vector<double> & A)=0;
    virtual  vector<vector<double>> get_errors ()=0;
    virtual  vector<vector<double>> get_weights ()=0;
    virtual string getType ()=0;

};

#endif // LAYER_H
