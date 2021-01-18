// MIT License
// Copyright (c) 2019 Vincent Runge

#include "OmegaSN.h"
#include "Costs.h"


#include <algorithm> // std::reverse // std::min // std::max
#include<iostream>
#include <stdlib.h>

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

OmegaSN::OmegaSN(std::vector< double >& values, double firstdata, unsigned int nbSeg, unsigned int n)
{
  nbStates = values.size();
  nbSegments = nbSeg;
  states = new double[nbStates];
  for(unsigned int i = 0; i < nbStates; i++){states[i] = values[i];}

  ///
  /// MATRIX INITIALIZATION
  ///
  S12P = new double*[3]; ///matrix of vectors S1, S2 and SP
  for(unsigned int i = 0; i < 3; i++){S12P[i] = new double[n + 1];}

  Q = new double**[nbSegments];
  lastChpt = new unsigned int**[nbSegments];
  lastIndState = new unsigned int**[nbSegments];

  for(unsigned int i = 0; i < nbSegments; i++)
  {
    Q[i] = new double*[nbStates];
    lastChpt[i] = new unsigned int*[nbStates];
    lastIndState[i] = new unsigned int*[nbStates];
    for(unsigned int j = 0; j < nbStates; j++)
    {
      Q[i][j] = new double[n + 1];
      lastChpt[i][j] = new unsigned int[n + 1];
      lastIndState[i][j] = new unsigned int[n + 1];
    }
  }
}



OmegaSN::~OmegaSN()
{
  delete[] states;
  states = NULL;
  for(unsigned int i = 0; i < 3; i++){delete[] S12P[i];}
  delete[] S12P;
  S12P = NULL;

  for(unsigned int i = 0; i < nbSegments; i++)
  {
    for(unsigned int j = 0; j < nbStates; j++)
    {
      delete[] Q[i][j];
      delete[] lastChpt[i][j];
      delete[] lastIndState[i][j];
    }
    delete[] Q[i];
    delete[] lastChpt[i];
    delete[] lastIndState[i];
  }
  delete[] Q;
  delete[] lastChpt;
  delete[] lastIndState;
  Q = NULL;
  lastChpt = NULL;
  lastIndState = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > OmegaSN::GetChangepoints() const {return(changepoints);}
std::vector< double > OmegaSN::GetParameters() const {return(parameters);}
double OmegaSN::GetGlobalCost() const {return(globalCost);}
double OmegaSN::GetPruning() const {return(pruning);}


//####### preprocessing #######////####### preprocessing #######////####### preprocessing #######//
//####### preprocessing #######////####### preprocessing #######////####### preprocessing #######//

double** OmegaSN::preprocessing(std::vector< double >& data) const
{
  unsigned int n = data.size();
  S12P[0][0] = 0;
  S12P[1][0] = 0;
  S12P[2][0] = 0;
  for(unsigned int i = 1; i < (n + 1); i++){S12P[0][i] = S12P[0][i-1] + data[i-1];}
  for(unsigned int i = 1; i < (n + 1); i++){S12P[1][i] = S12P[1][i-1] + (data[i-1] * data[i-1]);}
  for(unsigned int i = 1; i < (n + 1); i++){S12P[2][i] = S12P[2][i-1] + i * data[i-1];}
  return(S12P);
}

//####### Q0init #######////####### Q0init #######////####### Q0init #######//
//####### Q0init #######////####### Q0init #######////####### Q0init #######//

void OmegaSN::Q0init(std::vector< double >& data) const
{
  unsigned int n = data.size();
  unsigned int ZERO = 0;
  Costs cost;
  double minCost;
  double temp_cost;
  unsigned int temp_indState = 0;

  /// fill Q[0] = 1 segment only : first 2 columns
  for(unsigned int j = 0; j < nbStates; j++)
  {
    Q[0][j][0] = 0; ///inutile?
    Q[0][j][1] = (data[0] - states[j])*(data[0] - states[j]); ///inutile?
    lastIndState[0][j][0] = j;
    lastIndState[0][j][1] = j;
    lastChpt[0][j][0] = 0;
    lastChpt[0][j][1] = 1;
  }
  /// fill Q[0] = 1 segment only
  for(unsigned int k = 2; k < (n + 1 - nbSegments + 1); k++)
  {
    for(unsigned int j = 0; j < nbStates; j++)
    {
      minCost = INFINITY;
      for(unsigned int i = 0; i < nbStates; i++)
      {
        temp_cost = cost.slopeCost(states[i], states[j], ZERO, k, S12P[0][0], S12P[0][k], S12P[1][0], S12P[1][k], S12P[2][0], S12P[2][k]);
        if(temp_cost < minCost){minCost = temp_cost; temp_indState = i;}
      }
      Q[0][j][k] = minCost;
      lastIndState[0][j][k] = temp_indState;
      lastChpt[0][j][k] = 1;
    }
  }
}

//####### algoNULL #######////####### algoNULL #######////####### algoNULL #######//
//####### algoNULL #######////####### algoNULL #######////####### algoNULL #######//
//no pruning

void OmegaSN::algoNULL(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;
  Costs cost;
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  /// PREPROCESSING
  S12P = preprocessing(data);
  Q0init(data);

  /// ALGO
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    for(unsigned int nb = 1; nb < nbSegments; nb++)
    {
      for(unsigned int v = 0; v < p; v++)
      {
        ///// EXPLORE MATRIX size p*T
        temp_Q = INFINITY;
        temp_indState = 0;
        temp_chpt = 0;

        for(unsigned int t = nb + 1; t < T; t++)
        {
          for(unsigned int u = 0; u < p; u++) /////explore column of states
          {
            temp_cost = Q[nb-1][u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]);
            if(temp_Q > temp_cost)
            {
              temp_Q = temp_cost;
              temp_indState = u;
              temp_chpt = t;
            }
          }
        }
        ///// Write response
        Q[nb][v][T] = temp_Q;
        lastChpt[nb][v][T] = temp_chpt;
        lastIndState[nb][v][T] = temp_indState;
        //std::cout << "nb " << nb << " v "<< v << " T " << T << " Q " << temp_Q << std::endl;
      }
    }
  }
  pruning = 1;

}


//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//no pruning

void OmegaSN::algoISOTONIC(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;
  Costs cost;
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  /// PREPROCESSING
  S12P = preprocessing(data);
  Q0init(data);

  /// ALGO
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    for(unsigned int nb = 1; nb < nbSegments; nb++)
    {
      for(unsigned int v = 0; v < p; v++)
      {
        ///// EXPLORE MATRIX size p*T
        temp_Q = INFINITY;
        temp_indState = 0;
        temp_chpt = 0;

        for(unsigned int t = nb + 1; t < T; t++)
        {
          for(unsigned int u = 0; u < (v + 1); u++) /////explore column of states
          {
            temp_cost = Q[nb-1][u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]);
            if(temp_Q > temp_cost)
            {
              temp_Q = temp_cost;
              temp_indState = u;
              temp_chpt = t;
            }
          }
        }
        ///// Write response
        Q[nb][v][T] = temp_Q;
        lastChpt[nb][v][T] = temp_chpt;
        lastIndState[nb][v][T] = temp_indState;
        //std::cout << "nb " << nb << " v "<< v << " T " << T << " Q " << temp_Q << " --- " << temp_indState << " && " << temp_chpt << std::endl;
      }
    }
  }
  pruning = 1;

}


//####### backtracking #######////####### backtracking #######////####### backtracking #######//
//####### backtracking #######////####### backtracking #######////####### backtracking #######//

void OmegaSN::backtracking(unsigned int n)
{
  unsigned int p = nbStates;
  double temp_Q = Q[nbSegments - 1][0][n];
  unsigned int temp_indState = 0;

  for(unsigned int v = 1; v < p; v++)
  {
    if(Q[nbSegments - 1][v][n] < temp_Q)
    {
      temp_Q = Q[nbSegments - 1][v][n];
      temp_indState = v;
    }
  }

  globalCost = Q[nbSegments - 1][temp_indState][n]; ///fill the globalCost OmegaSN variable
  unsigned int temp_chpt = n;
  unsigned int nb = nbSegments - 1;

// /*
//std::cout << " u " << temp_indState << " t " << temp_chpt << std::endl;

  while(temp_chpt > 1)
  {
    changepoints.push_back(temp_chpt);
    parameters.push_back(states[temp_indState]);

    temp_chpt = lastChpt[nb][temp_indState][temp_chpt];
    temp_indState = lastIndState[nb][temp_indState][changepoints[changepoints.size()-1]];
    nb = nb - 1;
    //std::cout << " u " << temp_indState << " t " << temp_chpt << std::endl;
  }

  changepoints.push_back(1);
  parameters.push_back(states[temp_indState]);

  // */

  ///reverse the vector
  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(parameters.begin(), parameters.end());
}
