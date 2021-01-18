// MIT License
// Copyright (c) 2019 Vincent Runge

#include "OmegaOP.h"
#include "Costs.h"


#include <algorithm> // std::reverse // std::min // std::max
#include<iostream>
#include <stdlib.h>

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

OmegaOP::OmegaOP(std::vector< double >& values, double firstdata, double beta, unsigned int n)
{
  nbStates = values.size();
  states = new double[nbStates];
  for(unsigned int i = 0; i < nbStates; i++){states[i] = values[i];}

  /// MATRIX INITIALIZATION
  S12P = new double*[3]; ///matrix of vectors S1, S2 and SP
  Q = new double*[nbStates]; ///matrix of costs
  lastChpt = new unsigned int*[nbStates]; ///matrix of best last changepoints
  lastIndState = new unsigned int*[nbStates]; ///matrix of starting states for the best last segment

  for(unsigned int i = 0; i < 3; i++){S12P[i] = new double[n + 1];}
  for(unsigned int i = 0; i < nbStates; i++){Q[i] = new double[n + 1];}
  for(unsigned int i = 0; i < nbStates; i++){lastChpt[i] = new unsigned int[n + 1];}
  for(unsigned int i = 0; i < nbStates; i++){lastIndState[i] = new unsigned int[n + 1];}

  /// FILL FIRST COLUMN in Q, lastIndState,lastChpt
  for(unsigned int i = 0; i < nbStates; i++)
  {
    Q[i][0] = 0;
    Q[i][1] = (firstdata - states[i])*(firstdata - states[i]);
    lastIndState[i][0] = i;
    lastIndState[i][1] = i;
    lastChpt[i][0] = 0;
    lastChpt[i][1] = 1;
  }
  penalty = beta;
}



OmegaOP::~OmegaOP()
{
  delete[] states;
  states = NULL;
  for(unsigned int i = 0; i < 3; i++){delete[] S12P[i];}
  for(unsigned int i = 0; i < nbStates; i++){delete[] Q[i];}
  for(unsigned int i = 0; i < nbStates; i++){delete[] lastChpt[i];}
  for(unsigned int i = 0; i < nbStates; i++){delete[] lastIndState[i];}
  delete[] S12P;
  S12P = NULL;
  delete[] Q;
  Q = NULL;
  delete[] lastChpt;
  lastChpt = NULL;
  delete[] lastIndState;
  lastIndState = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > OmegaOP::GetChangepoints() const {return(changepoints);}
std::vector< double > OmegaOP::GetParameters() const {return(parameters);}
double OmegaOP::GetGlobalCost() const {return(globalCost);}
double OmegaOP::GetPruning() const {return(pruning);}


//####### preprocessing #######////####### preprocessing #######////####### preprocessing #######//
//####### preprocessing #######////####### preprocessing #######////####### preprocessing #######//

double** OmegaOP::preprocessing(std::vector< double >& data) const
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


// algo // algoChannel // algoPruning // algoPruningPELT
// backtracking
// algoISOTONIC // algoUNIMODAL // algoOUTLIER

//####### algo #######////####### algo #######////####### algo #######//
//####### algo #######////####### algo #######////####### algo #######//

void OmegaOP::algo(std::vector< double >& data)
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

  ///
  /// ALGO
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      ///// EXPLORE MATRIX size p*T
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 1; t < T; t++)
      {
        for(unsigned int u = 0; u < p; u++) /////explore column of states
        {
          temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]) + penalty;
          //std::cout <<" t2 "<<  t << " v2 " << states[u] << " t3 "<< T <<" v3 " << states[v] << " cost  " << temp_cost;
          if(temp_Q > temp_cost)
          {
            temp_Q = temp_cost;
            temp_indState = u;
            temp_chpt = t;
          }
          //std::cout << std::endl;
        }

      }
      ///// Write response
      Q[v][T] = temp_Q;
      lastChpt[v][T] = temp_chpt;
      lastIndState[v][T] = temp_indState;
    }
  }
  pruning = 1;

}


//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//
//####### algoChannel #######////####### algoChannel #######////####### algoChannel #######//

void OmegaOP::algoChannel(std::vector< double >& data)
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

  //////////////////////////////////////////////////
  /// CHANNEL INFORMATION
  /// u1 / u2 = "min / max" in each column of Q
  unsigned int* u1 = new unsigned int[n + 1];
  unsigned int* u2 = new unsigned int[n + 1];
  unsigned int theStart;
  unsigned int theEnd;
  double theV = 0;
  unsigned int indexTheV = 0;
  unsigned int nbPosition = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    /// FILL u1 and u2 vectors in position T-1
    theStart = 0;
    while(theStart < (p - 1) && Q[theStart][T - 1] > Q[theStart + 1][T - 1])
      {theStart = theStart + 1;}
    u1[T-1] = theStart;
    theEnd = p - 1;
    while(theEnd > 0 && Q[theEnd][T - 1] > Q[theEnd - 1][T - 1])
      {theEnd = theEnd - 1;}
    u2[T-1] = theEnd;

    for(unsigned int v = 0; v < p; v++)
    {
      ///// EXPLORE MATRIX size p*T
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 1; t < T; t++)
      {
        ///// FIND the minimum of the cost in start state
        if(t < (T-1))
        {
          theV = cost.vhat(states[v], t, T, S12P[0][t], S12P[0][T], S12P[2][t], S12P[2][T]);
          indexTheV = cost.closestStateIndex(theV, states, p);
        }
        else
        {
          indexTheV = u1[T-1]; // if t = T-1 cost.slopeCost does not depend on u
        }

        /// explore values (in column of states) ONLY between min(u1[t],indexTheV) and max(u2[t],indexTheV)
        for(unsigned int u = std::min(u1[t],indexTheV); u < std::max(u2[t],indexTheV) + 1; u++)
        {
          nbPosition = nbPosition + 1; //we explore +1 position
          temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]) + penalty;
          //std::cout <<" t2 "<<  t << " v2 " << states[u] << " t3 "<< T <<" v3 " << states[v] << " cost  " << temp_cost;
          if(temp_Q > temp_cost)
          {
            temp_Q = temp_cost;
            temp_indState = u;
            temp_chpt = t;
          }
          //std::cout << std::endl;
        }
      }

      ///// Write response
      Q[v][T] = temp_Q;
      lastChpt[v][T] = temp_chpt;
      lastIndState[v][T] = temp_indState;
    }
  }

  pruning = 2.0*nbPosition/(1.0*p*p*n*(n-1)); //nbPosition seen / nbPosition total

  delete [] u1;
  u1 = NULL;
  delete [] u2;
  u2 = NULL;
}


//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//
//####### algoPruning #######////####### algoPruning #######////####### algoPruning #######//

void OmegaOP::algoPruning(std::vector< double >& data)
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

  /// list for positions (u,t)
  std::list< unsigned int>* u_pos = new std::list< unsigned int>[p];
  std::list< unsigned int>* t_pos = new std::list< unsigned int>[p];
  std::list<unsigned int>::iterator u_it;
  std::list<unsigned int>::iterator t_it;

  //////////////////////////////////////////////////
  /// PRUNING INFORMATION
  ///
  double S;
  double** coeffsAG = new double*[n + 1]; ///matrix of coefficients: Ap, Am, Gp, Gm
  for(unsigned int i = 0; i < (n + 1); i++){coeffsAG[i] = new double[4];}
  cost.fillCoeffsAG(coeffsAG, S12P[0], n);
  unsigned int nbPosition = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      /////
      ///// Add last column to explore
      /////
      for(unsigned int w = 0; w < p; w++)
      {
        u_pos[v].push_back(w);
        t_pos[v].push_back(T-1);
      }

      ///// EXPLORE ""MATRIX"" size p*T (this is a list of positions (u,t))
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      u_it = u_pos[v].begin();
      t_it = t_pos[v].begin();
      while(t_it != t_pos[v].end())
      {
        nbPosition = nbPosition + 1;
        temp_cost = Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S12P[0][*t_it], S12P[0][T], S12P[1][*t_it], S12P[1][T], S12P[2][*t_it], S12P[2][T]) + penalty;
        if(temp_Q > temp_cost)
        {
          temp_Q = temp_cost;
          temp_indState = *u_it;
          temp_chpt = *t_it;
        }
        ++u_it;
        ++t_it;
      }

      /////
      ///// Write response
      /////
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;

      /////
      ///// PRUNING STEP
      /////
      u_it = u_pos[v].begin();
      t_it = t_pos[v].begin();
      while(t_it != t_pos[v].end())
      {
        ///inequality test
        if(Q[*u_it][*t_it] + cost.slopeCost(states[*u_it], states[v], *t_it, T, S12P[0][*t_it], S12P[0][T], S12P[1][*t_it], S12P[1][T], S12P[2][*t_it], S12P[2][T]) > Q[v][T])
        {
          ///linear regression
          S = S12P[2][T] - S12P[2][*t_it] - *t_it*(S12P[0][T] - S12P[0][*t_it]);
          if(cost.pruningTest(states[*u_it], states[v], *t_it, T, n, S, coeffsAG[T][0], coeffsAG[T][1], coeffsAG[T][2], coeffsAG[T][3]))
          {
            u_it = u_pos[v].erase(u_it);
            t_it = t_pos[v].erase(t_it);
          }
          else{++u_it; ++t_it;}
        }
        else{++u_it; ++t_it;}
      }

    }
  }

  pruning = 2.0*nbPosition/(1.0*p*p*n*(n-1)); //nbPosition seen / nbPosition total

  for(unsigned int i = 0; i < (n + 1); i++){delete(coeffsAG[i]);}
  delete [] coeffsAG;
  coeffsAG = NULL;
  delete [] t_pos;
  t_pos = NULL;
  delete [] u_pos;
  u_pos = NULL;
}

//####### backtracking #######////####### backtracking #######////####### backtracking #######//
//####### backtracking #######////####### backtracking #######////####### backtracking #######//

void OmegaOP::backtracking(unsigned int n)
{
  unsigned int p = nbStates;
  double temp_Q = Q[0][n];
  unsigned int temp_indState = 0;

  for(unsigned int v = 1; v < p; v++)
  {
    if(Q[v][n] < temp_Q)
    {
      temp_Q = Q[v][n];
      temp_indState = v;
    }
  }

  globalCost = Q[temp_indState][n]; ///fill the globalCost OmegaOP variable
  unsigned int temp_chpt = n;

  while(temp_chpt > 1)
  {
    changepoints.push_back(temp_chpt);
    parameters.push_back(states[temp_indState]);

    temp_chpt = lastChpt[temp_indState][temp_chpt];
    temp_indState = lastIndState[temp_indState][changepoints[changepoints.size()-1]];
  }

  changepoints.push_back(1);
  parameters.push_back(states[temp_indState]);

  ///reverse the vector
  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(parameters.begin(), parameters.end());
}

//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//####### algoISOTONIC #######////####### algoISOTONIC #######////####### algoISOTONIC #######//
//channel pruning

void OmegaOP::algoISOTONIC(std::vector< double >& data)
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

  //////////////////////////////////////////////////
  /// CHANNEL INFORMATION
  /// u1 / u2 = "min / max" in each column of Q
  unsigned int* u1 = new unsigned int[n + 1];
  unsigned int* u2 = new unsigned int[n + 1];
  unsigned int theStart;
  unsigned int theEnd;
  double theV = 0;
  unsigned int indexTheV = 0;
  unsigned int nbPosition = 0;

  ///
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    /// FILL u1 and u2 vectors in position T-1
    theStart = 0;
    while(theStart < (p - 1) && Q[theStart][T - 1] > Q[theStart + 1][T - 1])
    {theStart = theStart + 1;}
    u1[T-1] = theStart;
    theEnd = p - 1;
    while(theEnd > 0 && Q[theEnd][T - 1] > Q[theEnd - 1][T - 1])
    {theEnd = theEnd - 1;}
    u2[T-1] = theEnd;

    for(unsigned int v = 0; v < p; v++)
    {
      ///// EXPLORE MATRIX size p*T
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 1; t < T; t++)
      {
        ///// FIND the minimum of the cost in start state
        if(t < (T-1))
        {
          theV = cost.vhat(states[v], t, T, S12P[0][t], S12P[0][T], S12P[2][t], S12P[2][T]);
          indexTheV = cost.closestStateIndex(theV, states, p);
        }
        else
        {
          indexTheV = u1[T-1]; // if t = T-1 cost.slopeCost does not depend on u
        }

        /// explore values (in column of states) with add. constaint u <= v
        for(unsigned int u = std::min(std::min(u1[t],indexTheV), v); u < std::min(std::max(u2[t],indexTheV), v) + 1; u++) /////explore column of states
        {
          nbPosition = nbPosition + 1; //we explore +1 position
          temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]) + penalty;
          //std::cout <<" t2 "<<  t << " v2 " << states[u] << " t3 "<< T <<" v3 " << states[v] << " cost  " << temp_cost;
          if(temp_Q > temp_cost)
          {
            temp_Q = temp_cost;
            temp_indState = u;
            temp_chpt = t;
          }
          //std::cout << std::endl;
        }
      }

      ///// Write response
      Q[v][T] = temp_Q;
      lastChpt[v][T] = temp_chpt;
      lastIndState[v][T] = temp_indState;
    }
  }

  pruning = 2.0*nbPosition/(1.0*p*p*n*(n-1)); //nbPosition seen / nbPosition total

  delete [] u1;
  u1 = NULL;
  delete [] u2;
  u2 = NULL;
}

//####### algoUNIMODAL #######////####### algoUNIMODAL #######////####### algoUNIMODAL #######//
//####### algoUNIMODAL #######////####### algoUNIMODAL #######////####### algoUNIMODAL #######//
// NO PRUNING

void OmegaOP::algoUNIMODAL(std::vector< double >& data)
{
  unsigned int n = data.size();
  unsigned int p = nbStates;
  Costs cost;
  double temp_cost = 0;
  double temp_Q = -1;
  int temp_chpt = -1;
  unsigned int temp_indState = 0;

  int** SLOPE = new int*[p];
  for(unsigned int i = 0; i < p; i++){SLOPE[i] = new int[n];}

  /// PREPROCESSING
  S12P = preprocessing(data);

  ///
  /// ALGO
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      ///// EXPLORE MATRIX size p*T
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 1; t < T; t++)
      {
        for(unsigned int u = 0; u < p; u++) /////explore column of states
        {
          if(!(u < v && SLOPE[u][t] == -1))
          {
            temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]) + penalty;
            if(temp_Q > temp_cost)
            {
              temp_Q = temp_cost;
              temp_indState = u;
              temp_chpt = t;
            }
          }
        }

      }
      /////
      ///// Write response
      /////
      Q[v][T] = temp_Q;
      lastIndState[v][T] = temp_indState;
      lastChpt[v][T] = temp_chpt;

      if(temp_indState > v){SLOPE[v][T] = -1;}
      if(temp_indState < v){SLOPE[v][T] = 1;}
      if(temp_indState == v){if(SLOPE[temp_indState][temp_chpt] == -1){SLOPE[v][T] = -1;}}
    }
  }

  pruning = 1;

  for(unsigned int i = 0; i < p; i++){delete(SLOPE[i]);}
  delete [] SLOPE;
  SLOPE = NULL;

}



//####### algoSMOOTHING #######////####### algoSMOOTHING #######////####### algoSMOOTHING #######//
//####### algoSMOOTHING #######////####### algoSMOOTHING #######////####### algoSMOOTHING #######//
// NO PRUNING

void OmegaOP::algoSMOOTHING(std::vector< double >& data, double minAngle)
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

  ///
  /// ALGO
  /// states u to v -> time position t to T
  /// explore in (u,t) for fixed (v,T)
  ///
  for(unsigned int T = 2; T < (n + 1); T++)
  {
    for(unsigned int v = 0; v < p; v++)
    {
      ///// EXPLORE MATRIX size p*T
      temp_Q = INFINITY;
      temp_indState = 0;
      temp_chpt = 0;

      for(unsigned int t = 1; t < T; t++)
      {
        for(unsigned int u = 0; u < p; u++) /////explore column of states
        {
          if(cost.angleTest(lastChpt[u][t], t, T, states[lastIndState[u][t]], states[u], states[v], minAngle))
          {
            temp_cost = Q[u][t] + cost.slopeCost(states[u], states[v], t, T, S12P[0][t], S12P[0][T], S12P[1][t], S12P[1][T], S12P[2][t], S12P[2][T]) + penalty;
            //std::cout <<" t2 "<<  t << " v2 " << states[u] << " t3 "<< T <<" v3 " << states[v] << " cost  " << temp_cost;
            if(temp_Q > temp_cost)
            {
              temp_Q = temp_cost;
              temp_indState = u;
              temp_chpt = t;
            }
          }
          //std::cout << std::endl;
        }

      }
      ///// Write response
      Q[v][T] = temp_Q;
      lastChpt[v][T] = temp_chpt;
      lastIndState[v][T] = temp_indState;
    }
  }
  pruning = 1;

}
