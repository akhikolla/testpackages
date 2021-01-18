// MIT License
// Copyright (c) 2019 Vincent Runge

#include<Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include<math.h>
#include"OmegaOP.h"
#include"OmegaSN.h"
#include"op2d.h" //Marco Pascucci code

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List slopeOPtransfer(std::vector<double> data, std::vector<double> states, double penalty, std::string constraint = "null", double minAngle = 0, std::string type = "channel")
{
  OmegaOP omega = OmegaOP(states, data[0], penalty, data.size());
  //DIFFERENT PRUNING + NO CONSTRAINT
  if(type == "null" && constraint == "null"){omega.algo(data);}
  if(type == "channel" && constraint == "null"){omega.algoChannel(data);}
  if(type == "pruning" && constraint == "null"){omega.algoPruning(data);}

  //DIFFERENT CONSTRAINTS
  if(constraint == "isotonic"){omega.algoISOTONIC(data);}
  if(constraint == "unimodal"){omega.algoUNIMODAL(data);}
  if(constraint == "smoothing"){omega.algoSMOOTHING(data, minAngle);}

  omega.backtracking(data.size());

  /// RETURN
  List res = List::create(
    _["changepoints"] = omega.GetChangepoints(),
    _["parameters"] = omega.GetParameters(),
    _["globalCost"] = omega.GetGlobalCost(),
    _["pruningPower"] = omega.GetPruning()
  );
  return res;
}


// [[Rcpp::export]]
List slopeSNtransfer(std::vector<double> data, std::vector<double> states, unsigned int nbSegments, std::string constraint = "null")
{
  OmegaSN omega = OmegaSN(states, data[0], nbSegments, data.size());
  //DIFFERENT CONSTRAINTS
  if(constraint == "null"){omega.algoNULL(data);}
  if(constraint == "isotonic"){omega.algoISOTONIC(data);}

  omega.backtracking(data.size());

  /// RETURN
  List res = List::create(
    _["changepoints"] = omega.GetChangepoints(),
    _["parameters"] = omega.GetParameters(),
    _["globalCost"] = omega.GetGlobalCost(),
    _["pruningPower"] = omega.GetPruning()
  );
  return res;
}


////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// linear fit OP //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// MIT License
// Copyright (c) 2019 Marco Pascucci
//' @title linearOP
//' @description
//' An optimal partitioning algorithm with a linear fit for each segment
//' @name linearOP
//' @param x a vector (see data)
//' @param data a vector defining the data points (x[i], data[i])
//' @param penalty the penalty for introducing a new segment
//' @param cc a boolean to impose a continuity constraint
//'
//' @export
// [[Rcpp::export]]
List linearOP(std::vector<double> x, std::vector<double> data, double penalty, bool cc = false)
{
  if (x.size() != data.size()) {
    stop("x and y must have the same length.");
  }

  PeltResult<double,double> pr = op2D(x,data,penalty);

  for(unsigned int i = 0; i < pr.cp.size(); i++){
    pr.cp[i] = pr.cp[i] + 1;
  }

  // RETURN
  List res = List::create(
    _["cp_indexes"] = pr.cp,
    _["x"] = pr.x,
    _["y"] = pr.y,
    _["globalCost"] = pr.cost
  );

  return(res);
}
