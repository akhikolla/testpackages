/*
 *  Rcpp.cpp
 *  PCMBaseCpp
 *
 * Copyright 2017,2018 Venelin Mitov
 *
 * This file is part of PCMBaseCpp: A C++ backend for calculating the likelihood of phylogenetic comparative models.
 *
 * PCMBaseCpp is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU General Public License as
 * published by the Free Software Foundation.
 *
 * PCMBaseCpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with PCMBaseCpp.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */
#include <RcppArmadillo.h>

#include<vector>
#include<string>
#include<sstream>

#include "QuadraticPolyWhite.h"
#include "QuadraticPolyBM.h"
#include "QuadraticPolyBM1D.h"
#include "QuadraticPolyOU.h"
#include "QuadraticPolyOU1D.h"
#include "QuadraticPolyJOU.h"
#include "QuadraticPolyDOU.h"
#include "QuadraticPolyMixedGaussian.h"
#include "QuadraticPolyMixedGaussian1D.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

void R_init_PCMBaseCpp(DllInfo *info) {
   /* Register routines, allocate resources. */
   R_registerRoutines(info, NULL, NULL, NULL, NULL);
   R_useDynamicSymbols(info, TRUE);
}

void R_unload_PCMBaseCpp(DllInfo *info) {
   /* Release resources. */
}

using namespace PCMBaseCpp;
using namespace std;

SPLITT::Tree<uint, double>* CreatePCMBaseCppTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::Tree<uint, double>(br_0, br_1, t);
}

RCPP_MODULE(PCMBaseCpp__Tree) {
  Rcpp::class_<SPLITT::Tree<uint, double> > ( "PCMBaseCpp__Tree" )
  .factory<Rcpp::List const&>( &CreatePCMBaseCppTree )
  .property("num_nodes", &SPLITT::Tree<uint, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &SPLITT::Tree<uint, double>::FindChildren )
  .method("OrderNodes", &SPLITT::Tree<uint, double>::OrderNodes )
  ;
}

SPLITT::OrderedTree<uint, double>* CreatePCMBaseCppOrderedTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::OrderedTree<uint, double>(br_0, br_1, t);
}


RCPP_MODULE(PCMBaseCpp__OrderedTree) {
  Rcpp::class_<SPLITT::Tree<uint, double> > ( "PCMBaseCpp__Tree" )
  .factory<Rcpp::List const&>( &CreatePCMBaseCppTree )
  .property("num_nodes", &SPLITT::Tree<uint, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<uint, double>::FindIdOfParent )
  .method("FindChildren", &SPLITT::Tree<uint, double>::FindChildren )
  .method("OrderNodes", &SPLITT::Tree<uint, double>::OrderNodes )
  ;
  Rcpp::class_<SPLITT::OrderedTree<uint, double> >( "PCMBaseCpp__OrderedTree" )
    .derives<SPLITT::Tree<uint, double> > ( "PCMBaseCpp__Tree" )
    .factory<Rcpp::List const&>( &CreatePCMBaseCppOrderedTree )
    .property("num_levels", &SPLITT::OrderedTree<uint, double>::num_levels )
    .property("num_parallel_ranges_prune", &SPLITT::OrderedTree<uint, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &SPLITT::OrderedTree<uint, double>::ranges_id_visit )
    .property("ranges_id_prune", &SPLITT::OrderedTree<uint, double>::ranges_id_prune )
  ;
}

struct ParsedRObjects {
  double threshold_SV;
  double threshold_EV;
  double threshold_skip_singular;
  double threshold_Lambda_ij;
  double NA_double_;
  
  bool skip_singular;
  bool transpose_Sigma_x;
  arma::mat const& X;
  arma::cube VE;
  Rcpp::List pcListInt;
  std::vector<arma::uvec> Pc;
  SPLITT::uvec br_0;
  SPLITT::uvec br_1; 
  SPLITT::vec t;
  
  uint RModel;
  std::vector<arma::uword> regimes;
  std::vector<arma::u8> jumps;
  
  SPLITT::uint num_tips;
  SPLITT::uint num_branches;
  SPLITT::uvec tip_names;
  
  ParsedRObjects(
    arma::mat const& X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo):
    
    threshold_SV(static_cast<double>(metaInfo["PCMBase.Threshold.SV"])),
    threshold_EV(static_cast<double>(metaInfo["PCMBase.Threshold.EV"])),
    threshold_skip_singular(static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"])), 
    threshold_Lambda_ij(static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"])),
    NA_double_(static_cast<double>(metaInfo["NA_double_"])),
    skip_singular(static_cast<int>(metaInfo["PCMBase.Skip.Singular"])),
    transpose_Sigma_x(static_cast<int>(metaInfo["PCMBase.Transpose.Sigma_x"])),
    X(X),
    VE(Rcpp::as<arma::cube>(metaInfo["VE"])),
    pcListInt(Rcpp::as<Rcpp::List>(metaInfo["pcListInt"])), 
    Pc(Rcpp::as<arma::uword>(metaInfo["M"])),
    RModel(Rcpp::as<uint>(metaInfo["RModel"])),
    regimes(Rcpp::as<vector<arma::uword> >(metaInfo["r"])),
    jumps(Rcpp::as<vector<arma::u8> >(metaInfo["xi"])), 
    num_tips(Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size()),
    tip_names(SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips)) {
    
    if(threshold_SV < 0) {
      ostringstream os;
      os<<"Rcpp.cpp:ParsedRObjects:: The argument threshold_SV should be non-negative real number.";
      throw invalid_argument(os.str());
    }
    if(threshold_EV < 0) {
      ostringstream os;
      os<<"Rcpp.cpp:ParsedRObjects:: The argument threshold_EV should be non-negative real number.";
      throw invalid_argument(os.str());
    }
    if(threshold_Lambda_ij < 0) {
      ostringstream os;
      os<<"ERR:03825:PCMBaseCpp:Rcpp.cpp:ParsedRObjects:: The argument threshold_Lambda_ij should be non-negative double.";
      throw invalid_argument(os.str());
    }
    
    for(arma::uword i = 0; i < Pc.size(); ++i) {
      Pc[i] = Rcpp::as<arma::uvec>(pcListInt[i]);
    }
    
    arma::umat branches = tree["edge"];
    br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
    br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
    t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
    num_branches = branches.n_rows;
    
    using namespace std;
    
    if(regimes.size() != branches.n_rows) {
      ostringstream os;
      os<<"ERR:03821:PCMBaseCpp:Rcpp.cpp:ParsedRObjects:: The slot r in metaInfo has different length ("<<regimes.size()<<
        ") than the number of edges ("<<branches.n_rows<<").";
      throw logic_error(os.str());
    }
    
    if(jumps.size() != branches.n_rows) {
      ostringstream os;
      os<<"ERR:03822:PCMBaseCpp:Rcpp.cpp:ParsedRObjects:: The slot jumps in trees has different length ("<<jumps.size()<<
        ") than the number of edges ("<<branches.n_rows<<").";
      throw logic_error(os.str());
    }
  }
};

QuadraticPolyWhite* CreateQuadraticPolyWhite(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyWhite::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
  }
  
  typename QuadraticPolyWhite::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, pObjs.Pc, pObjs.RModel, 
      std::vector<std::string>(), 
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyWhite(pObjs.br_0, pObjs.br_1, lengths, data);
}

  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyWhite::TreeType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyWhite::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolyWhite) {
    Rcpp::class_<QuadraticPolyWhite::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyWhite_Tree" )
    .property("num_nodes", &QuadraticPolyWhite::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolyWhite::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolyWhite::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolyWhite::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolyWhite::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolyWhite::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolyWhite::TreeType>( "PCMBaseCpp__QuadraticPolyWhite_OrderedTree" )
      .derives<QuadraticPolyWhite::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyWhite_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolyWhite::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolyWhite::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolyWhite::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolyWhite::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolyWhite::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolyWhite::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyWhite_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolyWhite::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolyWhite::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolyWhite::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyWhite_ParallelPruning" )
      .derives<QuadraticPolyWhite::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyWhite_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolyWhite::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolyWhite::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolyWhite::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolyWhite::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolyWhite::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolyWhite::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolyWhite::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolyWhite>( "PCMBaseCpp__QuadraticPolyWhite" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyWhite)
      .method( "TraverseTree", &QuadraticPolyWhite::TraverseTree )
      .method( "StateAtNode", &QuadraticPolyWhite::StateAtNode )
      .property( "tree", &QuadraticPolyWhite::tree )
      .property( "algorithm", &QuadraticPolyWhite::algorithm )
    ;
  }

QuadraticPolyBM* CreateQuadraticPolyBM(
    arma::mat const& X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
    
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyBM::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
  }
  
  typename QuadraticPolyBM::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, pObjs.Pc, pObjs.RModel, 
      std::vector<std::string>(), 
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyBM(pObjs.br_0, pObjs.br_1, lengths, data);
}

//RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyBM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyBM::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyBM) {
  Rcpp::class_<QuadraticPolyBM::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyBM_Tree" )
  .property("num_nodes", &QuadraticPolyBM::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyBM::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyBM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyBM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyBM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyBM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyBM::TreeType>( "PCMBaseCpp__QuadraticPolyBM_OrderedTree" )
    .derives<QuadraticPolyBM::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyBM_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyBM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyBM::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyBM::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyBM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyBM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyBM::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyBM_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyBM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyBM::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyBM::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyBM_ParallelPruning" )
    .derives<QuadraticPolyBM::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyBM_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyBM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyBM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyBM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyBM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyBM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyBM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyBM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyBM>( "PCMBaseCpp__QuadraticPolyBM" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyBM)
    .method( "TraverseTree", &QuadraticPolyBM::TraverseTree )
    .method( "StateAtNode", &QuadraticPolyBM::StateAtNode )
    .property( "tree", &QuadraticPolyBM::tree )
    .property( "algorithm", &QuadraticPolyBM::algorithm )
  ;
}


QuadraticPolyBM1D* CreateQuadraticPolyBM1D(
    arma::mat const& X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyBM1D::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
  }
  
  typename QuadraticPolyBM1D::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, 
      pObjs.RModel, 
      std::vector<std::string>(), 
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyBM1D(pObjs.br_0, pObjs.br_1, lengths, data);
}

//RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyBM1D::TreeType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyBM1D::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolyBM1D) {
    Rcpp::class_<QuadraticPolyBM1D::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyBM1D_Tree" )
    .property("num_nodes", &QuadraticPolyBM1D::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolyBM1D::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolyBM1D::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolyBM1D::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolyBM1D::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolyBM1D::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolyBM1D::TreeType>( "PCMBaseCpp__QuadraticPolyBM1D_OrderedTree" )
      .derives<QuadraticPolyBM1D::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyBM1D_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolyBM1D::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolyBM1D::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolyBM1D::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolyBM1D::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolyBM1D::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolyBM1D::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyBM1D_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolyBM1D::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolyBM1D::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolyBM1D::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyBM1D_ParallelPruning" )
      .derives<QuadraticPolyBM1D::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyBM1D_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolyBM1D::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolyBM1D::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolyBM1D::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolyBM1D::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolyBM1D::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolyBM1D::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolyBM1D::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolyBM1D>( "PCMBaseCpp__QuadraticPolyBM1D" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyBM1D)
      .method( "TraverseTree", &QuadraticPolyBM1D::TraverseTree )
      .method( "StateAtNode", &QuadraticPolyBM1D::StateAtNode )
      .property( "tree", &QuadraticPolyBM1D::tree )
      .property( "algorithm", &QuadraticPolyBM1D::algorithm )
    ;
  }


QuadraticPolyOU* CreateQuadraticPolyOU(
    arma::mat const& X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyOU::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
  }
  
  typename QuadraticPolyOU::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, pObjs.Pc, 
      pObjs.RModel, std::vector<std::string>(), 
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyOU(pObjs.br_0, pObjs.br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyOU) {
  Rcpp::class_<QuadraticPolyOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyOU_Tree" )
  .property("num_nodes", &QuadraticPolyOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyOU::TreeType>( "PCMBaseCpp__QuadraticPolyOU_OrderedTree" )
    .derives<QuadraticPolyOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyOU_ParallelPruning" )
    .derives<QuadraticPolyOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyOU>( "PCMBaseCpp__QuadraticPolyOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyOU)
    .method( "TraverseTree", &QuadraticPolyOU::TraverseTree )
    .method( "StateAtNode", &QuadraticPolyOU::StateAtNode )
    .property( "tree", &QuadraticPolyOU::tree )
    .property( "algorithm", &QuadraticPolyOU::algorithm )
  ;
}

QuadraticPolyOU1D* CreateQuadraticPolyOU1D(
    arma::mat const& X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyOU1D::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
  }
  
  typename QuadraticPolyOU1D::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, 
      pObjs.RModel, std::vector<std::string>(), 
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  return new QuadraticPolyOU1D(pObjs.br_0, pObjs.br_1, lengths, data);
}


  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyOU1D::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolyOU1D) {
    Rcpp::class_<QuadraticPolyOU1D::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyOU1D_Tree" )
    .property("num_nodes", &QuadraticPolyOU1D::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolyOU1D::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolyOU1D::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolyOU1D::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolyOU1D::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolyOU1D::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolyOU1D::TreeType>( "PCMBaseCpp__QuadraticPolyOU1D_OrderedTree" )
      .derives<QuadraticPolyOU1D::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyOU1D_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolyOU1D::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolyOU1D::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolyOU1D::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolyOU1D::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolyOU1D::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolyOU1D::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyOU1D_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolyOU1D::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolyOU1D::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolyOU1D::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyOU1D_ParallelPruning" )
      .derives<QuadraticPolyOU1D::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyOU1D_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolyOU1D::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolyOU1D::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolyOU1D::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolyOU1D::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolyOU1D::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolyOU1D::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolyOU1D::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolyOU1D>( "PCMBaseCpp__QuadraticPolyOU1D" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyOU1D)
      .method( "TraverseTree", &QuadraticPolyOU1D::TraverseTree )
      .method( "StateAtNode", &QuadraticPolyOU1D::StateAtNode )
      .property( "tree", &QuadraticPolyOU1D::tree )
      .property( "algorithm", &QuadraticPolyOU1D::algorithm )
    ;
  }




QuadraticPolyJOU* CreateQuadraticPolyJOU(
    arma::mat const& X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyJOU::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
    lengths[i].jump_ = pObjs.jumps[i];
  }
  
  typename QuadraticPolyJOU::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, pObjs.Pc, 
      pObjs.RModel, std::vector<std::string>(), 
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyJOU(pObjs.br_0, pObjs.br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyJOU::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyJOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyJOU) {
  Rcpp::class_<QuadraticPolyJOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyJOU_Tree" )
  .property("num_nodes", &QuadraticPolyJOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyJOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyJOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyJOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyJOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyJOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyJOU::TreeType>( "PCMBaseCpp__QuadraticPolyJOU_OrderedTree" )
    .derives<QuadraticPolyJOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyJOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyJOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyJOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyJOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyJOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyJOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyJOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyJOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyJOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyJOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyJOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyJOU_ParallelPruning" )
    .derives<QuadraticPolyJOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyJOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyJOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyJOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyJOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyJOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyJOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyJOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyJOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyJOU>( "PCMBaseCpp__QuadraticPolyJOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyJOU)
    .method( "TraverseTree", &QuadraticPolyJOU::TraverseTree )
    .method( "StateAtNode", &QuadraticPolyJOU::StateAtNode )
    .property( "tree", &QuadraticPolyJOU::tree )
    .property( "algorithm", &QuadraticPolyJOU::algorithm )
  ;
}

QuadraticPolyDOU* CreateQuadraticPolyDOU(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyDOU::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
  }
  
  typename QuadraticPolyDOU::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, pObjs.Pc, 
      pObjs.RModel, std::vector<std::string>(), 
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyDOU(pObjs.br_0, pObjs.br_1, lengths, data);;
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyDOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyDOU) {
  Rcpp::class_<QuadraticPolyDOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyDOU_Tree" )
  .property("num_nodes", &QuadraticPolyDOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyDOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyDOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyDOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyDOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyDOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyDOU::TreeType>( "PCMBaseCpp__QuadraticPolyDOU_OrderedTree" )
    .derives<QuadraticPolyDOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyDOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyDOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyDOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyDOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyDOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyDOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyDOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyDOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyDOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyDOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyDOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyDOU_ParallelPruning" )
    .derives<QuadraticPolyDOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyDOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyDOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyDOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyDOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyDOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyDOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyDOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyDOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyDOU>( "PCMBaseCpp__QuadraticPolyDOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyDOU)
    .method( "TraverseTree", &QuadraticPolyDOU::TraverseTree )
    .method( "StateAtNode", &QuadraticPolyDOU::StateAtNode )
    .property( "tree", &QuadraticPolyDOU::tree )
    .property( "algorithm", &QuadraticPolyDOU::algorithm )
  ;
}


QuadraticPolyMixedGaussian* CreateQuadraticPolyMixedGaussian(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo,
    std::vector<std::string> const& regimeModels) {
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyMixedGaussian::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
    lengths[i].jump_ = pObjs.jumps[i];
  }
  
  typename QuadraticPolyMixedGaussian::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, pObjs.Pc, 
      pObjs.RModel, 
      regimeModels,
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyMixedGaussian(pObjs.br_0, pObjs.br_1, lengths, data);
}

  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyMixedGaussian::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolyMixedGaussian) {
    Rcpp::class_<QuadraticPolyMixedGaussian::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyMixedGaussian_Tree" )
    .property("num_nodes", &QuadraticPolyMixedGaussian::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolyMixedGaussian::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolyMixedGaussian::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolyMixedGaussian::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolyMixedGaussian::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolyMixedGaussian::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::TreeType>( "PCMBaseCpp__QuadraticPolyMixedGaussian_OrderedTree" )
      .derives<QuadraticPolyMixedGaussian::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyMixedGaussian_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolyMixedGaussian::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolyMixedGaussian::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolyMixedGaussian::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolyMixedGaussian::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolyMixedGaussian::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyMixedGaussian_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolyMixedGaussian::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolyMixedGaussian::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyMixedGaussian_ParallelPruning" )
      .derives<QuadraticPolyMixedGaussian::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyMixedGaussian_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolyMixedGaussian::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolyMixedGaussian::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolyMixedGaussian::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolyMixedGaussian::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolyMixedGaussian::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolyMixedGaussian::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolyMixedGaussian::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian>( "PCMBaseCpp__QuadraticPolyMixedGaussian" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyMixedGaussian)
      .method( "TraverseTree", &QuadraticPolyMixedGaussian::TraverseTree )
      .method( "StateAtNode", &QuadraticPolyMixedGaussian::StateAtNode )
      .property( "tree", &QuadraticPolyMixedGaussian::tree )
      .property( "algorithm", &QuadraticPolyMixedGaussian::algorithm )
    ;
  }

QuadraticPolyMixedGaussian1D* CreateQuadraticPolyMixedGaussian1D(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo,
    std::vector<std::string> const& regimeModels) {
  
  ParsedRObjects pObjs(X, tree, model, metaInfo);
  
  vector<typename QuadraticPolyMixedGaussian1D::LengthType> lengths(pObjs.num_branches);
  
  for(arma::uword i = 0; i < pObjs.num_branches; ++i) {
    lengths[i].length_ = pObjs.t[i];
    lengths[i].regime_ = pObjs.regimes[i] - 1;
    lengths[i].jump_ = pObjs.jumps[i];
  }
  
  typename QuadraticPolyMixedGaussian1D::DataType data(
      pObjs.tip_names, pObjs.X, pObjs.VE, 
      pObjs.RModel, 
      regimeModels,
      pObjs.threshold_SV, pObjs.threshold_EV, 
      pObjs.threshold_skip_singular, pObjs.skip_singular,
      pObjs.transpose_Sigma_x,
      pObjs.threshold_Lambda_ij,
      pObjs.NA_double_);
  
  return new QuadraticPolyMixedGaussian1D(pObjs.br_0, pObjs.br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyMixedGaussian1D::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolyMixedGaussian1D) {
    Rcpp::class_<QuadraticPolyMixedGaussian1D::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyMixedGaussian1D_Tree" )
    .property("num_nodes", &QuadraticPolyMixedGaussian1D::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolyMixedGaussian1D::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolyMixedGaussian1D::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolyMixedGaussian1D::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolyMixedGaussian1D::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolyMixedGaussian1D::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian1D::TreeType>( "PCMBaseCpp__QuadraticPolyMixedGaussian1D_OrderedTree" )
      .derives<QuadraticPolyMixedGaussian1D::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyMixedGaussian1D_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolyMixedGaussian1D::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolyMixedGaussian1D::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolyMixedGaussian1D::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolyMixedGaussian1D::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolyMixedGaussian1D::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian1D::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyMixedGaussian1D_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolyMixedGaussian1D::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolyMixedGaussian1D::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian1D::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyMixedGaussian1D_ParallelPruning" )
      .derives<QuadraticPolyMixedGaussian1D::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyMixedGaussian1D_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolyMixedGaussian1D::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolyMixedGaussian1D::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolyMixedGaussian1D::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolyMixedGaussian1D::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolyMixedGaussian1D::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolyMixedGaussian1D::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolyMixedGaussian1D::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian1D>( "PCMBaseCpp__QuadraticPolyMixedGaussian1D" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyMixedGaussian1D)
      .method( "TraverseTree", &QuadraticPolyMixedGaussian1D::TraverseTree )
      .method( "StateAtNode", &QuadraticPolyMixedGaussian1D::StateAtNode )
      .property( "tree", &QuadraticPolyMixedGaussian1D::tree )
      .property( "algorithm", &QuadraticPolyMixedGaussian1D::algorithm )
    ;
  }
