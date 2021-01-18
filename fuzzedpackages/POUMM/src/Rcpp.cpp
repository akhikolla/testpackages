/**
 *  Rcpp.cpp
 *  POUMM
 *
 * Copyright 2015-2019 Venelin Mitov
 *
 * This file is part of the R-package POUMM: The Phylogenetic Ornstein-Uhlenbeck Mixed Model
 *
 * POUMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * POUMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with POUMM.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */

#include <Rcpp.h>

#include "AbcPOUMM.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]

// BEGIN: Needed for r-devel (R 3.4)
void R_init_POUMM(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}


void R_unload_POUMM(DllInfo *info) {
  /* Release resources. */
}
// END Needed for r-devel (R 3.4)

using namespace SPLITT;


ParallelPruningAbcPOUMM* CreateParallelPruningAbcPOUMM(
    Rcpp::List const& tree, vec const& z, vec const& se) {
  Rcpp::IntegerMatrix branches = tree["edge"];
  uvec parents(branches.column(0).begin(), branches.column(0).end());
  uvec daughters(branches.column(1).begin(), branches.column(1).end());
  vec t = Rcpp::as<vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  uvec tip_names = Seq(uint(1), num_tips);
  typename ParallelPruningAbcPOUMM::DataType data(tip_names, z, se);
  return new ParallelPruningAbcPOUMM(parents, daughters, t, data);
}

RCPP_EXPOSED_CLASS_NODECL(ParallelPruningAbcPOUMM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(ParallelPruningAbcPOUMM::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(ParallelPruningAbcPOUMM::AlgorithmType)

RCPP_MODULE(POUMM_AbcPOUMM) {
  Rcpp::class_<ParallelPruningAbcPOUMM::TreeType::Tree> ( "POUMM_Tree" )
  .property("num_nodes", &ParallelPruningAbcPOUMM::TreeType::Tree::num_nodes )
  .property("num_tips", &ParallelPruningAbcPOUMM::TreeType::Tree::num_tips )
  .method("LengthOfBranch", &ParallelPruningAbcPOUMM::TreeType::Tree::LengthOfBranch )
  .method("FindNodeWithId", &ParallelPruningAbcPOUMM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &ParallelPruningAbcPOUMM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &ParallelPruningAbcPOUMM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &ParallelPruningAbcPOUMM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::TreeType>( "POUMM_OrderedTree" )
    .derives<ParallelPruningAbcPOUMM::TreeType::Tree> ( "POUMM_Tree" )
    .method("RangeIdPruneNode", &ParallelPruningAbcPOUMM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &ParallelPruningAbcPOUMM::TreeType::RangeIdVisitNode )
    .property("num_levels", &ParallelPruningAbcPOUMM::TreeType::num_levels )
    .property("num_parallel_ranges_prune", &ParallelPruningAbcPOUMM::TreeType::num_parallel_ranges_prune )
    .property("ranges_id_visit", &ParallelPruningAbcPOUMM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &ParallelPruningAbcPOUMM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::AlgorithmType::ParentType>( "POUMM_TraversalAlgorithm" )
    .property( "VersionOPENMP", &ParallelPruningAbcPOUMM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &ParallelPruningAbcPOUMM::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::AlgorithmType> ( "POUMM_ParallelPruning" )
    .derives<ParallelPruningAbcPOUMM::AlgorithmType::ParentType>( "POUMM_TraversalAlgorithm" )
    .method( "ModeAutoStep", &ParallelPruningAbcPOUMM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &ParallelPruningAbcPOUMM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &ParallelPruningAbcPOUMM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &ParallelPruningAbcPOUMM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &ParallelPruningAbcPOUMM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &ParallelPruningAbcPOUMM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &ParallelPruningAbcPOUMM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::TraversalSpecificationType> ( "POUMM_PruningSpec" )
  .field( "a", &ParallelPruningAbcPOUMM::TraversalSpecificationType::a )
  .field( "b", &ParallelPruningAbcPOUMM::TraversalSpecificationType::b )
  .field( "c", &ParallelPruningAbcPOUMM::TraversalSpecificationType::c )
  .field( "se", &ParallelPruningAbcPOUMM::TraversalSpecificationType::se )
  .field( "z", &ParallelPruningAbcPOUMM::TraversalSpecificationType::z )
  .field( "alpha", &ParallelPruningAbcPOUMM::TraversalSpecificationType::alpha )
  .field( "sigma2", &ParallelPruningAbcPOUMM::TraversalSpecificationType::sigma2 )
  .field( "sigmae2", &ParallelPruningAbcPOUMM::TraversalSpecificationType::sigmae2 )
  .field( "theta", &ParallelPruningAbcPOUMM::TraversalSpecificationType::theta )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM>( "POUMM_AbcPOUMM" )
  .factory<Rcpp::List const&, vec const&, vec const&>(&CreateParallelPruningAbcPOUMM)
  .method( "DoPruning", &ParallelPruningAbcPOUMM::TraverseTree )
  .property( "tree", &ParallelPruningAbcPOUMM::tree )
  .property( "spec", &ParallelPruningAbcPOUMM::spec )
  .property( "algorithm", &ParallelPruningAbcPOUMM::algorithm )
  ;
}
