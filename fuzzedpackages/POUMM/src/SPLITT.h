/*
 *  SPLITT.h
 *  SPLITT
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of SPLITT: a generic C++ library for Serial and Parallel
 * Lineage Traversal of Trees.
 *
 * SPLITT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * SPLITT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SPLITT.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */

#ifndef SPLITT_SPLITT_H_
#define SPLITT_SPLITT_H_

#include <algorithm>
#include <vector>
#include <array>
#include <math.h>
#include <sstream>
#include <limits>
#include <numeric>
#include <chrono>
#include <unordered_map>
#include <mutex>
#include <condition_variable>

#ifdef _OPENMP

// Need to decide wheter to use '#pragma omp for' or '#pragma omp for simd'
#if _OPENMP >= 201307  // OMP 4.0 or higher

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#else // #if _OPENMP >= 201307

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD  /*_Pragma("omp simd")*/

#endif // _OPENMP >= 201307

#include <omp.h>

#else // #ifdef _OPENMP

// the preprocessor directives should simply be ignored at compile-time
#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // #ifdef _OPENMP


//' @name SPLITT
//' @title SPLITT: A generic C++ library for Serial and Parallel Lineage Traversal of Trees
//' 
//' @description All basic types, functions and classes defined in SPLITT.
//' 
//' %\if{html}{\figure{UmlDiagram8.pdf}{options: width=100\%}}
//' %\if{latex}{\figure{UmlDiagram8.pdf}{options: width=15cm}}
//' 
//' @section Basic types:
//' \describe{
//' \item{\link[=SPLITT::uint]{uint}}{}
//' \item{\link[=SPLITT::uvec]{uvec}}{}
//' \item{\link[=SPLITT::vec]{vec}}{}
//' \item{\link[=SPLITT::bvec]{bvec}}{}
//' }
//'
//' @section Global constants:
//' \describe{
//' \item{\link[=SPLITT::G_NA_UINT]{G_NA_UINT}}{}
//' \item{\link[=SPLITT::G_EMPTY_UVEC]{G_EMPTY_UVEC}}{}
//' \item{\link[=SPLITT::G_PI]{G_PI}}{}
//' }
//' @section Global functions:
//' \describe{
//' \item{\link[=SPLITT::SortIndices]{SortIndices}}{}
//' \item{\link[=SPLITT::At]{At}}{}
//' \item{\link[=SPLITT::Match]{Match}}{}
//' \item{\link[=SPLITT::Seq]{Seq}}{}
//' \item{\link[=SPLITT::IsNA]{IsNA}}{}
//' \item{\link[=SPLITT::NotIsNA]{NotIsNA}}{}
//' }
//' 
//' @section Classes:
//' \describe{
//' \item{\link[=SPLITT::TraversalSpecification]{TraversalSpecification}}{}
//' \item{\link[=SPLITT::TraversalTask]{TraversalTask}}{}
//' \item{\link[=SPLITT::TraversalTaskLightweight]{TraversalTaskLightweight}}{}
//' \item{\link[=SPLITT::Tree]{Tree}}{}
//' \item{\link[=SPLITT::OrderedTree]{OrderedTree}}{}
//' \item{\link[=SPLITT::ThreadExceptionHandler]{ThreadExceptionHandler}}{}
//' }
//' 
//' // to enable generation of man-pages use Rcpp::export below (no spaces between :'s).
//' [[Rcpp : : export]]
namespace SPLITT{


/*******************************************************************************
 
 Basic types

******************************************************************************/


//' @name SPLITT::uint
//' @backref src/SPLITT.h
//' @title An alias for the \code{'unsigned int'} basic type.
//' @description 
//' \code{typedef unsigned int uint;}
//' @family basic types
typedef unsigned int uint;

//' @name SPLITT::uvec
//' @backref src/SPLITT.h
//' @title a vector of \code{\link[=SPLITT::uint]{uint}}'s.
//' @description 
//' \code{typedef }\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<uint> uvec;}
//' @family basic types
typedef std::vector<uint> uvec;

//' @name SPLITT::vec
//' @backref src/SPLITT.h
//' @title a vector of \code{double}s.
//' 
//' @description 
//' \code{typedef }\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<double> vec;}
//' 
//' @family basic types
//' @seealso \link{SPLITT}
typedef std::vector<double> vec;

//' @name SPLITT::bvec
//' @backref src/SPLITT.h
//' @title a vector of \code{bool}s.
//' 
//' @description 
//' \code{typedef }\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<bool> bvec;}
//' 
//' @family basic types
//' @seealso \link{SPLITT}
typedef std::vector<bool> bvec;

/*******************************************************************************
 
 Global constants

******************************************************************************/
#ifndef SPLITT_CONSTANTS_H
#define SPLITT_CONSTANTS_H

//' @name SPLITT::G_PI
//' @backref src/SPLITT.h
//' 
//' @title The mathematical constant Pi. 
//' 
//' @description Currently defined as
//' \code{const double G_PI = 3.1415926535897932385;}
//' 
//' @family global constants
//' 
//' @seealso \code{\link{SPLITT}}
const double G_PI = 3.1415926535897932385;

//' @name SPLITT::G_EMPTY_UVEC
//' @backref src/SPLITT.h
//' @title A global constant for the empty \code{\link[=SPLITT::uvec]{uvec}};
//' @description  
//' Currently defined as
//' \code{const uvec G_EMPTY_UVEC;}
//' 
//' @family global constants
//' @seealso \link{SPLITT}
const uvec G_EMPTY_UVEC;

//' @name SPLITT::G_NA_UINT
//' @backref src/SPLITT.h
//' @title A global constant for the integer NA.
//' 
//' @description 
//' Currently defined as
//' \code{const uint G_NA_UINT = std::numeric_limits<uint>::max();}
//' @family global constants
//' @seealso \link{SPLITT}
const uint G_NA_UINT = std::numeric_limits<uint>::max();

#endif // SPLITT_CONSTANTS_H

/*******************************************************************************
 
 Global functions

******************************************************************************/

//' @name SPLITT::SortIndices
//' @backref src/SPLITT.h
//' @title Indices in a vector ordered in ascending order of the corresponding values.
//' 
//' @description
//' \code{template <class VectorClass>
//'   inline std::vector<uint> SortIndices(VectorClass const& v);} 
//'   
//'   This is a template function. The template argument \code{VectorClass} must 
//'   be an index-able class, such as 
//'   \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector<T>}},
//'   where the type of the elements, \code{T}, must have operator <.
//' @param v a \code{const&} to a \code{VectorClass} object.
//' @return an \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::uint]{uint}>}. 
//' 
//' @family global functions
//' @seealso \link{SPLITT}
template <class VectorClass> inline std::vector<uint> SortIndices(VectorClass const& v) {
  // initialize original index locations
  std::vector<uint> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indices based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

//' @name SPLITT::At
//' @backref src/SPLITT.h
//' @title A sub-vector with the elements at given positions.
//' 
//' @description 
//' This is a template function with two overwrites:
//' 
//' \code{
//' template<class VectorValues, class VectorPositions>
//' inline VectorValues At(VectorValues const& v, VectorPositions const& positions);}
//' 
//' \code{
//' template<class VectorValues>
//' inline VectorValues At(VectorValues const& v, \link[=SPLITT::bvec]{bvec} const& mask);}
//'   
//' The template argument \code{VectorValues} must 
//'   be an index-able class, such as std::vector; the template argument 
//'   VectorPositions must have a forward iterator.
//' @param v a \code{const&} to a \code{VectorValues} object. 
//' @param positions a \code{const&} to \code{VectorPositions} object. The elements in positions
//'   must be convertible to indices in v.
//' @param mask a \code{\link[=SPLITT::bvec]{bvec} const&} of the same length as
//' \code{v} with elements equal to \code{true} at the positions to be returned.
//' 
//' @return a \code{VectorValues} object.
//' @family global functions
//' @seealso \link{SPLITT}
template<class VectorValues, class VectorPositions>
inline VectorValues At(VectorValues const& v, VectorPositions const& positions) {
  VectorValues sub;
  sub.resize(positions.size());

  size_t sub_i = 0;
  for(auto pit = positions.begin(); pit != positions.end(); pit++,sub_i++){
    sub[sub_i] = v[*pit];
  }
  return sub;
}

template<class VectorValues>
inline VectorValues At(VectorValues const& v, bvec const& mask) {
  if(mask.size() != v.size()) {
    throw std::length_error("ERR:01001:SPLITT:SPLITT.h:At:: bool vector mask should have the same length as v.");
  }

  size_t res_size = 0;
  for(auto b : mask) if(b) ++res_size;

  VectorValues sub(res_size);

  size_t sub_i = 0;
  for(uint i = 0; i < v.size(); i++){
    if(mask[i]) sub[sub_i++] = v[i];
  }
  return sub;
}

//' @name SPLITT::Match
//' @backref src/SPLITT.h
//' @title Match the first occurrences of elements in a vector of integers.
//' 
//' @description
//' This is a template function optimized for searching integer types (internal 
//' use only).
//' 
//' \code{
//' template<class VectorValues, class PosType>
//' inline std::vector<PosType> Match(VectorValues const& x, 
//' VectorValues const& table, PosType const& NA);}
//'
//' For each element of \code{x} return the index of its first occurence in 
//' \code{table} or \code{NA} if the element is not found in \code{table} or is 
//' equal to \code{NA}. It is assumed that \code{x} does not have duplicated 
//' elements or \code{NA} elements.
//' 
//' The template argument \code{VectorValues} must 
//'   be an index-able class, such as 
//'   \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<T>},
//'   where \code{T} must be an integer type.
//'   The template argument \code{PosType} must be an integer type that can be used
//'   as index type in \code{VectorValues}.
//' @param x \code{VectorValues const&} of values to search for in \code{table}. 
//' This vector should not have duplicates or \code{NA}s.
//' @param table \code{VectorValues const&} of (possibly duplicated) values to 
//' matched against the elements in \code{x}. 
//' @param NA \code{PosType const&} an integer type that can be used as index-type
//' in \code{x} and \code{table}.
//' @return an \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<PosType>}.
//' @family global functions
//' @seealso \link{SPLITT}
template<class VectorValues, class PosType>
inline std::vector<PosType> Match(
    VectorValues const& x, VectorValues const& table, PosType const& NA) {
  auto minmax_x = std::minmax_element(x.begin(), x.end());
  std::vector<PosType> index(*minmax_x.second - *minmax_x.first + 1, NA);
  for(PosType i = 0; i < table.size(); ++i) {
    if(table[i] >= *minmax_x.first && table[i] <= *minmax_x.second &&
       index[table[i] - *minmax_x.first] == NA) {
      index[table[i] - *minmax_x.first] = i;
    }
  }

  std::vector<PosType> positions(x.size());
  for(size_t i = 0; i < x.size(); ++i) {
    positions[i] = index[x[i] - *minmax_x.first];
  }
  return positions;
}

//' @name SPLITT::Seq
//' @backref src/SPLITT.h
//' @title Generate a sequence of consecutive integers.
//' 
//' @description
//' This is a template function.
//' 
//' \code{
//' template<class T> inline std::vector<T> Seq(T const& first, T const& last);}
//' 
//' Currently the function is implemented as a wrapper for 
//' \href{https://en.cppreference.com/w/cpp/algorithm/iota}{\code{std::iota}}\code{<T>}.
//' @param first,last \code{T const&} first and last elements in the sequence.
//' @return 
//' a \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<T>} object.
//' @family global functions
//' @seealso \link{SPLITT}
template<class T>
inline std::vector<T> Seq(T const& first, T const& last) {
  std::vector<T> res(last-first+1);
  std::iota(res.begin(), res.end(), first);
  return res;
}

//' @name SPLITT::IsNA
//' @backref src/SPLITT.h
//' @title Check a vector for \code{NA}s.
//' 
//' @description
//' This is a template function with a specification for the type \link[=SPLITT::uint]{uint}.
//' 
//' \code{
//' template<class T> inline bvec IsNA(std::vector<T> const& x, T const& NA);
//' 
//' inline bvec IsNA(uvec const& x);}
//' 
//' The element type \code{T} must define an \code{==} operator. 
//' 
//' @param x a \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<T> const&} 
//' of elements to be checked. 
//' @param NA \code{T const&} value specifying NA (can be any value).
//' 
//' @return a \code{\link[=SPLITT::bvec]{bvec}} of the same length as \code{x}
//' with \code{false} everywhere, except for the positions where there are \code{NA}s.
//' @family global functions
//' @seealso \link{SPLITT}
template<class T>
inline bvec IsNA(std::vector<T> const& x, T const& NA) {
  bvec res(x.size(), true);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA) res[i] = true;
  }
  return res;
}

inline bvec IsNA(uvec const& x) {
  return IsNA(x, G_NA_UINT);
}

//' @name SPLITT::NotIsNA
//' @backref src/SPLITT.h
//' @title Check a vector for \code{NA}s.
//' 
//' @description
//' This is a template function with a specification for the type \link[=SPLITT::uint]{uint}.
//' 
//' \code{
//' template<class T> inline bvec NotIsNA(std::vector<T> const& x, T const& NA);
//' 
//' inline bvec NotIsNA(uvec const& x);}
//' 
//' The element type \code{T} must define an \code{==} operator. 
//' 
//' @param x a \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<T> const&} 
//' of elements to be checked. 
//' @param NA \code{T const&} value specifying NA (can be any value).
//' 
//' @return a \code{\link[=SPLITT::bvec]{bvec}} of the same length as \code{x}
//' with \code{true} everywhere, except for the positions where there are \code{NA}s.
//' @family global functions
//' @seealso \link{SPLITT}
template<class T>
inline bvec NotIsNA(std::vector<T> const& x, T const& NA) {
  bvec res(x.size(), true);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA) res[i] = false;
  }
  return res;
}

inline bvec NotIsNA(uvec const& x) {
  return NotIsNA(x, G_NA_UINT);
}





/*******************************************************************************
 
                                     Classes
 
 ******************************************************************************/

//' 
//' @name SPLITT::TraversalSpecification
//' @backref src/SPLITT.h
//' 
//' @title An abstract base class for specifying node traversal operations.
//' @description This is an abstract class (not intended for instantiation).
//' The user must provide a TraversalSpecificationImplementation class implementing
//' this class' methods as described below. It is
//' recommended to inherit from this class, but it is not obligatory, because 
//' this is is not checked during compilation. Use the documentation of this 
//' class as a guide through the steps of writing a tree traversal specification
//' based on \code{\link{SPLITT}}.
//' @seealso \code{\link{SPLITT}}
//' 
template<class Tree> class TraversalSpecification {
protected:
  // A reference to a Tree available for inheriting classes
  Tree const& ref_tree_;
  // A protected constructor that initializes the tree-reference. This constructor
  // must be called explicitly in the initalization list of inheriting class constructors.
  TraversalSpecification(Tree const& tree): ref_tree_(tree) {}
public:
  // public typedefs. These typedefs must be provided by an implementation class.
  // 1. typedef Tree TreeType;
  // 2. typedef PostOrderTraversal<ImlementationClass> AlgorithmType;
  // 3. typedef ImplementationSpecificParameterType ParameterType;
  // 4. typedef ImplementationSpecificDataType DataType;
  // 5. typedef ImplementationSpecificStateType StateType;
  
  
  // The following methods must be present any implementation
  // 6. constructor: will be called by a TraversalTask object; Here, it is
  // commented out, because the DataType is not known.
  // ImplementationClassName(TreeType & tree, DataType & input_data) :
  //   TraversalSpecification(tree) {
  //     implementation specific initialization using the tree and the input_data.
  // }
  
  
  // The following methods get called by the TraversalAlgorithm implementation:
  
  // 7. Setting the model parameters prior to starting the pruning procedure on the tree.
  // This method is called by the TraversalTask.TraverseTree(ParamterType const&, uint mode)
  // method. The method declaration is commented out because ParameterType is not known
  // and must be specified by the implementing class.
  // void SetParameter(ParameterType const& par);
  
  // 8. InitNode(i) is called on each node in the tree right after SetParameter and
  // before any of the VisitNode and PruneNode has been called. There is no predefined
  // order of the calls to InitNode and they may be executed in parallel. Therefore, only
  // node-specific data initialization, including the length of the branch
  // leading to node i, can take place in this method.
  void InitNode(uint i) {}
  
  
  // 9. VisitNode(i) is called on each tip or internal node (EXCLUDING THE ROOT),
  // in the tree after PruneNode has been called on each descendant of i.
  // The method is the perfect place to calculate the state of node i using the
  // pre-calculated states of its descendants. Although, it is guaranteed
  // that VisitNode(i) is called before VisitNode(i_parent), this method SHOULD NOT BE USED
  // FOR ALTERING THE STATE of i_parent, because this would conflict with
  // a concurrent execution of VisitNode on a sibling of i (see also PruneNode).
  void VisitNode(uint i) {}
  
  // 10. PruneNode(i, i_parent) is called on each tip or internal node (EXCLUDING THE ROOT)
  // after VisitNode(i) and in sync with PruneNode(k, i_parent), for any sibling k of i.
  // Thus, it is safe to use PruneNode to update the state of i_parent.
  void PruneNode(uint i, uint i_parent) {}
  
  // 11. StateType StateAtRoot() is called after PruneNode has been called on each
  // direct descendant of the root node. If necessary, VisitNode(i_root) can be called
  // here, in order to calculate the final state of the root. The value returned by this
  // function is also returned by the TraversalTask.TraverseTree(ParameterType const& par, uint mode)
  // method.
};

// 12. After the class TraversalSpecificationImplementation has been defined it is
// time to specify the TraversalTask template. This is not obligatory but can be very
// convinient for creating TraversalTask objects with the user specific implementation
// and to call their TraverseTree method.
// typedef TraversalTask<TraversalSpecificationImplementation> > MyTraversalTask;

//' @name SPLITT::TraversalTask
//' @backref src/SPLITT.h
//' 
//' @title A composite of a tree, a traversal specification and a traversal algorithm.
//' 
//' @description
//' \code{
//' template<class TraversalSpecification> class TraversalTask;
//' }
//' @section Template Arguments:
//' \itemize{
//' \item{class TraversalSpecification}{see \code{\link[=SPLITT::TraversalSpecification]{TraversalSpecification}}.}
//' }
//' @section Public Methods:
//' \describe{
//' \item{\link[=SPLITT::TraversalTask::TraversalTaskLightweight]{TraversalTask}}{}
//' \item{\link[=SPLITT::TraversalTask::TraverseTree]{TraverseTree}}{}
//' \item{\link[=SPLITT::TraversalTask::tree]{tree}}{}
//' \item{\link[=SPLITT::TraversalTask::spec]{spec}}{}
//' \item{\link[=SPLITT::TraversalTask::algorithm]{algorithm}}{}
//' }
//' @seealso \link{SPLITT::TraversalSpecification}  
//' @seealso \link{SPLITT} 
template<class TraversalSpecification>
class TraversalTask {
public:
  typedef TraversalSpecification TraversalSpecificationType;
  typedef typename TraversalSpecification::TreeType TreeType;
  typedef typename TraversalSpecification::AlgorithmType AlgorithmType;
  typedef typename AlgorithmType::ModeType ModeType;
  typedef typename TreeType::NodeType NodeType;
  typedef typename TreeType::LengthType LengthType;
  typedef typename TraversalSpecificationType::DataType DataType;
  typedef typename TraversalSpecificationType::ParameterType ParameterType;
  typedef typename TraversalSpecificationType::StateType StateType;
  
  TraversalTask(
    std::vector<NodeType> const& branch_start_nodes,
    std::vector<NodeType> const& branch_end_nodes,
    std::vector<LengthType> const& branch_lengths,
    DataType const& data):
    tree_(branch_start_nodes, branch_end_nodes, branch_lengths),
    spec_(tree_, data),
    algorithm_(tree_, spec_) {}
  
  StateType TraverseTree(ParameterType const& par, uint mode) {
    spec_.SetParameter(par);
    algorithm_.TraverseTree(static_cast<ModeType>(mode));
    return spec_.StateAtRoot();
  }
  
  StateType StateAtNode(uint i) {
    return spec_.StateAtNode(i);
  }
  
  TreeType & tree() {
    return tree_;
  }
  TraversalSpecification & spec() {
    return spec_;
  }
  AlgorithmType & algorithm() {
    return algorithm_;
  }
  
protected:
  TreeType tree_;
  TraversalSpecification spec_;
  AlgorithmType algorithm_;
};

//' @name SPLITT::TraversalTaskLightweight
//' @backref src/SPLITT.h
//' 
//' @title A lighter TraversalTask class which gets a reference to an already 
//' constructed tree.
//' 
//' @description
//' \code{
//' template<class TraversalSpecification> class TraversalTaskLightweight;
//' }
//' @section Template Arguments:
//' \itemize{
//' \item{class TraversalSpecification}{see \code{\link[=SPLITT::TraversalSpecification]{TraversalSpecification}}.}
//' }
//' @section Public Methods:
//' \describe{
//' \item{\link[=SPLITT::TraversalTaskLightweight::TraversalTaskLightweight]{TraversalTaskLightweight}}{}
//' \item{\link[=SPLITT::TraversalTaskLightweight::TraverseTree]{TraverseTree}}{}
//' \item{\link[=SPLITT::TraversalTaskLightweight::spec]{spec}}{}
//' \item{\link[=SPLITT::TraversalTaskLightweight::algorithm]{algorithm}}{}
//' }
//' @seealso \link{SPLITT::TraversalSpecification}  
//' @seealso \link{SPLITT} 
template<class TraversalSpecification>
class TraversalTaskLightweight {
public:
  typedef TraversalSpecification TraversalSpecificationType;
  typedef typename TraversalSpecification::TreeType TreeType;
  typedef typename TraversalSpecification::AlgorithmType AlgorithmType;
  typedef typename AlgorithmType::ModeType ModeType;
  typedef typename TreeType::NodeType NodeType;
  typedef typename TreeType::LengthType LengthType;
  typedef typename TraversalSpecificationType::DataType DataType;
  typedef typename TraversalSpecificationType::ParameterType ParameterType;
  typedef typename TraversalSpecificationType::StateType StateType;
  
  TraversalTaskLightweight(
    TreeType const& tree,
    DataType const& data):
    tree_(tree),
    spec_(tree_, data),
    algorithm_(tree_, spec_) {}
  
  StateType TraverseTree(ParameterType const& par, uint mode) {
    spec_.SetParameter(par);
    algorithm_.TraverseTree(static_cast<ModeType>(mode));
    return spec_.StateAtRoot();
  }
  
  TraversalSpecification & spec() {
    return spec_;
  }
  AlgorithmType & algorithm() {
    return algorithm_;
  }
  
protected:
  TreeType const& tree_;
  TraversalSpecification spec_;
  AlgorithmType algorithm_;
};


//' @name SPLITT::Tree
//' 
//' @title Base class \code{Tree}.
//' 
//' @description A generic C++ template class defining the data structure and 
//' basic operations with a tree. 
//'   
//' \code{template<class Node, class Length>class Tree;}
//' 
//' @section Template Arguments:
//' \itemize{
//' \item{class Node}{see \code{\link[=SPLITT::Tree::NodeType]{NodeType}}.}
//' \item{class Length}{see \code{\link[=SPLITT::Tree::LengthType]{LengthType}}.}
//' }
//' @section Public Methods:
//' \describe{
//' \item{\link[=SPLITT::Tree::Tree]{Tree}}{}
//' \item{\link[=SPLITT::Tree::num_tips]{num_tips}}{}
//' \item{\link[=SPLITT::Tree::num_nodes]{num_nodes}}{}
//' \item{\link[=SPLITT::Tree::BranchLengths]{BranchLengths}}{}
//' \item{\link[=SPLITT::Tree::FindChildren]{FindChildren}}{}
//' \item{\link[=SPLITT::Tree::FindIdOfNode]{FindIdOfNode}}{}
//' \item{\link[=SPLITT::Tree::FindIdOfParent]{FindIdOfParent}}{}
//' \item{\link[=SPLITT::Tree::FindNodeWithId]{FindNodeWithId}}{}
//' \item{\link[=SPLITT::Tree::HasBranchLengths]{HasBranchLengths}}{}
//' \item{\link[=SPLITT::Tree::LengthOfBranch]{LengthOfBranch}}{}
//' \item{\link[=SPLITT::Tree::SetBranchLengths]{SetBranchLengths}}{}
//' \item{\link[=SPLITT::Tree::SetLengthOfBranch]{SetLengthOfBranch}}{}
//' \item{\link[=SPLITT::Tree::OrderNodesPosType]{OrderNodesPosType}}{}
//' \item{\link[=SPLITT::Tree::OrderNodes]{OrderNodes}}{}
//' }
//' 
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}}
//' @seealso \link{SPLITT}
template<class Node, class Length>
class Tree {
public:
//' @name SPLITT::Tree::NodeType 
//' @title Abstract type for nodes in the tree.
//' @description A public typedef in class \code{\link[=SPLITT::Tree]{Tree}}. A synonym for the 
//'   template argument \code{Node}. Defines a 
//'   hash-able type such as \code{int} or 
//'   \href{http://en.cppreference.com/w/cpp/string/basic_string}{\code{std::string}}. 
//'   Specifically, this should be a type for which 
//'   \href{http://en.cppreference.com/w/cpp/utility/hash}{\code{std::hash}}
//'   specialization does exist. This is the application-specific node-type. 
//'   The branches in the tree are defined as couples of nodes 
//'   <branch_start_node, branch_end_node>.
//' @details
//'   During the construction of \code{\link[=SPLITT::Tree]{Tree}} object, the 
//'   nodes are assigned \code{unsigned int} ids from 0 to M-1 (M being the 
//'   number of nodes in the tree). 
//'   
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
typedef Node NodeType;
  
//' @name SPLITT::Tree::LengthType
//'
//' @title Abstract type for the branch-lengths in a \code{\link[=SPLITT::Tree]{Tree}}.
//' 
//' @description
//' \code{typedef Length LengthType;}
//' 
//' A public typedef in class \code{\link[=SPLITT::Tree]{Tree}}. 
//'   A synonym for the template argument Length. Defines a type that can be 
//'   associated with a branch. Can be a basic type, e.g. \code{double}, but also
//'   a composite of several attributes on a branch, such as a \code{double} 
//'   length and an \code{int} color.
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
typedef Length LengthType;

private:
  // private default constructor
  Tree() {}

protected:
  uint num_tips_;
  uint num_nodes_;
  uvec id_parent_;

  typedef std::unordered_map<NodeType, uint> MapType;
  MapType map_node_to_id_;
  std::vector<NodeType> map_id_to_node_;
  std::vector<LengthType> lengths_;
  std::vector<uvec> id_child_nodes_;

  void init_id_child_nodes() {
    id_child_nodes_ = std::vector<uvec>(this->num_nodes() - this->num_tips());

    // fill child vectors
    for(uint i = 0; i < this->num_nodes() - 1; i++) {
      id_child_nodes_[this->FindIdOfParent(i) - this->num_tips()].push_back(i);
    }
  }

public:
//' @name SPLITT::Tree::Tree
//' 
//' @title Constructor for class \code{\link[=SPLITT::Tree]{Tree}}.
//' 
//' @description 
//' \code{
//' Tree(std::vector<NodeType> const& branch_start_nodes,
//' std::vector<NodeType> const& branch_end_nodes,
//' std::vector<LengthType> const& branch_lengths);}
//' 
//' Constructs the tree object given a list of branches. The list of branches
//' is specified from the corresponding elements in the three vectors passed as
//' arguments. 
//' 
//' @param branch_start_nodes 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const&}: 
//'   starting node for every branch in the tree.
//' @param branch_end_nodes 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const&}:
//'   ending node for every branch in the tree; must be the same length as 
//'   \code{branch_start_nodes}.
//' @param branch_lengths 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const&}: 
//' lengths associated with the branches. Pass an empty vector for \code{branch_lengths} 
//' for a tree without branch lengths (i.e. only a topology).
//'
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
Tree(std::vector<NodeType> const& branch_start_nodes,
       std::vector<NodeType> const& branch_end_nodes,
       std::vector<LengthType> const& branch_lengths) {

    if(branch_start_nodes.size() != branch_end_nodes.size()) {
      std::ostringstream oss;
      oss<<"ERR:01011:SPLITT:SPLITT.h:Tree::"<<
        " branch_start_nodes and branch_end_nodes should be the same size, but were "
         <<branch_start_nodes.size()<<" and "<<
      branch_end_nodes.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    // There should be exactly num_nodes_ = number-of-branches + 1
    // distinct nodes
    // This is because each branch can be mapped to its ending node. The +1
    // corresponds to the root node, to which no branch points.
    this->num_nodes_ = branch_start_nodes.size() + 1;


    // we distinguish three types of nodes:
    enum NodeRole { ROOT, INTERNAL, TIP };

    // initially, we traverse the list of branches and order the nodes as they
    // appear during the traversal from 0 to num_nodes_-1. This order is done by
    // incrementing node_id_temp.
    uint node_id_temp = 0;

    std::vector<NodeRole> node_types(num_nodes_, ROOT);

    this->map_id_to_node_.resize(num_nodes_);

    this->map_node_to_id_.reserve(num_nodes_);

    uvec branch_starts_temp(branch_start_nodes.size(), G_NA_UINT);
    uvec branch_ends_temp(branch_start_nodes.size(), G_NA_UINT);
    uvec ending_at(num_nodes_, G_NA_UINT);

    std::vector<typename MapType::iterator> it_map_node_to_id_;
    it_map_node_to_id_.reserve(num_nodes_);


    for(uint i = 0; i < branch_start_nodes.size(); ++i) {
      if(branch_start_nodes[i] == branch_end_nodes[i]) {
        std::ostringstream oss;
        oss<<"ERR:01012:SPLITT:SPLITT.h:Tree:: Found a branch with the same start and end node ("<<
          branch_start_nodes[i]<<"). Not allowed. ";
        throw std::logic_error(oss.str());
      }


      auto it1 = map_node_to_id_.insert(
        std::pair<NodeType, uint>(branch_start_nodes[i], node_id_temp));

      if(it1.second) {
        // node encountered for the first time and inserted in the map_node_to_id_
        map_id_to_node_[node_id_temp] = branch_start_nodes[i];
        if(node_types[node_id_temp] == TIP) {
          node_types[node_id_temp] = INTERNAL;
        }
        branch_starts_temp[i] = node_id_temp;
        it_map_node_to_id_.push_back(it1.first);
        node_id_temp++;
      } else {
        // node encountered in a previous branch
        if(node_types[it1.first->second] == TIP) {
          // the previous encounter of the node was as a branch-end
          node_types[it1.first->second] = INTERNAL;
        } else {
          // do nothing
        }
        branch_starts_temp[i] = it1.first->second;
      }


      auto it2 = map_node_to_id_.insert(std::pair<NodeType, uint>(branch_end_nodes[i], node_id_temp));

      if(it2.second) {
        // node encountered for the first time and inserted in the map_node_to_id
        map_id_to_node_[node_id_temp] = branch_end_nodes[i];

        if(node_types[node_id_temp] == ROOT) {
          // not known if the node has descendants, so we set its type to TIP.
          node_types[node_id_temp] = TIP;
        }
        branch_ends_temp[i] = node_id_temp;
        ending_at[node_id_temp] = i;
        it_map_node_to_id_.push_back(it2.first);
        node_id_temp++;

      } else {
        // node has been previously encountered
        if(ending_at[it2.first->second] != G_NA_UINT) {
          std::ostringstream oss;
          oss<<"ERR:01013:SPLITT:SPLITT.h:Tree:: Found at least two branches ending at the same node ("<<
            it2.first->first<<"). Check for cycles or repeated branches. ";
          throw std::logic_error(oss.str());
        } else {
          if(node_types[it2.first->second] == ROOT) {
            // the previous enounters of the node were as branch-start -> set
            // the node's type to INTERNAL, because we know for sure that it
            // has descendants.
            node_types[it2.first->second] = INTERNAL;
          }
          branch_ends_temp[i] = it2.first->second;
          ending_at[it2.first->second] = i;
        }
      }
    }

    if(map_node_to_id_.size() != num_nodes_) {
      std::ostringstream oss;
      oss<<"ERR:01014:SPLITT:SPLITT.h:Tree:: The number of distinct nodes ("<<map_node_to_id_.size()<<
        ") should equal the number-of-branches+1 ("<<num_nodes_<<").";
      throw std::logic_error(oss.str());
    }

    auto num_roots = count(node_types.begin(), node_types.end(), ROOT);
    if(num_roots != 1) {
      std::ostringstream oss;
      oss<<"ERR:01015:SPLITT:SPLITT.h:Tree:: There should be exactly one ROOT node, but "<<num_roots<<
        " were found. Check for cycles or for multiple trees.";
      throw std::logic_error(oss.str());
    }

    this->num_tips_ = count(node_types.begin(), node_types.end(), TIP);
    if(num_tips_ == 0) {
      std::ostringstream oss;
      oss<<"ERR:01016:SPLITT:SPLITT.h:Tree:: There should be at least one TIP node, but none"<<
        " was found. Check for cycles.";
      throw std::logic_error(oss.str());
    }

    // assign new ids according to the following convention:
    // tips are numbered from 0 to num_tips_ - 1;
    // internal nodes are numbered from num_tips_ to num_nodes_ - 2;
    // root is numbered num_nodes_ - 1;
    std::vector<uint> node_ids(num_nodes_, G_NA_UINT);
    uint tip_no = 0, internal_no = num_tips_;
    for(uint i = 0; i < num_nodes_; i++) {
      if(node_types[i] == TIP) {
        node_ids[i] = tip_no;
        tip_no ++;
      } else if(node_types[i] == INTERNAL) {
        node_ids[i] = internal_no;
        internal_no++;
      } else {
        // Here node_types[i] == ROOT should be true
        node_ids[i] = num_nodes_ - 1;
      }
      //map_node_to_id_[map_id_to_node_[i]] = node_ids[i];
      it_map_node_to_id_[i]->second = node_ids[i];
    }


    this->map_id_to_node_ = At(map_id_to_node_, SortIndices(node_ids));

    this->id_parent_ = uvec(num_nodes_ - 1);

    if(branch_lengths.size() == num_nodes_ - 1) {
      this->lengths_ = std::vector<LengthType>(num_nodes_ - 1);
    } else if(branch_lengths.size() != 0) {
      std::ostringstream oss;
      oss<<"ERR:01017:SPLITT:SPLITT.h:Tree:: branch_lengths should be either empty or of size num_nodes_-1 ("<<
        num_nodes_-1<<") but is "<<branch_lengths.size()<<"."<<std::endl;
      throw std::invalid_argument(oss.str());
    }


    if(HasBranchLengths()) {
      for(uint i = 0; i < num_nodes_ - 1; i++) {
        uint branch_start_i = node_ids[branch_starts_temp[i]];
        uint branch_end_i = node_ids[branch_ends_temp[i]];
        id_parent_[branch_end_i] = branch_start_i;
        lengths_[branch_end_i] = branch_lengths[i];
      }
    } else {
      for(uint i = 0; i < num_nodes_ - 1; i++) {
        uint branch_start_i = node_ids[branch_starts_temp[i]];
        uint branch_end_i = node_ids[branch_ends_temp[i]];
        id_parent_[branch_end_i] = branch_start_i;
      }
    }

    init_id_child_nodes();
  }

//' @name SPLITT::Tree::num_nodes
//' @title Number of nodes in a tree.
//'  
//' @description 
//' \code{\link[=SPLITT::uint]{uint} Tree::num_nodes() const;}
//'  
//' @return a \code{\link[=SPLITT::uint]{uint}}: the numbef of nodes in the tree, 
//' including tips, internal nodes and the root.
//' 
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
uint num_nodes() const {
    return num_nodes_;
  }

//' @name SPLITT::Tree::num_tips
//' 
//' @title Number of tips (a.k.a. leaves) in the tree
//' @description
//' \code{\link[=SPLITT::uint]{uint} \link[=SPLITT::Tree]{Tree}::num_nodes() const;}
//' 
//' @return \code{\link[=SPLITT::uint]{uint}}, the number of tips in the tree.
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
uint num_tips() const {
    return num_tips_;
  }
  
//' @name SPLITT::Tree::HasBranchLengths
//' 
//' @title Does a tree has lengths associated with its branches?
//' @description
//' \code{bool HasBranchLengths() const;}
//' 
//' @return \code{bool}, \code{true} if the tree has branch lengths.
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
bool HasBranchLengths() const {
    return lengths_.size() == id_parent_.size();
  }

//' @name SPLITT::Tree::LengthOfBranch
//' 
//' @title Get the length of a branch ending at node with id \code{i}.
//' @description
//' \code{\link[=SPLITT::Tree::LengthType]{LengthType} const& \link[=SPLITT::Tree]{Tree}::LengthOfBranch(uint i) const;}
//' 
//' @param i \code{\link[=SPLITT::uint]{uint}}; the id of the end-node for the branch
//' 
//' @return \code{\link[=SPLITT::Tree::LengthType]{LengthType}}; the length associated with the branch ending at node \code{i}.
//'
//' @family public methods in SPLITT::Tree
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
LengthType const& LengthOfBranch(uint i) const {
    if(i >= lengths_.size()) {
      std::ostringstream oss;
      oss<<"ERR:01021:SPLITT:SPLITT.h:LengthOfBranch:: i is beyond the size of the lengths_ vector."<<
        "Check i and that the tree has branches."<<std::endl;
    }
    return lengths_[i];
  }
  
//' @name SPLITT::Tree::BranchLengths
//' 
//' @title Get a reference to the internal vector of branch lengths.
//' @description
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const& BranchLengths() const;}
//' @return \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const&}; 
//'   a const reference to the internally stored vector of branch lengths, in the order of the end-node ids.
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
std::vector<LengthType> const& BranchLengths() const {
    return lengths_;
  }

//' @name SPLITT::Tree::SetLengthOfBranch
//' 
//' @title Set the length of a branch ending at node with id \code{i} to a given \code{value}.
//' @description
//' \code{void SetLengthOfBranch(\link[=SPLITT::uint]{uint} i, \link[=SPLITT::Tree::LengthType]{LengthType} const& value);}
//' @param i \code{\link[=SPLITT::uint]{uint}}; the id of the end-node of the branch;
//' @param value \code{\link[=SPLITT::Tree::LengthType]{LengthType} const&}; the new value to set.
//' @return \code{void}
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
void SetLengthOfBranch(uint i, LengthType const& value) {
    if(!HasBranchLengths()) {
      std::ostringstream oss;
      oss<<"ERR:01031:SPLITT:SPLITT.h:SetLengthOfBranch:: Trying to set a branch length on a tree without branch lengths. "<<
        "Use a SetBranchLengths method to add branch lengths first."<<std::endl;
      throw std::logic_error(oss.str());
    } else if(i >= lengths_.size()) {
      std::ostringstream oss;
      oss<<"i should be smaller than "<<lengths_.size()<<" but was "<<i<<std::endl;
      throw std::out_of_range(oss.str());
    } else {
      lengths_[i] = value;
    }
  }

//' @name SPLITT::Tree::SetBranchLengths
//'
//' @title Set new branch lengths to selected or all branches in a tree. 
//' @description This method has two overwrites:
//' \describe{
//' \item{1. Set the lengths of the branches in the order given by their application-specific end-nodes:}{
//' \code{void SetBranchLengths(}
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const& nodes_branch_ends,} 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const& lengths);}
//' }
//' \item{2. Set a new internally stored vector of branch lengths:}{
//' \code{void SetBranchLengths(}
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const& lengths);}
//' }
//' }
//' 
//' 
//' If the tree has no branch lengths, the supplied arguments should
//'   be of length M-1, where M is the total number of nodes in the tree (-1, 
//'   because there is no branch leading to the root).
//' @param nodes_branch_ends 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const&}: 
//'   a const reference to a new vector of branch lengths, in the order of the 
//'   nodes in \code{nodes_branch_ends}.
//' @param lengths 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::LengthType]{LengthType}> const&}: 
//' \describe{
//' \item{For overwrite 1.}{a const reference to a new vector of branch lengths, 
//' in the order of the nodes in \code{nodes_branch_ends};} 
//' \item{For overwrite 2.}{a const reference to a new vector of branch lengths, 
//' in the order of the end-node ids. In this case (2.), the vector should be of 
//' length M-1, where M is the number of nodes in the tree.}
//' }
//' @return \code{void}
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
void SetBranchLengths(std::vector<NodeType> const& nodes_branch_ends,
                        std::vector<LengthType> const& lengths) {
    if(nodes_branch_ends.size() != lengths.size()) {
      throw std::invalid_argument("ERR:01051:SPLITT:SPLITT.h:SetBranchLengths:: The vectors nodes_branch_ends and lengths should be the same size.");
    }
    if( !HasBranchLengths() ) {
      if(nodes_branch_ends.size() != num_nodes_ - 1) {
        std::ostringstream oss;
        oss<<"ERR:01052:SPLITT:SPLITT.h:SetBranchLengths:: Trying to set branch lengths on a tree without such."<<
          " In this case, the vectors nodes_branch_ends and lengths should have an"<<
            "element for each branch but their size is "<<
              nodes_branch_ends.size()<<" (should be "<<num_nodes_ - 1<<")."<<std::endl;
        throw std::invalid_argument(oss.str());
      }
      lengths_ = vec(num_nodes_ - 1);
    }
    std::vector<bool> visited(num_nodes_ - 1, false);
    for(uint i = 0; i < nodes_branch_ends.size(); ++i) {
      uint id = FindIdOfNode(nodes_branch_ends[i]);
      if(i == G_NA_UINT || i == num_nodes_ - 1) {
        std::ostringstream oss;
        oss<<"ERR:01053:SPLITT:SPLITT.h:SetBranchLengths:: No branch ends at node identified as "<<id<<
          ". Check that nodes_branch_ends correspond to tips or internal nodes (excluding the root)"<<std::endl;
        throw std::logic_error(oss.str());
      } else if(visited[id]) {
        std::ostringstream oss;
        oss<<"ERR:01054:SPLITT:SPLITT.h:SetBranchLengths:: Trying to set the length of the same branch twice. Check nodes_branch_ends for duplicates."<<std::endl;
        throw std::logic_error(oss.str());
      }
      visited[id] = true;
      lengths_[id] = lengths[i];
    }

  }

  void SetBranchLengths(std::vector<LengthType> const& lengths) {
    if(lengths.size() != 0 && lengths.size() != num_nodes_ - 1) {
      std::ostringstream oss;
      oss<<"ERR:01041:SPLITT:SPLITT.h:SetBranchLengths:: lengths should be either empty or of size num_nodes_-1 ("<<
        num_nodes_-1<<") but is "<<lengths.size()<<"."<<std::endl;
    } else {
      lengths_ = lengths;
    }
  }
  
//' @name SPLITT::Tree::FindNodeWithId
//' @title Get the node with the specified id.
//' @description
//' \code{\link[=SPLITT::Tree::NodeType]{NodeType} const& 
//' FindNodeWithId(\link[=SPLITT::uint]{uint} id) const;}
//' 
//' @param id \code{\link[=SPLITT::uint]{uint}} the id of the node (should be 
//' between 0 and M-1, where M is the number of nodes in hte tree).
//' 
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
NodeType const& FindNodeWithId(uint id) const {
    return map_id_to_node_[id];
  }

//' @name SPLITT::Tree::FindIdOfNode
//' @title Get the internally stored id of a node.
//' 
//' @description
//' \code{\link[=SPLITT::uint]{uint} FindIdOfNode(\link[=SPLITT::Tree::NodeType]{NodeType} const& node) const;}
//' 
//' @param node \code{\link[=SPLITT::Tree::NodeType]{NodeType} const&}; the node;
//' 
//' @return an \code{\link[=SPLITT::uint]{uint}} from 0 to M-1, where M is the 
//' number of nodes in the tree.
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
uint FindIdOfNode(NodeType const& node) const {
    auto it = map_node_to_id_.find(node);
    if(it == map_node_to_id_.end()) {
      return G_NA_UINT;
    } else {
      return it->second;
    }
  }

//' @name SPLITT::Tree::FindIdOfParent
//' @title Get the parent id of a node with id id_child.
//' @description
//' \code{uint FindIdOfParent(uint id_child) const;}
//' 
//' @param id_child \code{\link[=SPLITT::uint]{uint}}, the id of the child node. 
//' 
//' @return an \code{\link[=SPLITT::uint]{uint}} from 0 to M-1, where M is the 
//' number of nodes in the tree.
//'   
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
uint FindIdOfParent(uint id_child) const {
    return this->id_parent_[id_child];
  }

//' @name SPLITT::Tree::FindChildren
//' @title Get a vector with the ids of the children of a node.
//' @description
//' \code{\link[=SPLITT::uvec]{uvec} const& FindChildren(uint i) const;}
//' 
//' @param i \code{\link[=SPLITT::uint]{uint}}: the id of a node. 
//' 
//' @return a \code{uvec const&}: a const reference to an internally stored vector 
//' of the ids of the children of \code{i}. If \code{i} is a tip, 
//'   \code{\link[=SPLITT::G_EMPTY_UVEC]{G_EMPTY_UVEC}} is returned. The returned 
//'   vector reference is valid as long as the tree object is not destroyed.
//'   
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
uvec const& FindChildren(uint i) const {
    if(i < this->num_tips()) {
      return G_EMPTY_UVEC;
    } else if(i - this->num_tips() < id_child_nodes_.size()) {
      return id_child_nodes_[i - this->num_tips()];
    } else {
      throw std::invalid_argument("ERR:01061:SPLITT:SPLITT.h:FindChildren:: i must be smaller than the number of nodes.");
    }
  }

//' @name SPLITT::Tree::OrderNodes
//' @title Reorder a vector of nodes
//' @description
//' \code{
//' \link[=SPLITT::uvec]{uvec} OrderNodes(}\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const& nodes) const;}
//' 
//' @param nodes \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const&}; 
//'   a vector of application-specific nodes.
//' @return \code{\link[=SPLITT::uvec]{uvec}}; a vector of positions in \code{nodes} in the order of 
//'   their internally stored ids.
//'   
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
uvec OrderNodes(std::vector<NodeType> const& nodes) const {
    return OrderNodesPosType(nodes, G_NA_UINT);
  }
  

//' @name SPLITT::Tree::OrderNodesPosType
//' @title Reorder a vector of nodes (generic w.r.t. the position type).
//' @section Template Arguments:
//' \describe{
//' \item{PosType}{an integer type for the positions in the returned vector.}}
//' 
//' @description
//' \code{
//' template<class PosType>} 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<PosType> OrderNodesPosType(}\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const& nodes, PosType const& NA) const;}
//' 
//' @param nodes \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<PosType> OrderNodesPosType(}\href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::Tree::NodeType]{NodeType}> const&}; 
//'   a vector of application-specific nodes.
//' @param NA \code{PosType const&}: \code{NA} value used mainly for the purpose of template-specification.
//' 
//' @return \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<PosType>}: 
//' a vector of positions in nodes in the order of their internally stored ids.
//' the element-type of the returned vector can be specified as a template argument 
//' or dedcued from the type of \code{NA}.
//' 
//' @family public methods in SPLITT::Tree 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
template<class PosType>
  std::vector<PosType> OrderNodesPosType(std::vector<NodeType> const& nodes, PosType const& NA) const {
    uvec ids(nodes.size());
    for(uint i = 0; i < nodes.size(); ++i) {
      auto it = this->map_node_to_id_.find(nodes[i]);
      if(it == this->map_node_to_id_.end()) {
        std::ostringstream oss;
        oss<<"ERR:01071:SPLITT:SPLITT.h:OrderNodesPosType:: At least one of the nodes is not present in the tree ("<<
          nodes[i]<<").";
        throw std::invalid_argument(oss.str());
      } else {
        ids[i] = it->second;
      }
    }
    std::vector<PosType> m = Match(Seq(uint(0), this->num_nodes_ - 1), ids, NA);
    return At(m, NotIsNA(m, NA));
  }
};

//' @name SPLITT::OrderedTree
//' 
//' @title A generic tree class optimized for parallel ordered traversal.
//' 
//' @description 
//' This is a template class inheriting from \code{\link[=SPLITT::Tree]{Tree}}.
//' 
//' \code{
//' template<class Node, class Length>class OrderedTree: public Tree<Node, Length>}
//' 
//' @section Template Arguments:
//' \describe{
//' \item{class Node}{see \code{\link[=SPLITT::OrderedTree::NodeType]{NodeType}}.}
//' \item{class Length}{see \code{\link[=SPLITT::OrderedTree::LengthType]{LengthType}}.}
//' }
//' @section Public Methods:
//' \describe{
//' \item{\code{\link[=SPLITT::OrderedTree::num_levels]{num_levels}}}{}
//' \item{\code{\link[=SPLITT::OrderedTree::num_parallel_ranges_prune]{num_parallel_ranges_prune}}}{}
//' \item{\code{\link[=SPLITT::OrderedTree::ranges_id_visit]{ranges_id_visit}}}{}
//' \item{\code{\link[=SPLITT::OrderedTree::RangeIdVisitNode]{RangeIdVisitNode}}}{}
//' \item{\code{\link[=SPLITT::OrderedTree::ranges_id_prune]{ranges_id_prune}}}{}
//' \item{\code{\link[=SPLITT::OrderedTree::RangeIdPruneNode]{RangeIdPruneNode}}}{}
//' }
//' 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
template<class Node, class Length>
class OrderedTree: public Tree<Node, Length> {
public:
//' @name SPLITT::OrderedTree::NodeType 
//' @title Abstract type for nodes in the tree.
//' @description A public typedef in class \code{\link[=SPLITT::OrderedTree]{OrderedTree}}. A synonym for the 
//'   template argument \code{Node}. Defines a 
//'   hash-able type such as \code{int} or 
//'   \href{http://en.cppreference.com/w/cpp/string/basic_string}{\code{std::string}}. 
//'   Specifically, this should be a type for which 
//'   \href{http://en.cppreference.com/w/cpp/utility/hash}{\code{std::hash}}
//'   specialization does exist. This is the application-specific node-type. 
//'   The branches in the tree are defined as couples of nodes 
//'   <branch_start_node, branch_end_node>.
//' @details
//'   During the construction of \code{\link[=SPLITT::OrderedTree]{OrderedTree}} object, the 
//'   nodes are assigned \code{unsigned int} ids from 0 to M-1 (M being the 
//'   number of nodes in the tree). 
//'   
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}}
//' @seealso \link{SPLITT} 
  typedef Node NodeType;
  
//' @name SPLITT::OrderedTree::LengthType
//'
//' @title Abstract type for the branch-lengths in a \code{\link[=SPLITT::OrderedTree]{OrderedTree}}.
//' 
//' @description
//' \code{typedef Length LengthType;}
//' 
//' A public typedef in class \code{\link[=SPLITT::OrderedTree]{OrderedTree}}. 
//'   A synonym for the template argument Length. Defines a type that can be 
//'   associated with a branch. Can be a basic type, e.g. \code{double}, but also
//'   a composite of several attributes on a branch, such as a \code{double} 
//'   length and an \code{int} color.
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}}
//' @seealso \link{SPLITT} 
  typedef Length LengthType;
  
private:
  // default constructor;
  OrderedTree() {}

protected:
  uvec ranges_id_visit_;
  uvec ranges_id_prune_;

public:

//' @name SPLITT::OrderedTree::OrderedTree
//' 
//' @title Constructor for class \code{\link[=SPLITT::OrderedTree]{OrderedTree}}.
//' 
//' @description 
//' \code{
//' OrderedTree(std::vector<NodeType> const& branch_start_nodes,
//' std::vector<NodeType> const& branch_end_nodes,
//' std::vector<LengthType> const& branch_lengths);}
//' 
//' Constructs the tree object given a list of branches. The list of branches
//'   is specified from the corresponding elements in the three vectors passed as
//'   arguments. Creates the internal data-objects needed for ordered traversal 
//'   of the nodes in the tree. 
//' 
//' @param branch_start_nodes 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::OrderedTree::NodeType]{NodeType}> const&}: 
//'   starting node for every branch in the tree.
//' @param branch_end_nodes 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::OrderedTree::NodeType]{NodeType}> const&}:
//'   ending node for every branch in the tree; must be the same length as 
//'   \code{branch_start_nodes}.
//' @param branch_lengths 
//' \href{http://en.cppreference.com/w/cpp/container/vector}{\code{std::vector}}\code{<\link[=SPLITT::OrderedTree::LengthType]{LengthType}> const&}: 
//' lengths associated with the branches. Pass an empty vector for \code{branch_lengths} 
//' for a tree without branch lengths (i.e. only a topology).
//'
//' @family public methods in SPLITT::OrderedTree
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}} \code{\link[=SPLITT::Tree]{Tree}} \code{\link[=SPLITT::Tree::Tree]{Tree::Tree()}}
//' @seealso \link{SPLITT} 
  OrderedTree(
    std::vector<NodeType> const& branch_start_nodes,
    std::vector<NodeType> const& branch_end_nodes,
    std::vector<LengthType> const& branch_lengths):
  Tree<NodeType, LengthType>(branch_start_nodes, branch_end_nodes, branch_lengths),
  ranges_id_visit_(1, 0),
  ranges_id_prune_(1, 0) {

    // insert a fictive branch leading to the root of the tree.
    uvec branch_ends = Seq(uint(0), this->num_nodes_ - 1);

    uvec num_children_remaining(this->num_nodes_, 0);
    for(uint i : this->id_parent_) num_children_remaining[i]++;

    // start by pruning the tips of the tree
    uvec tips_this_level = Seq(uint(0), this->num_tips_ - 1);

    uvec default_pos_vector;
    default_pos_vector.reserve(2);
    std::vector<uvec> pos_of_parent(this->num_nodes_ - this->num_tips_,
                                    default_pos_vector);

    uvec order_branches;
    order_branches.reserve(this->num_nodes_);

    while(tips_this_level[0] != this->num_nodes_ - 1) {
      // while the root has not become a tip itself
      ranges_id_visit_.push_back(
        ranges_id_visit_[ranges_id_visit_.size() - 1] + tips_this_level.size());

      // unique parents at this level
      uvec parents_this_level;
      parents_this_level.reserve(tips_this_level.size());

      uvec tips_next_level;
      tips_next_level.reserve(tips_this_level.size() / 2);

      for(uint i = 0; i < tips_this_level.size(); i++) {
        uint i_parent = this->id_parent_[tips_this_level[i]];
        if(pos_of_parent[i_parent - this->num_tips_].empty()) {
          parents_this_level.push_back(i_parent);
        }
        pos_of_parent[i_parent - this->num_tips_].push_back(i);
      }

      uint num_parents_remaining = parents_this_level.size();
      while( num_parents_remaining ) {

        uint num_parent_updates = 0;
        for(auto i_parent: parents_this_level) {
          if(!pos_of_parent[i_parent - this->num_tips_].empty()) {
            uint i = pos_of_parent[i_parent - this->num_tips_].back();

            num_parent_updates ++;
            order_branches.push_back(tips_this_level[i]);

            num_children_remaining[i_parent]--;
            if(num_children_remaining[i_parent] == 0) {
              tips_next_level.push_back(i_parent);
            }
            pos_of_parent[i_parent - this->num_tips_].pop_back();
            if(pos_of_parent[i_parent - this->num_tips_].empty()) {
              num_parents_remaining--;
            }
          }
        }

        ranges_id_prune_.push_back(
          ranges_id_prune_[ranges_id_prune_.size() - 1] + num_parent_updates);
      }

      tips_this_level = tips_next_level;
    }

    if(this->HasBranchLengths()) {
      this->lengths_ = At(this->lengths_, order_branches);
    }

    uvec id_old = order_branches;
    id_old.push_back(this->num_nodes_ - 1);

    this->id_parent_ = Match(At(this->id_parent_, order_branches),
                             id_old, G_NA_UINT);

    // update maps
    std::vector<NodeType> map_id_to_node(this->num_nodes_);
    for (uint i = 0; i < this->num_nodes_; i++) {
      map_id_to_node[i] = this->map_id_to_node_[id_old[i]];
      this->map_node_to_id_[map_id_to_node[i]] = i;
    }

    std::swap(this->map_id_to_node_, map_id_to_node);

    this->init_id_child_nodes();
  }

//' @name SPLITT::OrderedTree::num_levels
//' 
//' @title Number of levels (ranges) of parallel \code{VisitNode} operations 
//' during post-order traversal.
//' 
//' @description 
//' \code{
//' \link[=SPLITT::uint]{uint} num_levels() const;}
//' 
//' During range-based post-order traversal, levels represent groups of nodes
//' that can be visited independent from one another and, therefore, in parallel.
//' In a balanced tree the number of levels is in the order of O(log2(N)), where N
//' is the number of tips in the tree. Hence, parallelization can be efficient.
//' In a strongly unbalanced tree, the number of levels is in the order of O(N), 
//' and the parallelization of the \code{VisitNode} operation cannot improve the 
//' speed. 
//' 
//' @family public methods in SPLITT::OrderedTree
//' @seealso \code{\link[=SPLITT::Tree]{Tree}} 
//'  \code{\link[=SPLITT::OrderedTree]{OrderedTree}}
//'  \code{\link[=SPLITT::OrderedTree::num_parallel_ranges_prune]{num_parallel_ranges_prune}}
//' @seealso \link{SPLITT} 
  uint num_levels() const {
    return ranges_id_visit_.size() - 1;
  }

//' @name SPLITT::OrderedTree::num_parallel_ranges_prune
//'  
//' @title Number of parallel ranges of nodeduring  
//' 
//' @description 
//' \code{
//' \link[=SPLITT::uint]{uint} num_levels() const;}
//' 
//' During post-order traversal, a node is pruned from its parent after it gets 
//' visited. This is the so called operation \code{PruneNode}. Prune-ranges 
//' represent groups of nodes that can be pruned independent from one another 
//' and, therefore, in parallel. Conceptually prune-ranges are similar to levels 
//' of (parallel) \code{VisitNode} operations (see also
//' \link[=SPLITT::OrderedTree::num_levels]{num_levels}()).
//' In a balanced tree the number of levels is in the order of O(log2(N)), where N
//' is the number of tips in the tree and, therefore, parallelization can be 
//' beneficial. In a strongly unbalanced tree, the number of levels is in the 
//' order of O(N), and the parallelization of the \code{VisitNode} operation
//' cannot improve the speed. 
//' 
//' @family public methods in SPLITT::OrderedTree
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}} \code{\link[=SPLITT::Tree]{Tree}} \code{\link[=SPLITT::Tree::Tree]{Tree::Tree()}}
//' @seealso \link{SPLITT} 
  uint num_parallel_ranges_prune() const {
    return ranges_id_prune_.size() - 1;
  }

//' @name SPLITT::OrderedTree::ranges_id_visit
//'  
//' @title Internally stored vector of start ids for each VisitNode-range 
//' 
//' @description 
//' \code{
//' \link[=SPLITT::uvec]{uvec} const& ranges_id_visit() const;}
//' 
//' @family public methods in SPLITT::OrderedTree
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}} 
//' @seealso \link{SPLITT} 
  uvec const& ranges_id_visit() const {
    return ranges_id_visit_;
  }

//' @name SPLITT::OrderedTree::RangeIdVisitNode
//'  
//' @title First and last id (0-based node indices) of the VisitNode-range at a level
//' 
//' @param i_level the level (0-based) of the VisitNode-range
//' 
//' @description 
//' \href{http://en.cppreference.com/w/cpp/container/array}{\code{std::array}}\code{<\link[=SPLITT::uint]{uint}, 2>}\code{ RangeIdVisitNode(\link[=SPLITT::uint]{uint} i_level) const;}
//' 
//' @family public methods in SPLITT::OrderedTree
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}} 
//' @seealso \link{SPLITT} 
  std::array<uint, 2> RangeIdVisitNode(uint i_level) const {
    // double braces required by C++11 standard,
    // http://en.cppreference.com/w/cpp/container/array
    return std::array<uint, 2> {{ranges_id_visit_[i_level],
                                 ranges_id_visit_[i_level+1] - 1}};
  }

//' @name SPLITT::OrderedTree::ranges_id_prune
//'  
//' @title Internally stored vector of start ids for each PruneNode-range 
//' 
//' @description 
//' \code{
//' \link[=SPLITT::uvec]{uvec} const& ranges_id_prune() const;}
//' 
//' @family public methods in SPLITT::OrderedTree
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}} 
//' @seealso \link{SPLITT} 
  uvec const& ranges_id_prune() const {
    return ranges_id_prune_;
  }

  
//' @name SPLITT::OrderedTree::RangeIdPruneNode
//'  
//' @title First and last id (0-based node indices) of the PruneNode-range at a step
//' 
//' @param i_step the step (0-based) of the PruneNode-range
//' 
//' @description 
//' \href{http://en.cppreference.com/w/cpp/container/array}{\code{std::array}}\code{<\link[=SPLITT::uint]{uint}, 2>}\code{ RangeIdPruneNode(\link[=SPLITT::uint]{uint} i_step) const;}
//' 
//' @family public methods in SPLITT::OrderedTree
//' @seealso \code{\link[=SPLITT::OrderedTree]{OrderedTree}} 
//' @seealso \link{SPLITT} 
  std::array<uint, 2> RangeIdPruneNode(uint i_step) const {
    return std::array<uint, 2> {{ranges_id_prune_[i_step],
                                 ranges_id_prune_[i_step+1] - 1}};
  }
};

enum PostOrderMode {
  AUTO = 0,
  SINGLE_THREAD_LOOP_POSTORDER = 10,
  SINGLE_THREAD_LOOP_PRUNES = 11,
  SINGLE_THREAD_LOOP_VISITS = 12,
  MULTI_THREAD_LOOP_PRUNES = 21,
  MULTI_THREAD_LOOP_VISITS = 22,
  MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES = 23,
  MULTI_THREAD_VISIT_QUEUE = 24,
  MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION = 25,
  HYBRID_LOOP_PRUNES = 31,
  HYBRID_LOOP_VISITS = 32,
  HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES = 33
};

inline std::ostream& operator<< (std::ostream& os, PostOrderMode mode) {
  switch(mode) {
  case PostOrderMode::AUTO: os<<"AUTO"; break;
  case PostOrderMode::SINGLE_THREAD_LOOP_POSTORDER: os<<"SINGLE_THREAD_LOOP_POSTORDER"; break;
  case PostOrderMode::SINGLE_THREAD_LOOP_PRUNES: os<<"SINGLE_THREAD_LOOP_PRUNES"; break;
  case PostOrderMode::SINGLE_THREAD_LOOP_VISITS: os<<"SINGLE_THREAD_LOOP_VISITS"; break;
  case PostOrderMode::MULTI_THREAD_LOOP_PRUNES: os<<"MULTI_THREAD_LOOP_PRUNES"; break;
  case PostOrderMode::MULTI_THREAD_LOOP_VISITS: os<<"MULTI_THREAD_LOOP_VISITS"; break;
  case PostOrderMode::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES: os<<"MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES"; break;
  case PostOrderMode::MULTI_THREAD_VISIT_QUEUE: os<<"MULTI_THREAD_VISIT_QUEUE"; break;
  case PostOrderMode::MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION: os<<"MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION"; break;
  case PostOrderMode::HYBRID_LOOP_PRUNES: os<<"HYBRID_LOOP_PRUNES"; break;
  case PostOrderMode::HYBRID_LOOP_VISITS: os<<"HYBRID_LOOP_VISITS"; break;
  case PostOrderMode::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES: os<<"HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES"; break;
  };
  return os<< static_cast<int>(mode);
}

template<class TreeType> class VisitQueue {
  std::mutex mutex_;
  std::condition_variable has_a_new_node_;
  
  TreeType const& ref_tree_;
  uvec queue_;
  uvec::iterator it_queue_begin;
  uvec::iterator it_queue_end;
  uvec num_non_visited_children_;
public:
  
  // non-thread safe (called in single-thread mode)
  void Init(uvec const& num_children) {
    std::copy(num_children.begin(), num_children.end(),
              num_non_visited_children_.begin());
    it_queue_begin = queue_.begin();
    it_queue_end = queue_.begin() + ref_tree_.num_tips();
    std::iota(it_queue_begin, it_queue_end, 0);
  }
  
  bool IsTemporarilyEmpty() const {
    return it_queue_begin == it_queue_end && it_queue_end < queue_.end();
  }
  
  // thread-safe
  uint NextInQueue() {
    std::unique_lock<std::mutex> lock(mutex_);
    
    while( IsTemporarilyEmpty() ) {
      has_a_new_node_.wait(lock);
    }
    
    if(it_queue_begin < it_queue_end) {
      uint res = *it_queue_begin;
      ++it_queue_begin;
      return res;
    } else if(it_queue_begin == queue_.end()) {
      // algorithm thread should stop here. all waiting threads should be notified,
      // since no other elements will be inserted in the queue.
      has_a_new_node_.notify_all();
      return ref_tree_.num_nodes();
    } else {
      // algorithm thread continues to check for new node to visit
      // should never execute this
      //std::cout<<"Error returning G_NA_UINT from VisitQueue."<<std::endl;
      return G_NA_UINT;
    }
  }
  
  // thread-safe
  // if the parent of i becomes visit-able, it gets inserted in the
  // queue.
  void RemoveVisitedNode(uint i) {
    std::unique_lock<std::mutex> lock(mutex_);
    
    uint i_parent = ref_tree_.FindIdOfParent(i);
    num_non_visited_children_[i_parent - ref_tree_.num_tips()]--;
    if(num_non_visited_children_[i_parent - ref_tree_.num_tips()] == 0) {
      *it_queue_end = i_parent;
      *it_queue_end++;
      has_a_new_node_.notify_one();
    }
  }
  
  // non-thread-safe. should call Init() before using.
  VisitQueue(TreeType const& tree):
  ref_tree_(tree),
  queue_(tree.num_nodes()),
  it_queue_begin(queue_.begin()),
  it_queue_end(queue_.begin()),
  num_non_visited_children_(tree.num_nodes() - tree.num_tips()) {}
  
  // Copy initialization (non-thread-safe)
  VisitQueue(const VisitQueue& other): ref_tree_(other.ref_tree_) {
    auto other_begin = other.queue_.begin();
    queue_ = other.queue_;
    it_queue_begin = queue_.begin() + (other.it_queue_begin  - other_begin);
    it_queue_end = queue_.begin() + (other.it_queue_end  - other_begin);
    num_non_visited_children_ = other.num_non_visited_children_;
  }
};

//' @name SPLITT::TraversalAlgorithm
//' 
//' @title Base-class for parallel tree traversal implementations.
//' 
//' @description 
//' This is a template class inheriting from \code{\link[=SPLITT::Tree]{Tree}}.
//' 
//' \code{
//' template<class Node, class Length>class OrderedTree: public Tree<Node, Length>}
//' 
//' @section Template Arguments:
//' \describe{
//' \item{class Node}{see \code{\link[=SPLITT::OrderedTree::NodeType]{NodeType}}.}
//' \item{class Length}{see \code{\link[=SPLITT::OrderedTree::LengthType]{LengthType}}.}
//' }
//' @section Public Methods:
//' \describe{
//' \item{\code{TraversalAlgorithm}}{}
//' }
//' 
//' @seealso \code{\link[=SPLITT::Tree]{Tree}}
//' @seealso \link{SPLITT} 
template<class TraversalSpecification>
class TraversalAlgorithm {
  
public:
  typedef typename TraversalSpecification::TreeType TreeType;

  TreeType const& ref_tree_;
  TraversalSpecification& ref_spec_;

  uvec num_children_;
  VisitQueue<TreeType> visit_queue_;

public:
  TraversalAlgorithm(TreeType const& tree, TraversalSpecification& spec):
  ref_tree_(tree),
  ref_spec_(spec),
  num_children_(tree.num_nodes() - tree.num_tips()),
  visit_queue_(tree) {
    for(uint i = tree.num_tips(); i < tree.num_nodes(); i++) {
      num_children_[i - tree.num_tips()] = tree.FindChildren(i).size();
    }
  }

  uint NumOmpThreads() const {
#ifdef _OPENMP
    return omp_get_max_threads();
#else 
    return 1;
#endif  // #ifdef _OPENMP
  }

  uint VersionOPENMP() const {
#ifdef _OPENMP
    return _OPENMP;
#else
    return 0;
#endif
  }
};

//' @name SPLITT::ThreadExceptionHandler
//' 
//' @title An internal class for thread-safe exception hadling within parallel sections.
//' 
//' @description This class is used as a wrapper for InitNode,VistNode and PruneNode 
//' function calls within parallel TraverseTree executions. It is inspired from this 
//' \href{https://stackoverflow.com/questions/11828539/elegant-exceptionhandling-in-openmp}{stackoverflow discussion}. 
//' Synopsis:
//' \code{class ThreadExceptionHandler;}
//' 
//' @section Public Methods:
//' \describe{
//' \item{\link[=SPLITT::ThreadExceptionHandler::ThreadExceptionHandler]{\code{ThreadExceptionHandler}}}{}
//' \item{\link[=SPLITT::ThreadExceptionHandler::Run]{\code{Run}}}{}
//' \item{\link[=SPLITT::ThreadExceptionHandler::CaptureException]{\code{CaptureException}}}{}
//' \item{\link[=SPLITT::ThreadExceptionHandler::Rethrow]{\code{Rethrow}}}{}
//' }
//' @seealso \link[=SPLITT::TraversalAlgorithm]{TraversalAlgorithm}
//' @seealso \link[=SPLITT::PostOrderTraversal]{PostOrderTraversal}
//' @seealso \link[=SPLITT::PreOrderTraversal]{PreOrderTraversal}
//' @seealso \link{SPLITT} 
class ThreadExceptionHandler {
  std::exception_ptr ptr_;
  std::mutex         lock_;
public:
  ThreadExceptionHandler(): ptr_(nullptr) {}
  // Copy initialization (non-thread-safe)
  ThreadExceptionHandler(const ThreadExceptionHandler& other): ptr_(other.ptr_) {}
  
  template <typename Function, typename... Parameters>
  void Run(Function f, Parameters... params) {
    try {
      f(params...);
    }
    catch (...) {
      CaptureException();
    }
  }
  
  void CaptureException() { 
    std::unique_lock<std::mutex> guard(this->lock_);
    this->ptr_ = std::current_exception(); 
  }
  
  // If there is an exception it gets rethrown and the ptr_ is set to nullptr.
  void Rethrow() {
    if(this->ptr_) {
      std::exception_ptr ptr = this->ptr_;
      // reset ptr_ to nullptr so the handler object is cleaned and reusable
      // after the call to Rethrow().
      this->ptr_ = nullptr;
      std::rethrow_exception(ptr);
    }
  }
};

template<class TraversalSpecification>
class PostOrderTraversal: public TraversalAlgorithm<TraversalSpecification> {
  
  ThreadExceptionHandler exception_handler_;  
  
public:
  typedef TraversalAlgorithm<TraversalSpecification> ParentType;

  typedef PostOrderMode ModeType;

  PostOrderTraversal(typename TraversalSpecification::TreeType const& tree,
                     TraversalSpecification& spec): ParentType(tree, spec) { }

  void TraverseTree(ModeType mode) {
    switch(mode) {
    case ModeType::SINGLE_THREAD_LOOP_POSTORDER: TraverseTreeSingleThreadLoopPostorder(); break;
    case ModeType::SINGLE_THREAD_LOOP_PRUNES: TraverseTreeSingleThreadLoopPrunes(); break;
    case ModeType::SINGLE_THREAD_LOOP_VISITS: TraverseTreeSingleThreadLoopVisits(); break;
    case ModeType::MULTI_THREAD_LOOP_PRUNES: TraverseTreeMultiThreadLoopPrunes(); break;
    case ModeType::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES: TraverseTreeMultiThreadLoopVisitsThenLoopPrunes(); break;
    case ModeType::MULTI_THREAD_LOOP_VISITS: TraverseTreeMultiThreadLoopVisits(); break;
    case ModeType::MULTI_THREAD_VISIT_QUEUE: TraverseTreeMultiThreadVisitQueue(); break;
    case ModeType::MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION: TraverseTreeMultiThreadLoopPrunesNoException(); break;
    case ModeType::HYBRID_LOOP_PRUNES: TraverseTreeHybridLoopPrunes(); break;
    case ModeType::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES: TraverseTreeHybridLoopVisitsThenLoopPrunes(); break;
    case ModeType::HYBRID_LOOP_VISITS: TraverseTreeHybridLoopVisits(); break;
    default: TraverseTreeAuto();
    }
    exception_handler_.Rethrow();
  }
protected:
  uint current_step_tuning_ = 0;
  uint fastest_step_tuning_ = 0;

  double min_duration_tuning_ = std::numeric_limits<double>::max();
  std::vector<double> durations_tuning_;

  const uvec min_sizes_chunk_ = {8}; //, 4, 8, 16, 32};

  const std::vector<ModeType> choices_mode_auto_ = {
    ModeType::SINGLE_THREAD_LOOP_POSTORDER,
    ModeType::SINGLE_THREAD_LOOP_PRUNES,
    ModeType::SINGLE_THREAD_LOOP_VISITS,
    ModeType::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES,
    ModeType::MULTI_THREAD_LOOP_VISITS,
    ModeType::MULTI_THREAD_VISIT_QUEUE
  };

  const std::vector<ModeType> choices_hybrid_mode_auto_ = {
    ModeType::HYBRID_LOOP_PRUNES,
    ModeType::HYBRID_LOOP_VISITS,
    ModeType::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES
  };

public:
  bool IsTuning() const {
    return current_step_tuning_ < choices_mode_auto_.size() +
      min_sizes_chunk_.size() * choices_hybrid_mode_auto_.size();
  }


  std::string ModeAutoCurrent() const {
    std::ostringstream oss;
    oss<<ModeAuto();
    return oss.str();
  }

  std::string ModeAutoStep(uint step) const {
    std::ostringstream oss;
    oss<<ModeAuto(step);
    return oss.str();
  }

  ModeType ModeAuto() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return ModeAuto(step);
  }

  ModeType ModeAuto(uint step) const {
    if( step < choices_mode_auto_.size() ) {
      return choices_mode_auto_[step];
    } else {
      uint k = choices_hybrid_mode_auto_.size();
      uint l = step - choices_mode_auto_.size();
      return choices_hybrid_mode_auto_[(l/k) % k];
    }

  }

  uint IndexMinSizeChunkVisit() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    //return (step / min_sizes_chunk_.size()) % min_sizes_chunk_.size();
    return step % min_sizes_chunk_.size();
  }

  uint IndexMinSizeChunkPrune() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return step % min_sizes_chunk_.size();
  }

  uint min_size_chunk_visit() const {
    return min_sizes_chunk_[IndexMinSizeChunkVisit()];
  }

  uint min_size_chunk_prune() const {
    return min_sizes_chunk_[IndexMinSizeChunkPrune()];
  }

  uint fastest_step_tuning() const {
    return fastest_step_tuning_;
  }

  std::vector<double>  durations_tuning() const {
    return durations_tuning_;
  }

protected:
  void TraverseTreeAuto() {

    std::chrono::steady_clock::time_point start, end;
    double duration;

    ModeType mode = ModeAuto();

    if( IsTuning() ) {

      start = std::chrono::steady_clock::now();
      TraverseTree(mode);
      end = std::chrono::steady_clock::now();

      duration = std::chrono::duration<double, std::milli>(end - start).count();
      durations_tuning_.push_back(duration);
      if(duration < min_duration_tuning_) {
        min_duration_tuning_ = duration;
        fastest_step_tuning_ = current_step_tuning_;
      }
      current_step_tuning_++;

    } else {
      TraverseTree(mode);
    }
  }
  void TraverseTreeSingleThreadLoopPostorder() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }
    exception_handler_.Rethrow();

    for(uint i = 0; i < ParentType::ref_tree_.num_nodes() - 1; i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.VisitNode(i);
        ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
      });
    }
    exception_handler_.Rethrow();
  }

  void TraverseTreeSingleThreadLoopPrunes() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }
    exception_handler_.Rethrow();

    for(uint i_prune = 0;
        i_prune < ParentType::ref_tree_.num_parallel_ranges_prune();
        i_prune++) {
      auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);

    _PRAGMA_OMP_SIMD
      for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
        exception_handler_.Run([=]{
          ParentType::ref_spec_.VisitNode(i);
          ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
        });
      }
      exception_handler_.Rethrow();
    }
  }

  void TraverseTreeSingleThreadLoopVisits() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }
    exception_handler_.Rethrow();

    for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          if(i < ParentType::ref_tree_.num_tips()) {
            // i is a tip (only Visit)
            ParentType::ref_spec_.VisitNode(i);
          } else {
            // i is internal
            for(uint j: ParentType::ref_tree_.FindChildren(i)) {
              ParentType::ref_spec_.PruneNode(j, i);
            }
            ParentType::ref_spec_.VisitNode(i);
          }
        });
      }
      exception_handler_.Rethrow();
    }

    // VisitNode not called on the root node
    for(uint j: ParentType::ref_tree_.FindChildren(ParentType::ref_tree_.num_nodes() - 1)) {
      ParentType::ref_spec_.PruneNode(j, ParentType::ref_tree_.num_nodes() - 1);
    }
  }

  void TraverseTreeMultiThreadLoopVisitsThenLoopPrunes() {

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
    exception_handler_.Run([=]{
      ParentType::ref_spec_.InitNode(i);  
    });
  }
  exception_handler_.Rethrow();

  uint i_prune = 0;
  for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {

#pragma omp barrier

    auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          ParentType::ref_spec_.VisitNode(i);
        });
      }
      exception_handler_.Rethrow();

      uint num_branches_done = 0;

    while(num_branches_done != range_visit[1] - range_visit[0] + 1) {
#pragma omp barrier
      auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);

      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
          });
        }
        exception_handler_.Rethrow();

        num_branches_done +=  range_prune[1] - range_prune[0] + 1;
      ++i_prune;
    }
  }
}
  }

  void TraverseTreeMultiThreadLoopVisits() {
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }
    exception_handler_.Rethrow();

    for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          if(i < ParentType::ref_tree_.num_tips()) {
            // i is a tip (only Visit)
            ParentType::ref_spec_.VisitNode(i);
          } else {
            // i is internal
            for(uint j: ParentType::ref_tree_.FindChildren(i)) {
              ParentType::ref_spec_.PruneNode(j, i);
            }
            ParentType::ref_spec_.VisitNode(i);
          }
        });
      }
      exception_handler_.Rethrow();
    }
}
    // VisitNode not called on the root node
    for(uint j: ParentType::ref_tree_.FindChildren(ParentType::ref_tree_.num_nodes() - 1)) {
      ParentType::ref_spec_.PruneNode(j, ParentType::ref_tree_.num_nodes() - 1);
    }
  }

  void TraverseTreeMultiThreadVisitQueue() {
    ParentType::visit_queue_.Init(ParentType::num_children_);
#pragma omp parallel
{
  exception_handler_.Run([=]{
    while(true) {
      uint i = ParentType::visit_queue_.NextInQueue();
      if(i == G_NA_UINT) {
        continue;
      } else if(i == ParentType::ref_tree_.num_nodes()) {
        break;
      } else if(i < ParentType::ref_tree_.num_tips()) {
        // i is a tip (only Visit)
        ParentType::ref_spec_.InitNode(i);
        ParentType::ref_spec_.VisitNode(i);
        ParentType::visit_queue_.RemoveVisitedNode(i);
      } else if(i < ParentType::ref_tree_.num_nodes() - 1){
        // i is internal
        ParentType::ref_spec_.InitNode(i);
        uvec const& children = ParentType::ref_tree_.FindChildren(i);
        for(uint j: children) {
          ParentType::ref_spec_.PruneNode(j, i);
        }
        ParentType::ref_spec_.VisitNode(i);
        ParentType::visit_queue_.RemoveVisitedNode(i);
      } else {
        // i is the root
        ParentType::ref_spec_.InitNode(i);
        uvec const& children = ParentType::ref_tree_.FindChildren(i);
        for(uint j: children) {
          ParentType::ref_spec_.PruneNode(j, i);
        }
        // don't visit the root
      }
    }
  });
}
exception_handler_.Rethrow();
}

  void TraverseTreeMultiThreadLoopPrunes() {

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
    exception_handler_.Run([=]{
      ParentType::ref_spec_.InitNode(i);
    });
  }
  exception_handler_.Rethrow();

  for(uint i_prune = 0; i_prune < ParentType::ref_tree_.num_parallel_ranges_prune(); i_prune++) {
    auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);

    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
        exception_handler_.Run([=]{
          ParentType::ref_spec_.VisitNode(i);
          ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
        });
      }
      exception_handler_.Rethrow();
  }
}
  }

  void TraverseTreeMultiThreadLoopPrunesNoException() {
    
#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
    ParentType::ref_spec_.InitNode(i);
  }
  
  for(uint i_prune = 0; i_prune < ParentType::ref_tree_.num_parallel_ranges_prune(); i_prune++) {
    auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);
    
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
        ParentType::ref_spec_.VisitNode(i);
        ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
      }
  }
}
  }
  
  void TraverseTreeHybridLoopVisitsThenLoopPrunes() {
    uint min_size_chunk_visit = this->min_size_chunk_visit();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }
    exception_handler_.Rethrow();

    uint i_prune = 0;
  for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
    auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
#pragma omp barrier
    if(range_visit[1] - range_visit[0] + 1 >
         ParentType::NumOmpThreads() * min_size_chunk_visit) {
      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.VisitNode(i);
          });
        }
        exception_handler_.Rethrow();
    } else if(tid == 0) {
      // only the master thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          ParentType::ref_spec_.VisitNode(i);
        });
      }
      exception_handler_.Rethrow();
    }

    if (tid == 0) {
      // only one (master) thread executes this
      uint num_branches_done = 0;
      while(num_branches_done != range_visit[1] - range_visit[0] + 1) {
        auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);
        _PRAGMA_OMP_SIMD
          for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
            exception_handler_.Run([=]{
              ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
            });
          }
          exception_handler_.Rethrow();

          num_branches_done +=  range_prune[1] - range_prune[0] + 1;
        ++i_prune;
      }
    }
  }
}
  }

  void TraverseTreeHybridLoopPrunes() {
    uint min_size_chunk_prune = this->min_size_chunk_prune();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }
    exception_handler_.Rethrow();


  for(uint i_prune = 0; i_prune < ParentType::ref_tree_.num_parallel_ranges_prune(); i_prune++) {
      auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);
#pragma omp barrier
      if (range_prune[1] - range_prune[0] + 1 >
            ParentType::NumOmpThreads() * min_size_chunk_prune) {
        _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.VisitNode(i);
            ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
          });
        }
        exception_handler_.Rethrow();
      } else if (tid == 0) {
        // only one (master) thread executes this
        _PRAGMA_OMP_SIMD
        for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.VisitNode(i);
            ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
          });
        }
        exception_handler_.Rethrow();
      }
    }
}
  }

  void TraverseTreeHybridLoopVisits() {
    uint min_size_chunk_visit = this->min_size_chunk_visit();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }
    exception_handler_.Rethrow();

  for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
    auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
#pragma omp barrier
    if(range_visit[1] - range_visit[0] + 1 >
         ParentType::NumOmpThreads() * min_size_chunk_visit) {
      _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          if(i < ParentType::ref_tree_.num_tips()) {
            // i is a tip (only Visit)
            ParentType::ref_spec_.VisitNode(i);
          } else if(i < ParentType::ref_tree_.num_nodes() - 1){
            // i is internal
            for(uint j: ParentType::ref_tree_.FindChildren(i)) {
              ParentType::ref_spec_.PruneNode(j, i);
            }
            ParentType::ref_spec_.VisitNode(i);
          }
        });
      }
      exception_handler_.Rethrow();
    } else if(tid == 0) {
      // only the master thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          if(i < ParentType::ref_tree_.num_tips()) {
            // i is a tip (only Visit)
            ParentType::ref_spec_.VisitNode(i);
          } else if(i < ParentType::ref_tree_.num_nodes() - 1){
            // i is internal
            for(uint j: ParentType::ref_tree_.FindChildren(i)) {
              ParentType::ref_spec_.PruneNode(j, i);
            }
            ParentType::ref_spec_.VisitNode(i);
          }
        });
      }
      exception_handler_.Rethrow();
    }
  }
}
    // VisitNode not called on the root
    for(uint j: ParentType::ref_tree_.FindChildren(ParentType::ref_tree_.num_nodes() - 1)) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.PruneNode(j, ParentType::ref_tree_.num_nodes() - 1);
      });
    }
    exception_handler_.Rethrow();
  }
};


enum PreOrderMode {
  PREORDER_AUTO = 0,
  PREORDER_SINGLE_THREAD_LOOP_PREORDER = 10,
  PREORDER_SINGLE_THREAD_LOOP_VISITS = 12,
  PREORDER_MULTI_THREAD_LOOP_VISITS = 22
};

inline std::ostream& operator<< (std::ostream& os, PreOrderMode mode) {
  switch(mode) {
  case PreOrderMode::PREORDER_AUTO: os<<"PREORDER_AUTO"; break;
  case PreOrderMode::PREORDER_SINGLE_THREAD_LOOP_PREORDER: os<<"PREORDER_SINGLE_THREAD_LOOP_PREORDER"; break;
  case PreOrderMode::PREORDER_SINGLE_THREAD_LOOP_VISITS: os<<"PREORDER_SINGLE_THREAD_LOOP_VISITS"; break;
  case PreOrderMode::PREORDER_MULTI_THREAD_LOOP_VISITS: os<<"PREORDER_MULTI_THREAD_LOOP_VISITS"; break;
  };
  return os<< static_cast<int>(mode);
}

template<class TraversalSpecification>
class PreOrderTraversal: public TraversalAlgorithm<TraversalSpecification> {

  typedef TraversalAlgorithm<TraversalSpecification> ParentType;

public:
  typedef PreOrderMode ModeType;

  PreOrderTraversal(typename TraversalSpecification::TreeType const& tree,
                     TraversalSpecification& spec): ParentType(tree, spec) { }

  void TraverseTree(ModeType mode) {
    switch(mode) {
    case ModeType::PREORDER_SINGLE_THREAD_LOOP_PREORDER: TraverseTreeSingleThreadLoopPreorder(); break;
    case ModeType::PREORDER_SINGLE_THREAD_LOOP_VISITS: TraverseTreeSingleThreadLoopVisits(); break;
    case ModeType::PREORDER_MULTI_THREAD_LOOP_VISITS: TraverseTreeMultiThreadLoopVisits(); break;
    default: TraverseTreeAuto();
    }
  }
protected:
  uint current_step_tuning_ = 0;
  uint fastest_step_tuning_ = 0;

  double min_duration_tuning_ = std::numeric_limits<double>::max();
  std::vector<double> durations_tuning_;

  const uvec min_sizes_chunk_ = {8}; //, 4, 8, 16, 32};

  const std::vector<ModeType> choices_mode_auto_ = {
    ModeType::PREORDER_SINGLE_THREAD_LOOP_PREORDER,
    ModeType::PREORDER_SINGLE_THREAD_LOOP_VISITS,
    ModeType::PREORDER_MULTI_THREAD_LOOP_VISITS
  };

public:
  bool IsTuning() const {
    return current_step_tuning_ < choices_mode_auto_.size();
  }


  std::string ModeAutoCurrent() const {
    std::ostringstream oss;
    oss<<ModeAuto();
    return oss.str();
  }

  std::string ModeAutoStep(uint step) const {
    std::ostringstream oss;
    oss<<ModeAuto(step);
    return oss.str();
  }

  ModeType ModeAuto() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return ModeAuto(step);
  }

  ModeType ModeAuto(uint step) const {
    return choices_mode_auto_[step%choices_mode_auto_.size()];
  }

  uint fastest_step_tuning() const {
    return fastest_step_tuning_;
  }

  std::vector<double>  durations_tuning() const {
    return durations_tuning_;
  }

protected:
  void TraverseTreeAuto() {

    std::chrono::steady_clock::time_point start, end;
    double duration;

    ModeType mode = ModeAuto();

    if( IsTuning() ) {

      start = std::chrono::steady_clock::now();
      TraverseTree(mode);
      end = std::chrono::steady_clock::now();

      duration = std::chrono::duration<double, std::milli>(end - start).count();
      durations_tuning_.push_back(duration);
      if(duration < min_duration_tuning_) {
        min_duration_tuning_ = duration;
        fastest_step_tuning_ = current_step_tuning_;
      }
      current_step_tuning_++;

    } else {
      TraverseTree(mode);
    }
  }

  void TraverseTreeSingleThreadLoopPreorder() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    for(uint i = ParentType::ref_tree_.num_nodes() - 1; ; i--) {
      ParentType::ref_spec_.VisitNode(i);
      if(i == 0) {
        break;
      }
    }
  }

  void TraverseTreeSingleThreadLoopVisits() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    ParentType::ref_spec_.VisitNode(ParentType::ref_tree_.num_nodes() - 1);

    for(uint i_level = ParentType::ref_tree_.num_levels(); i_level > 0; i_level--) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level - 1);
      _PRAGMA_OMP_SIMD
        for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
          ParentType::ref_spec_.VisitNode(i);
        }
    }
  }

  void TraverseTreeMultiThreadLoopVisits() {
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    ParentType::ref_spec_.VisitNode(ParentType::ref_tree_.num_nodes() - 1);

    for(uint i_level = ParentType::ref_tree_.num_levels(); i_level > 0; i_level--) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level - 1);
      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
          ParentType::ref_spec_.VisitNode(i);
        }
    }
}
  }

};

}
#endif // SPLITT_SPLITT_H_
