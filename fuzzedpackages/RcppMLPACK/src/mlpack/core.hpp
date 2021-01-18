/***
 * @file core.hpp
 *
 * Include all of the base components required to write MLPACK methods, and the
 * main MLPACK Doxygen documentation.
 *
 * This file is part of MLPACK 1.0.10.
 *
 * MLPACK is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * MLPACK is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details (LICENSE.txt).
 *
 * You should have received a copy of the GNU General Public License along with
 * MLPACK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __MLPACK_CORE_HPP
#define __MLPACK_CORE_HPP

#if _WIN64
#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif
#endif

/**
 * @mainpage MLPACK Documentation
 *
 * @section intro_sec Introduction
 *
 * MLPACK is an intuitive, fast, scalable C++ machine learning library, meant to
 * be a machine learning analog to LAPACK.  It aims to implement a wide array of
 * machine learning methods and function as a "swiss army knife" for machine
 * learning researchers.  The MLPACK development website can be found at
 * http://mlpack.org.
 *
 * MLPACK uses the Armadillo C++ matrix library (http://arma.sourceforge.net)
 * for general matrix, vector, and linear algebra support.  MLPACK also uses the
 * program_options, math_c99, and unit_test_framework components of the Boost
 * library; in addition, LibXml2 is used.
 *
 * @section howto How To Use This Documentation
 *
 * This documentation is API documentation similar to Javadoc.  It isn't
 * necessarily a tutorial, but it does provide detailed documentation on every
 * namespace, method, and class.
 *
 * Each MLPACK namespace generally refers to one machine learning method, so
 * browsing the list of namespaces provides some insight as to the breadth of
 * the methods contained in the library.
 *
 * To generate this documentation in your own local copy of MLPACK, you can
 * simply use Doxygen, from the root directory (@c /mlpack/trunk/ ):
 *
 * @code
 * $ doxygen
 * @endcode
 *
 * @section executables Executables
 *
 * MLPACK provides several executables so that MLPACK methods can be used
 * without any need for knowledge of C++.  These executables are all
 * self-documented, and that documentation can be accessed by running the
 * executables with the '-h' or '--help' flag.
 *
 * A full list of executables is given below:
 *
 * allkfn, allknn, emst, gmm, hmm_train, hmm_loglik, hmm_viterbi, hmm_generate,
 * kernel_pca, kmeans, lars, linear_regression, local_coordinate_coding, mvu,
 * nbc, nca, pca, radical, sparse_coding
 *
 * @section tutorial Tutorials
 *
 * A few short tutorials on how to use MLPACK are given below.
 *
 *  - @ref build
 *  - @ref matrices
 *  - @ref iodoc
 *  - @ref timer
 *  - @ref sample
 *  - @ref verinfo
 *
 * Tutorials on specific methods are also available.
 *
 *  - @ref nstutorial
 *  - @ref lrtutorial
 *  - @ref rstutorial
 *  - @ref dettutorial
 *  - @ref emst_tutorial
 *  - @ref kmtutorial
 *  - @ref fmkstutorial
 *
 * @section methods Methods in MLPACK
 *
 * The following methods are included in MLPACK:
 *
 *  - Decision Stump - mlpack::decision_stump::DecisionStump
 *  - Density Estimation Trees - mlpack::det::DTree
 *  - Euclidean Minimum Spanning Trees - mlpack::emst::DualTreeBoruvka
 *  - Gaussian Mixture Models (GMMs) - mlpack::gmm::GMM
 *  - Hidden Markov Models (HMMs) - mlpack::hmm::HMM
 *  - Kernel PCA - mlpack::kpca::KernelPCA
 *  - K-Means Clustering - mlpack::kmeans::KMeans
 *  - Least-Angle Regression (LARS/LASSO) - mlpack::regression::LARS
 *  - Local Coordinate Coding - mlpack::lcc::LocalCoordinateCoding
 *  - Locality-Sensitive Hashing - mlpack::neighbor::LSHSearch
 *  - Naive Bayes Classifier - mlpack::naive_bayes::NaiveBayesClassifier
 *  - Neighborhood Components Analysis (NCA) - mlpack::nca::NCA
 *  - Nonnegative Matrix Factorization (NMF) - mlpack::amf::AMF<>
 *  - Nystroem Method - mlpack::kernel::NystroemMethod
 *  - Perceptron - mlpack::perceptron::Perceptron
 *  - Principal Components Analysis (PCA) - mlpack::pca::PCA
 *  - QUIC-SVD - mlpack::svd::QUIC_SVD
 *  - RADICAL (ICA) - mlpack::radical::Radical
 *  - Regularized SVD - mlpack::svd::RegularizedSVD
 *  - Simple Least-Squares Linear Regression -
 *        mlpack::regression::LinearRegression
 *  - Sparse Autoencoding - mlpack::nn::SparseAutoencoder
 *  - Sparse Coding - mlpack::sparse_coding::SparseCoding
 *  - Tree-based neighbor search (AllkNN, AllkFN) -
 *        mlpack::neighbor::NeighborSearch
 *  - Tree-based range search - mlpack::range::RangeSearch
 *
 * @section remarks Final Remarks
 *
 * MLPACK contributors include:
 *
 *   - Ryan Curtin <gth671b@mail.gatech.edu>
 *   - James Cline <james.cline@gatech.edu>
 *   - Neil Slagle <nslagle3@gatech.edu>
 *   - Matthew Amidon <mamidon@gatech.edu>
 *   - Vlad Grantcharov <vlad321@gatech.edu>
 *   - Ajinkya Kale <kaleajinkya@gmail.com>
 *   - Bill March <march@gatech.edu>
 *   - Dongryeol Lee <dongryel@cc.gatech.edu>
 *   - Nishant Mehta <niche@cc.gatech.edu>
 *   - Parikshit Ram <p.ram@gatech.edu>
 *   - Rajendran Mohan <rmohan88@gatech.edu>
 *   - Trironk Kiatkungwanglai <trironk@gmail.com>
 *   - Patrick Mason <patrick.s.mason@gmail.com>
 *   - Chip Mappus <cmappus@gatech.edu>
 *   - Hua Ouyang <houyang@gatech.edu>
 *   - Long Quoc Tran <tqlong@gmail.com>
 *   - Noah Kauffman <notoriousnoah@gmail.com>
 *   - Guillermo Colon <gcolon7@mail.gatech.edu>
 *   - Wei Guan <wguan@cc.gatech.edu>
 *   - Ryan Riegel <rriegel@cc.gatech.edu>
 *   - Nikolaos Vasiloglou <nvasil@ieee.org>
 *   - Garry Boyer <garryb@gmail.com>
 *   - Andreas Löf <andreas.lof@cs.waikato.ac.nz>
 *   - Marcus Edel <marcus.edel@fu-berlin.de>
 *   - Mudit Raj Gupta <mudit.raaj.gupta@gmail.com>
 *   - Sumedh Ghaisas <sumedhghaisas@gmail.com>
 *   - Michael Fox <michaelfox99@gmail.com>
 *   - Ryan Birmingham <birm@gatech.edu>
 *   - Siddharth Agrawal <siddharth.950@gmail.com>
 *   - Saheb Motiani <saheb210692@gmail.com>
 *   - Yash Vadalia <yashdv@gmail.com>
 *   - Abhishek Laddha <laddhaabhishek11@gmail.com>
 *   - Vahab Akbarzadeh <v.akbarzadeh@gmail.com>
 *   - Andrew Wells <andrewmw94@gmail.com>
 *   - Zhihao Lou <lzh1984@gmail.com>
 */

// First, include all of the prerequisites.
#include <mlpack/prereqs.hpp>

// Now the core mlpack classes.
#include <mlpack/core/util/arma_traits.hpp>
#include <mlpack/core/util/string_util.hpp>
//#include <mlpack/core/util/cli.hpp>
//#include <mlpack/core/data/load.hpp>
//#include <mlpack/core/data/save.hpp>
#include <mlpack/core/data/normalize_labels.hpp>
#include <mlpack/core/math/clamp.hpp>
#include <mlpack/core/math/random.hpp>
#include <mlpack/core/math/lin_alg.hpp>
#include <mlpack/core/math/range.hpp>
#include <mlpack/core/math/round.hpp>
//#include <mlpack/core/util/save_restore_utility.hpp>
#include <mlpack/core/dists/discrete_distribution.hpp>
#include <mlpack/core/dists/gaussian_distribution.hpp>

// Include kernel traits.
#include <mlpack/core/kernels/kernel_traits.hpp>
#include <mlpack/core/kernels/linear_kernel.hpp>
#include <mlpack/core/kernels/polynomial_kernel.hpp>
#include <mlpack/core/kernels/cosine_distance.hpp>
#include <mlpack/core/kernels/gaussian_kernel.hpp>
#include <mlpack/core/kernels/epanechnikov_kernel.hpp>
#include <mlpack/core/kernels/hyperbolic_tangent_kernel.hpp>
#include <mlpack/core/kernels/laplacian_kernel.hpp>
#include <mlpack/core/kernels/pspectrum_string_kernel.hpp>
#include <mlpack/core/kernels/spherical_kernel.hpp>
#include <mlpack/core/kernels/triangular_kernel.hpp>

#endif

// Clean up unfortunate Windows preprocessor definitions, even if this file was
// already included.  Use std::min and std::max!
#ifdef _WIN32
  #ifdef min
    #undef min
  #endif

  #ifdef max
    #undef max
  #endif
#endif
