//
//  algorithm.h
//  ccdr_proj
//
//  Created by Bryon Aragam on 3/20/14.
//  Copyright (c) 2014-2015 Bryon Aragam. All rights reserved.
//

//------------------------------------------------------------------------------/
//
// IMPORTANT NOTES FOR RCPP COMPILATION:
//  1) Functions that return Rcpp objects (e.g. List, NumericVector, etc.) should
//      NOT be const -- right now this affects the get_R function
//  2) Need to employ 'using namespace std;' unless we want to add a whole lot of
//      Rcpp::
//
//------------------------------------------------------------------------------/

#ifndef algorithm_h
#define algorithm_h

#include <vector>
#include <iostream>
#include <math.h>
#include <time.h>  // for testing and profiling only

#include "defines.h"
//#include "Auxiliary.h"
#include "SparseBlockMatrix.h"
#include "PenaltyFunction.h"
#include "CCDrAlgorithm.h"
//#include "log.h" // moved to defines.h
#include "debug.h"

//------------------------------------------------------------------------------/
//   MAIN CCDR ALGORITHM CODE
//------------------------------------------------------------------------------/

//------------------------------------------------------------------------------/
//   GLOBAL VARIABLES
//
double ZERO_THRESH = 1e-12;
// const int MAX_CCS_ARRAY_SIZE = 4000; // upper bound on the array size used in checkCycleSparse

#ifdef _DEBUG_ON_
    int ccdinit_calls = 0, ccd_calls = 0, ccs_calls = 0, spu_calls = 0, spuV_calls = 0;
#endif
//------------------------------------------------------------------------------/

// old debug code used to be here

//------------------------------------------------------------------------------/
//   MAIN CCDR CODE
//------------------------------------------------------------------------------/

// prototype for gridCCDr
std::vector<SparseBlockMatrix> gridCCDr(const std::vector<double>& cors,    // multiple arrays concatenated as a vector
                                                                            // each array contains the correlations between predictors by removing rows where each node is under intervention
                                        SparseBlockMatrix betas,            // initial guess of beta matrix
                                        std::vector<double> sigmas,
                                        const std::vector<int>& nj,         // vector containing the number of times each node is free of intervention (to replace nn)
                                        const std::vector<int>& indexj,     // index vector to indicate the start position of the array for node j in 'cors'
                                        const std::vector<double>& aj,      // weight vector for penalty term p(beta_{ij})
                                        const std::vector<double>& lambdas, // vector containing the grid of regularization parameters to be tested
                                        const std::vector<int>& weights,    // additional weights for supporting white and black listing edges
                                        const std::vector<double>& params,  // vector containing user-defined parameters: {gamma, eps, maxIters, alpha}
                                        const int verbose                   // binary variable to specify whether or not to print progress reports
                                        );

// prototype for singleCCDr
SparseBlockMatrix singleCCDr(const std::vector<double>& cors,               // multiple arrays concatenated as a vector
                                                                            // each array contains the correlations between predictors by removing rows where each node is under intervention
                             SparseBlockMatrix betas,                       // initial guess of beta matrix
                             std::vector<double> sigmas,
                             const std::vector<int>& nj,                    // vector containing the number of times each node is free of intervention (to replace nn)
                             const std::vector<int>& indexj,                // index vector to indicate the start position of the array for node j in 'cors'
                             const std::vector<double>& aj,                 // weight vector for penalty term p(beta_{ij})
                             const double lambda,                           // value of regularization parameter
                             const std::vector<int>& weights,               // additional weights for supporting white and black listing edges
                             const std::vector<double>& params,             // vector containing user-defined parameters: {gamma, eps, maxIters, alpha}
                             const int verbose                              // binary variable to specify whether or not to print progress reports
);

// prototype for computeEdgeLoss
void computeEdgeLoss(const double betaUpdate,                               // proposed new value of beta_ab
                     const unsigned int a,                                  // initial node (i.e. update beta_ab)
                     const unsigned int b,                                  // terminal node (i.e. update beta_ab)
                     const double lambda,                                   // value of regularization parameter
                     const int njb,                                         // the b-th element in 'nj'; the number of times where node b is under intervention
                     const int b1,                                          // the b-th element in 'indexj'
                     const double ajb,                                      // the b-th element in 'aj'; the weight for penalty terms p(beta_{ab}) for a=1,...,n, a !=b
                     SparseBlockMatrix& betas,                              // current value of beta matrix
                     // const PenaltyFunction& pen,                         // penalty function
                     // since the introduction of weight on penalty terms, 'aj', the penalty parameters is no longer constant
                     // we will construct PenaltyFunction only when necessary; until then we just pass in parameters such as 'aj' and 'gammaMCP'
                     const double gammaMCP,                                 // for flexible penalty parameter
                     const std::vector<double>& cors,                       // multiple arrays concatenated as a vector
                                                                            // each array contains the correlations between predictors by removing rows where each node is under intervention
                     double S[],                                            // values of the loglikelihood function in the given block
                     const int verbose                                      // binary variable to specify whether or not to print progress reports
);

// prototype for concaveCDInit
void concaveCDInit(const double lambda,                                     // value of regularization parameter
                   const std::vector<int>& weights,                         // additional weights for supporting white and black listing edges
                   const std::vector<int>& nj,                              // vector containing the number of times each node is free of intervention (to replace nn)
                   const std::vector<int>& indexj,                          // index vector to indicate the start position of the array for node j in 'cors'
                   const std::vector<double>& aj,                           // weight vector for penalty term p(beta_{ij})
                   SparseBlockMatrix& betas,                                // current value of beta matrix
                   CCDrAlgorithm& alg,                                      // CCDrAlgorithm object for this run
                   const double gammaMCP,                                   // for flexible penalty parameter
                   const std::vector<double>& cors,                         // multiple arrays concatenated as a vector
                                                                            // each array contains the correlations between predictors by removing rows where each node is under intervention
                   const int verbose                                        // binary variable to specify whether or not to print progress reports
);

// prototype for concaveCD
void concaveCD(const double lambda,                                         // value of regularization parameter
               const std::vector<int>& weights,                             // additional weights for supporting white and black listing edges
               const std::vector<int>& nj,                                  // vector containing the number of times each node is free of intervention (to replace nn)
               const std::vector<int>& indexj,                              // index vector to indicate the start position of the array for node j in 'cors'
               const std::vector<double>& aj,                               // weight vector for penalty term p(beta_{ij})
               SparseBlockMatrix& betas,                                    // current value of beta matrix
               CCDrAlgorithm& alg,                                          // CCDrAlgorithm object for this run
               const double gammaMCP,                                       // for flexible penalty parameter (this replaces PenaltyFunction& 'pen')
               const std::vector<double>& cors,                             // multiple arrays concatenated as a vector
                                                                            // each array contains the correlations between predictors by removing rows where each node is under intervention
               const int verbose                                            // binary variable to specify whether or not to print progress reports
               );

//prototype for singleUpdate
double singleUpdate(const unsigned int a,                                   // initial node (i.e. update beta_ab)
                    const unsigned int b,                                   // terminal node (i.e. update beta_ab)
                    const double lambda,                                    // value of regularization parameter
                    const int njb,                                          // the b-th element in 'nj'; the number of times where node b is under intervention
                    const int b1,                                           // the b-th element in 'indexj'
                    const double ajb,                                       // the b-th element in 'aj'; the weight for penalty terms p(beta_{ab}) for a=1,...,n, a !=b
                    const SparseBlockMatrix& betas,                         // current value of beta matrix
                    const double gammaMCP,                                  // for flexible penalty parameter (this replaces PenaltyFunction& 'pen')
                    const std::vector<double>& cors,                        // multiple arrays concatenated as a vector
                                                                            // each array contains the correlations between predictors by removing rows where each node is under intervention
                    const int verbose                                       // binary variable to specify whether or not to print progress reports
);

//prototype for singleUpdateV
double singleUpdateV(const unsigned int a,                                  // initial node (i.e. update beta_ab)
                     const unsigned int b,                                  // terminal node (i.e. update beta_ab)
                     const double lambda,                                   // value of regularization parameter
                     const int njb,                                         // the b-th element in 'nj'; the number of times where node b is under intervention
                     const std::vector<int>& indexj,                        // index vector to indicate the start position of the array for node j in 'cors'
                     const std::vector<double>& aj,                         // weight vector for penalty term p(beta_{ij})
                     SparseBlockMatrix& betas,                              // current value of beta matrix
                     const double gammaMCP,                                 // for flexible penalty parameter (this replaces PenaltyFunction& 'pen')
                     const std::vector<double>& cors,                       // multiple arrays concatenated as a vector
                                                                            // each array contains the correlations between predictors by removing rows where each node is under intervention
                     double S[],                                            // for storing the values of S1, S2
                     const int verbose                                      // binary variable to specify whether or not to print progress reports
);

//prototype for checkCycleSparse
bool checkCycleSparse(const int node,                                       // number of nodes in graph (i.e. node = pp)
                      const SparseBlockMatrix& betas,                       // sparse matrix structure
                      int a,                                                // initial node
                      int b                                                 // terminal node
);

//
// gridCCDr
//
//   Computes a full array of CCDr estimates by progressing through a grid of regularization parameters, starting
//     with the first. By default, the algorithm picks a large value lambda_max such that the first estimate is
//     guaranteed to be zero, and then the value of lambda decreases as the algorithm proceeds, allowing more and
//     more edges into the model.
//
//   Output: A vector of SparseBlockMatrix objects, one estimate for each value of lambda in lambdas
//
//   NOTES:
//     -betas and lambdas can be anything to start with
//     -the C++ code enforces no defaults; these are all implemented in R
//     -it is very important that the params values are passed in the CORRECT ORDER: {gamma, eps, maxIters, alpha}
//
std::vector<SparseBlockMatrix> gridCCDr(const std::vector<double>& cors,
                                        SparseBlockMatrix betas,
                                        std::vector<double> sigmas,
                                        const std::vector<int>& nj,
                                        const std::vector<int>& indexj,
                                        const std::vector<double>& aj,
                                        const std::vector<double>& lambdas,
                                        const std::vector<int>& weights,
                                        const std::vector<double>& params,
                                        const int verbose
                                        ){
    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG2) << "Function call: gridCCDr";
    #endif

    int nlam = static_cast<int>(lambdas.size());    // how many values of lambda are in the supplied grid?
    double alpha = params[3];                       // value of alpha; needed to know when to terminate algorithm
    std::vector<SparseBlockMatrix> grid_betas;      // the vector of SBMs that will eventually be returned

    //
    // This function is simple: Simply call singleCCDr repeatedly for each value of lambda supplied
    //
    for(int l = 0; l < nlam; ++l){
        double lambda = lambdas[l]; // current value of lambda in the grid

        //--- VERBOSE ONLY ---//
        if(verbose){
            OUTPUT << "\nWorking on lambda = " << lambda << " [" << l+1 << "/" << nlam << "]";

            #ifdef _DEBUG_ON_
                FILE_LOG(logINFO) << "Working on lambda = " << lambda << " [" << l+1 << "/" << nlam << "]";
            #endif
        }
        //--------------------//

        // To save memory, simply overwrite the same object (betas)
        // After each call to singleCCDr, we push_back the estimated object to grid_betas so there is no loss of data
        betas = singleCCDr(cors, betas, sigmas, nj, indexj, aj, lambda, weights, params, verbose);
        grid_betas.push_back(betas);

        //--- VERBOSE ONLY ---//
        if(verbose){
            OUTPUT << " | " << betas.activeSetSize() << " || " << betas.recomputeActiveSetSize(true) << std::endl;

            #ifdef _DEBUG_ON_
                OUTPUT << std::endl << "Final beta matrix: " << std::endl;
                betas.print(10);

                FILE_LOG(logINFO) << "activeSetSize = " << betas.activeSetSize() << " / recomputeActiveSetSize = " << betas.recomputeActiveSetSize(true);
            #endif
        }
        //--------------------//

        // Once an estimate has been computed, the blocks vector is superfluous and so we get rid of it to save memory
        // This is negligible for small models, but in order to push the limits of the algorithm (i.e. pp in thousands+)
        //   we need to save memory wherever possible
        grid_betas[l].clearBlocks();

        if(betas.activeSetSize() >= alpha * betas.dim()){
            break;
        }
    }

    return grid_betas;
}

//
// singleCCDr
//
//   Computes a single CCDr estimate based on an initial guess (betas) and a single value of lambda. This is _NOT_
//     what is usually meant by the "CCDr algorithm". In particular, note that the use of a full grid, starting with
//     the zero matrix and lambda_max = sqrt(n) is pivotal to returning accurate results with CCDr. This function
//     mostly defined to encapsulate the behaviour of the algorithm for a single lambda.
//
//   Output: A single SparseBlockMatrix object, representing the estimate (Phi(lambda), R(lambda))
//
//   NOTES:
//     -betas and lambda can be anything to start with
//     -the C++ code enforces no defaults; these are all implemented in R
//     -it is very important that the params values are passed in the CORRECT ORDER: {gamma, eps, maxIters, alpha}
//
SparseBlockMatrix singleCCDr(const std::vector<double>& cors,
                             SparseBlockMatrix betas,
                             std::vector<double> sigmas,
                             const std::vector<int>& nj,
                             const std::vector<int>& indexj,
                             const std::vector<double>& aj,
                             const double lambda,
                             const std::vector<int>& weights,
                             const std::vector<double>& params,
                             const int verbose
                             ){
    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG2) << "Function call: singleCCDr";
        FILE_LOG(logDEBUG1) << "Number of nonzero entries: " << betas.activeSetSize();
    #endif

    // check if sigmas will be updated
    bool updateSigmasFlag = false;
    if(sigmas[0] < 0){ // < 0 => flag for updating
        updateSigmasFlag = true;
    } else{
        // set sigmas according to user input
        for(unsigned int j = 0; j < betas.dim(); ++j){
            betas.setSigma(j, sigmas[j]);
        }
    }

    //
    // Set parameters for algorithm
    //
    if(params.size() != 4){
        OUTPUT << "Parameter vector 'params' should have exactly four elements! Check your input." << std::endl;
    }

    double gammaMCP = params[0];  // set parameter for penalty function
    double eps = params[1];
    unsigned int maxIters = params[2];
    double alpha = params[3];

    // TODO: Implement randomization
    // For now we just ignore this argument by setting it = false
    bool randomize = false;

    //
    // Create some critical objects for the algorithm
    //
    CCDrAlgorithm CCDR = CCDrAlgorithm(maxIters,
                                       eps,
                                       alpha,
                                       betas.dim(),
                                       randomize,
                                       updateSigmasFlag,
                                       LINF // use Linf norm by default (could also use L1)
    );  // to keep track of the algorithm's progress
    // PenaltyFunction MCP = PenaltyFunction(gammaMCP);                        // to compute MCP function

    //
    // Begin the main part of the algorithm
    //

    // The active set hasn't actually changed, but this guarantees that as long as the first pass of
    //  concaveCDInit doesn't add too many edges, the algorithm will do at least one sweep over the
    //  initial active set to update the edge values
    CCDR.activeSetChanged();
    do{
        //
        // Once we have run a full sweep over all active blocks, reset the stop flags to be zero
        //  and do another full sweep using concaveCDInit. If the active set changes, we keep going,
        //  otherwise, we terminate.
        //
        CCDR.resetFlags();

        // This pass runs over all blocks
        concaveCDInit(lambda, weights, nj, indexj, aj, betas, CCDR, gammaMCP, cors, verbose);

        //
        // ADD EXTRA ALGORITHM CHECKS HERE IF NEEDED
        //

        // As long as new edges have been added and we have not exceeded the maximum number of allowed edges,
        //   continue with single parameter updates for all active edges
        if(CCDR.keepGoing()){
            // block for running the rest of the CD iterations over the given active set
            int iters = 1; // we already ran one pass to determine the active set
            while( CCDR.moar(iters)){
                concaveCD(lambda, weights, nj, indexj, aj, betas, CCDR, gammaMCP, cors, verbose);
                iters++;
            }
        }

        // we have finished a full sweep
        CCDR.addSweep();

    } while( CCDR.keepGoing());

#ifdef _DEBUG_ON_
    std::ostringstream final_out;
    final_out << "\n\n";
    final_out << "#####################################################\n";
    final_out << "#    Summary                                         \n";
    final_out << "# lambda = " << lambda << std::endl;
    final_out << "# Total number of calls to concaveCDInit: " << ccdinit_calls << std::endl;
    final_out << "# Total number of calls to concaveCD: " << ccd_calls << std::endl;
    final_out << "# Total number of calls to checkCycleSparse: " << ccs_calls << std::endl;
    final_out << "# Total number of calls to singleUpdate: " << spu_calls << std::endl;
    final_out << "# Total number of calls to singleUpdateV: " << spuV_calls << std::endl;
    final_out << "#####################################################\n";
    final_out << "\n\n";

    OUTPUT << final_out.str();
    FILE_LOG(logINFO) << final_out.str();
#endif

    return betas;
}

//
// concaveCDInit
//
//   Runs a single COMPLETE sweep over all possible sigmas and blocks in the model. This sweep is intended to check
//     if any new edges need to be added to the model. First each sigma (j=1,...,pp) is updated; then each edge is
//     updated (i=1,...,pp; j > i).
//
//   Input: Note that betas & alg are all passed (and hence updated) by reference (hence void)
//   Output: void
//
//   NOTES:
//     -by default, the order of the update is to iterate across rows, starting at the top
//     -later updates should allow the user to choose the direction of the updates:
//          *across rows
//          *across columns
//          *randomly
//     -we also update sigmas before betas: what is the effect of swapping these?
//
void concaveCDInit(const double lambda,
                   const std::vector<int>& weights,
                   const std::vector<int>& nj,
                   const std::vector<int>& indexj,
                   const std::vector<double>& aj,
                   SparseBlockMatrix& betas,
                   CCDrAlgorithm& alg,
                   const double gammaMCP,
                   const std::vector<double>& cors,
                   const int verbose
                   ){

    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG1) << "Function call: concaveCDInit";
        ccdinit_calls++;
    #endif

    alg.resetError(); // sets maxAbsError = 0

    double S[2] = {0, 0};   // to store the values of the loglikelihood when comparing edges in a block; use an array instead of a vector for efficiency (faster initialization)

    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG4) << "Computing sigmas...";
    #endif

    // number of nodes
    unsigned int pp = betas.dim();

    if(alg.updateSigmas()){
        //
        // Compute sigmas
        //   See Section 4.2.2. of the computational paper for the details of this calculation
        //
        for(unsigned int j = 0; j < pp; ++j){

            double c = 0;
            for(unsigned int l = 0; l < betas.rowsizes(j); ++l){
                unsigned int row = betas.row(j, l);

                if(j <= row){
                    c += betas.value(j, l) * cors[indexj[j]*pp*(pp+1)/2 + (j + row*(row+1)/2)]; // c += beta_ij * <xj,xi>
                }
                else{
                    c += betas.value(j, l) * cors[indexj[j]*pp*(pp+1)/2 + (row + j*(j+1)/2)];   // c += beta_ij * <xj,xi> (also)
                }
            }

            double s = 0.5 * (1.0 * c + sqrt(c * c + 4 * nj[j]));
            betas.setSigma(j, s);
        } // end for loop for sigmas
    }

    #ifdef _DEBUG_ON_
        std::ostringstream sigma_out;
        for(int jj = 0; jj < pp; ++jj){
            sigma_out << betas.sigma(jj) << " ";
        }
        FILE_LOG(logDEBUG1) << "Estimated sigmas: " << sigma_out.str();

        FILE_LOG(logDEBUG4) << "Computing betas...";
    #endif

    //
    // Main loop over all edges in model (i = 0...pp-1 and j > i)
    //

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // NOTE: Looping over rows first (or even columns first if we wanted to) has strong consequences later in the algorithm.
    //        For example, edges in the first row (resp. first column) are much more likely to be nonzero than later edges.
    //        Consider how to fix this, e.g. by RANDOMIZING the order of the updates.
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for(unsigned int i = 0; i < pp; ++i){
    	for(unsigned int j = i + 1; j < pp; ++j){

    	    int weightij = weights[j * pp + i];
    	    int weightji = weights[i * pp + j];

            double betaUpdateij = 0.0;
            if(weightij >= 0) betaUpdateij = singleUpdate(i, j, weightij * lambda, nj[j], indexj[j], aj[j], betas, gammaMCP, cors, verbose);
            double betaUpdateji = 0.0;
            if(weightji >= 0) betaUpdateji = singleUpdate(j, i, weightji * lambda, nj[i], indexj[i], aj[i], betas, gammaMCP, cors, verbose);
            // consider only pass in nj[j], indexj[j], aj[j]??
            bool hasCycleij = false, hasCycleji = false;

            if(fabs(betaUpdateij) > ZERO_THRESH){
                hasCycleij = checkCycleSparse(pp, betas, i, j);
            }

            if(fabs(betaUpdateji) > ZERO_THRESH){
                // If adding i->j induces a cycle, then j->i cannot induce a cycle, so we can skip checking in this case
                if(!hasCycleij){
                    hasCycleji = checkCycleSparse(pp, betas, j, i);
                }
            }

            if(hasCycleij){
                betaUpdateij = 0.0;
            } else if(hasCycleji){
                betaUpdateji = 0.0;
            } else{
                // single parameter update for beta_ji
                computeEdgeLoss(betaUpdateji, j, i, weightji * lambda, nj[i], indexj[i], aj[i], betas, gammaMCP, cors, S, verbose);
                double S1ji = S[0]; // Qi|betaji=0
                double S2ji = S[1]; // Qi|betaji=betaUpdate

            #ifdef _DEBUG_ON_
                FILE_LOG(logDEBUG4) << "S[0] = " << S[0] << ", S[1] = " << S[1];
            #endif

                // single parameter update for beta_ij
                computeEdgeLoss(betaUpdateij, i, j, weightij * lambda, nj[j], indexj[j], aj[j], betas, gammaMCP, cors, S, verbose);
                double S1ij = S[1]; // Qj|betaij=betaUpdate
                double S2ij = S[0]; // Qj|betaij=0

            #ifdef _DEBUG_ON_
                FILE_LOG(logDEBUG4) << "S[0] = " << S[0] << ", S[1] = " << S[1];
            #endif

                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // NOTE: What if (S1ji + S1ij) == (S2ji + S2ij)??? This case should be
                //       handled carefully: When the algorithm begins we are guaranteed
                //       to encounter this case.
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if((S1ji + S1ij) <= (S2ji + S2ij)){
                    // If S1 <= S2, then beta_ij gets updated and beta_ji = 0
                    betaUpdateji = 0.0;

                    // NOTE: Instead of updating the error here, update it after the SPUs have been computed
                    // and we are ready to update the sparse structure
                } else{
                    // If S2 < S1, then beta_ji gets updated and beta_ij = 0
                    betaUpdateij = 0.0;

                    // NOTE: Instead of updating the error here, update it after the SPUs have been computed
                    // and we are ready to update the sparse structure
                }

            }

//
// DEPRECATED: This if statement is a relic from some older code. It was a bug that has been fixed, but the code
//              is left here as a reminder to handle this carefully in the future.
// [FIXED] BUG: What if the current betaij is > 0, but betaUpdateij = 0?
//
//            if(fabs(betaUpdateij) > ZERO_THRESH || fabs(betaUpdateji) > ZERO_THRESH){

            #ifdef _DEBUG_ON_
                FILE_LOG(logDEBUG4) << "Sparse update for (" << i << ", " << j << "):";
                FILE_LOG(logDEBUG4) << "beta(" << i << ", " << j << ") = " << betaUpdateij;
                FILE_LOG(logDEBUG4) << "beta(" << j << ", " << i << ") = " << betaUpdateji;
            #endif

            // Sparse update for i->j
            unsigned int row = i, col = j;

            int found = betas.find(row, col); // potential bottleneck in the code!
            std::vector<double> err(2, 0);

            if(found >= 0){
                // if the block exists in the sparse matrix, update it's value
                //
                // NOTE: This fixes the issue wherein nonzero edges could not be zeroed out

                #ifdef _DEBUG_ON_
                    // check if we are removing the edge (i,j)
                    if(fabs(betas.findValue(row, col)) > ZERO_THRESH && fabs(betaUpdateij) < ZERO_THRESH){
//                            OUTPUT << "\n\n!!!!!!!!!!!!!!!!\n";
//                            OUTPUT << "concaveCDInit: Removing edge " << "(" << i << ", " << j << ") in model!" << std::endl;
//                            OUTPUT << "!!!!!!!!!!!!!!!!\n\n";

                        FILE_LOG(logWARNING) << "concaveCDInit: Removing edge " << "(" << i << ", " << j << ") in model!";
                    }

                    // check if we are removing the edge (j,i)
                    if(fabs(betas.findValue(col, row)) > ZERO_THRESH && fabs(betaUpdateji) < ZERO_THRESH){
//                            OUTPUT << "\n\n!!!!!!!!!!!!!!!!\n";
//                            OUTPUT << "concaveCDInit: Removing edge " << "(" << j << ", " << i << ") in model!" << std::endl;
//                            OUTPUT << "!!!!!!!!!!!!!!!!\n\n";

                        FILE_LOG(logWARNING) << "concaveCDInit: Removing edge " << "(" << j << ", " << i << ") in model!";
                    }
                #endif

                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // CHECKING THIS IS A BOTTLENECK IN THE CODE: Can we speed this up somehow?
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(fabs(betas.findValue(row, col)) > ZERO_THRESH && fabs(betaUpdateij) < ZERO_THRESH){
                    alg.activeSetChanged(); // since we removed an edge to the model, the active set has changed
                }
                if(fabs(betas.findValue(col, row)) > ZERO_THRESH && fabs(betaUpdateji) < ZERO_THRESH){
                    alg.activeSetChanged(); // since we removed an edge to the model, the active set has changed
                }

                err = betas.updateBlock(col, found, betaUpdateij, betaUpdateji);

                #ifdef _DEBUG_ON_
                    if(betas.dim() <= 5){
                        FILE_LOG(logDEBUG1) << printToFile(betas, 5);
                    }
                #endif
            } else{
                // only add a block if the update is nonzero
                if(fabs(betaUpdateij) > ZERO_THRESH || fabs(betaUpdateji) > ZERO_THRESH){
                    err = betas.addBlock(row, col, betaUpdateij, betaUpdateji);
                    alg.activeSetChanged(); // since we added an edge to the model, the active set has changed

                    #ifdef _DEBUG_ON_
                        if(betas.dim() <= 5){
                            FILE_LOG(logDEBUG1) << printToFile(betas, 5);
                        }
                    #endif

                }
            }

            //
            // Update the accumulated error
            //
            alg.updateError(err[0]);
            alg.updateError(err[1]);

            #ifdef _DEBUG_ON_
                FILE_LOG(logDEBUG4) << "activeSetLength = " << betas.activeSetSize();
                FILE_LOG(logDEBUG4) << "error = " << std::setprecision(4) << alg.getError();
            #endif

            // 04/05/14: This is the only place (so far) where activeSetSize() is used
            if(betas.activeSetSize() <= alg.edgeThreshold()){
                alg.belowThreshold();
            } else{
                return; // terminate the algorithm if threshold is met
            }

        } // end for j (over columns)

    } // end for i (over rows)

    return;

}

//
// concaveCD
//
//   Runs a single PARTIAL sweep over ONLY the active set (determined by concaveCDInit). This sweep is intended to
//     update the edge weights (beta_ij) using coordinate descent. First each sigma (j=1,...,pp) is updated; then
//     each active edge is updated (j=1,...,pp; k in rows[j]).
//
//   Input: Note that betas & alg are all passed (and hence updated) by reference (hence void)
//   Output: void
//
//   NOTES:
//     -by default, the order of the update is to iterate down columns
//          ***THIS IS THE OPPOSITE OF CONCAVECDINIT
//     -would allowing random order affect the results?
//     -since we are not adding any new edges, the order of sigmas/betas should not matter here
//
void concaveCD(const double lambda,
               const std::vector<int>& weights,
               const std::vector<int>& nj,
               const std::vector<int>& indexj,
               const std::vector<double>& aj,
               SparseBlockMatrix& betas,
               CCDrAlgorithm& alg,
               const double gammaMCP,
               const std::vector<double>& cors,
               const int verbose
               ){
    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG1) << "Function call: concaveCD";
        ccd_calls++;
    #endif

    alg.resetError(); // sets maxAbsError = 0

//    double S[2] = {0, 0};   // to store the values of the loglikelihood when comparing edges in a block; use an array instead of a vector for efficiency (faster initialization)

    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG4) << "Computing sigmas...";
    #endif

    // number of nodes
    unsigned int pp = betas.dim();

    if(alg.updateSigmas()){
        //
        // Compute sigmas
        //   See Section 4.2.2. for the details of this calculation
        //
        for(unsigned int j = 0; j < pp; ++j){
            double c = 0;
            for(unsigned int l = 0; l < betas.rowsizes(j); ++l){
                unsigned int row = betas.row(j, l);

                if(j <= row){
                    c += betas.value(j, l) * cors[indexj[j]*pp*(pp+1)/2 + (j + row*(row+1)/2)]; // c += beta_ij * <xj,xi>
                }
                else{
                    c += betas.value(j, l) * cors[indexj[j]*pp*(pp+1)/2 + (row + j*(j+1)/2)];   // c += beta_ij * <xj,xi> (also)
                }
            }

            double s = 0.5 * (1.0 * c + sqrt(c * c + 4 * nj[j]));
            betas.setSigma(j, s);
        } // end for loop for sigmas
    }

    #ifdef _DEBUG_ON_
        std::ostringstream sigma_out;
        for(int jj = 0; jj < pp; ++jj){
            sigma_out << betas.sigma(jj) << " ";
        }
        FILE_LOG(logDEBUG1) << "Estimated sigmas: " << sigma_out.str();

        FILE_LOG(logDEBUG4) << "Computing betas...";
    #endif

    for(unsigned int j = 0; j < pp; ++j){
    	for(unsigned int rowIdx = 0; rowIdx < betas.rowsizes(j); ++rowIdx){
            unsigned int i = betas.row(j, rowIdx); // get the row from the sparse structure

            // By iterating over every element in the sparse structure, we are
            //  actually hitting every block twice: we can correct for this by
            //  only updating edges with j > i (upper triangle)
            //
            // Note of course that both j -> i and i -> j get updated as a block
            //  so every edge does indeed end up getting updated
            if( j <= i) continue;

            int weightij = weights[j * pp + i];
    	    int weightji = weights[i * pp + j];

            // get the current values in the block
            double betakj = betas.value(j, rowIdx);
            double betajk = betas.getSiblingValue(j, rowIdx);

            // initialize the update values
            double betaUpdateij = 0.0;
            double betaUpdateji = 0.0;

            // only update the nonzero edge
            if(fabs(betakj) > ZERO_THRESH){
                if(weightij >= 0) betaUpdateij = singleUpdate(i, j, weightij * lambda, nj[j], indexj[j], aj[j], betas, gammaMCP, cors, verbose);
            } else if(fabs(betajk) > ZERO_THRESH){
                if(weightji >= 0) betaUpdateji = singleUpdate(j, i, weightji * lambda, nj[j], indexj[i], aj[i], betas, gammaMCP, cors, verbose);
            }

            //
            // Update the edge weights no matter what below -- if a block is "zeroed-out" this is ok
            //
            std::vector<double> err = betas.updateBlock(j, rowIdx, betaUpdateij, betaUpdateji);

            #ifdef _DEBUG_ON_
                if(betas.dim() <= 5){
                    FILE_LOG(logDEBUG1) << printToFile(betas, 5);
                }
            #endif

        } // end for rowIdx

    } // end for j

    return;

}

//
// singleUpdate
//
//   Compute the value of the single parameter update (SPU) for the edge between a -> b WITHOUT calculating
//     the likelihood function. Using this function presupposes we know which edge in a block (a->b or b->a)
//     needs to be updated. If we don't know (i.e. neither edge induces a cycle), we need to compute the
//     likelihood function using singleUpdateV.
//
//   Output: The value of the single parameter update
//
//   NOTES:
//     -See Sections 4.2.1 & 4.4 for a discussion of this calculation
//
double singleUpdate(const unsigned int a,
                    const unsigned int b,
                    const double lambda,
                    const int njb,
                    const int b1,
                    const double ajb,
                    const SparseBlockMatrix& betas,
                    const double gammaMCP,
                    const std::vector<double>& cors,
                    const int verbose
                    ){

    #ifdef _DEBUG_ON_
//        FILE_LOG(logDEBUG2) << "Function call: SingleUpdate(" << a << ", " << b << ") with lambda = " << lambda;
        spu_calls++;
    #endif

    unsigned int pp = betas.dim(); // for easier access
    double betaUpdate = 0; // initialize eventual return value

    //
    // res_ab = the value of the residual factor from the paper, given by
    //    \sum_h { x_hk * r_kj^(h) } = \rho_j*<xk,xj> - \sum_{i != k} \phi_ij <xi,xk>
    //
    // Here, b = j = col, a = i = row.
    //
    double res_ab = 0;

    // Get the value: \rho_j*<xk,xj>
    if(a <= b){
        res_ab = betas.sigma(b) * cors[b1*pp*(pp+1)/2 + (a + b*(b+1)/2)];
    }
    else{
        res_ab = betas.sigma(b) * cors[b1*pp*(pp+1)/2 + (b + a*(a+1)/2)];
    }

    // Subtract the terms \phi_ij <xi,xk>
    for(unsigned int i = 0; i < betas.rowsizes(b); ++i){
        unsigned int row = betas.row(b, i);
        if(row < a){
            res_ab -= cors[b1*pp*(pp+1)/2 + (row + a*(a+1)/2)] * betas.value(b, i);
        }
        else if(row > a){ // i=a is excluded
            res_ab -= cors[b1*pp*(pp+1)/2 + (a + row*(row+1)/2)] * betas.value(b, i);
        }
    }

    //
    // The SPU is given by S_gamma(res_ab, lambda), aka evaluating the threshold function
    //   associated with the penalty function at the residual factor res_ab given the fixed
    //   values of gamma and lambda.
    //
    PenaltyFunction pen = PenaltyFunction(gammaMCP / ajb);
    betaUpdate = pen.threshold(res_ab, ajb * lambda);

    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG2) << "Function call: singleUpdate(" << a << ", " << b << ") with lambda = " << lambda << "  /  res_ab = " << res_ab;
    #endif

    return betaUpdate;
}

//
// computeEdgeLoss
//
void computeEdgeLoss(const double betaUpdate,
                     const unsigned int a,
                     const unsigned int b,
                     const double lambda,
                     const int njb,
                     const int b1,
                     const double ajb,
                     SparseBlockMatrix& betas,
                     const double gammaMCP,
                     const std::vector<double>& cors,
                     double S[],
                     const int verbose){
    //
    // Now that we have computed the value of the SPU, we need to calculate its contribution
    //   to the penalized likelihood function
    //
    // We do this by separating the calculations into two natural components, the loss (i.e. negative
    //   log-likelihood) and the penalty. Each of these values is computed twice: Once with the new
    //   value of beta_ab = betaUpdate, and once assuming beta_ab = 0.
    //
    double loss = 0, penalty = 0;
    int oldIdx_ab = betas.find(a, b); // determine whether or not the edge a->b is already in the model

    // If the edge a->b is in the model, grab its value, otherwise set it to zero
    //   This variable is used later to restore betas to its original state
    double oldBeta_ab = (oldIdx_ab >= 0) ? betas.value(b, oldIdx_ab) : 0; // 7/31/14: used to be oldIdx_ab > 0; pretty sure this was a bug since this value can be zero and still valid

    // If the edge a->b is in the model, zero it out so we can compute the effect of eliminating this edge
    //   from the model
    if(oldIdx_ab >= 0) betas.setValue(b, oldIdx_ab, 0.0); // 7/31/14: used to be oldIdx_ab > 0; pretty sure this was a bug since this value can be zero and still valid

#ifdef _DEBUG_ON_
    // Old code needed only for debugging
    int oldIdx_ba = betas.find(b, a); // THIS VALUE IS NEVER USED IN THE MAIN CODE
    double oldBeta_ba = (oldIdx_ba >= 0) ? betas.value(b, oldIdx_ba) : 0; // THIS VALUE IS NEVER USED IN THE MAIN CODE

    FILE_LOG(logDEBUG3) << oldIdx_ab << " / " << oldIdx_ba;
    FILE_LOG(logDEBUG3) << oldBeta_ab << " / " << oldBeta_ba;
#endif

    // Compute the value of the loss
    unsigned int pp = betas.dim();
    loss = betas.sigma(b) * betas.sigma(b);
    for(unsigned int m = 0; m < betas.rowsizes(b); ++m){
        unsigned int row_m = betas.row(b, m);

        for(unsigned int n = 0; n < betas.rowsizes(b); ++n){
            unsigned int row_n = betas.row(b, n);

            if(row_m <= row_n){
                loss += cors[b1*pp*(pp+1)/2 + (row_m + row_n*(row_n+1)/2)] * betas.value(b, m) * betas.value(b, n);
            } else{
                loss += cors[b1*pp*(pp+1)/2 + (row_n + row_m*(row_m+1)/2)] * betas.value(b, m) * betas.value(b, n);
            }
        }

        if(row_m <= b){
            loss -= 2.0 * betas.sigma(b) * cors[b1*pp*(pp+1)/2 + (row_m + b*(b+1)/2)] * betas.value(b, m);
        } else{
            loss -= 2.0 * betas.sigma(b) * cors[b1*pp*(pp+1)/2 + (b + row_m*(row_m+1)/2)] * betas.value(b, m);
        }
    }

    // Compute the value of the penalty
    penalty = 0;
    PenaltyFunction pen = PenaltyFunction(gammaMCP / ajb);
    for(unsigned int i = 0; i < betas.rowsizes(b); ++i){
        penalty += pen.p(fabs(betas.value(b, i)), ajb * lambda);
    }
    //penalty -= pen.p(fabs(betas[(pp * b) + b]), lambda); // ignore contribution from beta_bb (should be zero anyway!!!)

    // S[0] = value of the likelihood with beta_ab = 0
    S[0] = -1.0 * njb* log(betas.sigma(b)) + 0.5 * loss + penalty;

    //
    // S[1] = value of the likelihood with beta_ab = betaUpdate
    //   Note that S[1] = S[0] + contribution of terms involving ONLY beta_ab
    //   Since we have already computed S[0], do not bother recalculating this value and simply
    //    add up the contributions from the terms involving beta_ab.
    //
    S[1] = S[0];

    //
    // If betaUpdate = 0, then the contribution from beta_ab is zero and there is nothing to add
    //
    if(fabs(betaUpdate) > ZERO_THRESH){
        // NO NEED TO RESET THIS VALUE, JUST USE BETAUPDATE DIRECTLY (see below)
        //betas.setValue(b, oldIdx_ab, betaUpdate);
        for(unsigned int i = 0; i < betas.rowsizes(b); ++i){
            unsigned int row = betas.row(b, i);

            if(row < a){
                S[1] += 2.0 * cors[b1*pp*(pp+1)/2 + (row + a*(a+1)/2)] * betas.value(b, i) * betaUpdate;
            }
            else if(i > a){ // case when i = a is handled separately below
                S[1] += 2.0 * cors[b1*pp*(pp+1)/2 + (a + row*(row+1)/2)] * betas.value(b, i) * betaUpdate;
            }
        }
        S[1] += cors[b1*pp*(pp+1)/2 + (a + a*(a+1)/2)] * betaUpdate * betaUpdate;

        if(a <= b){
            S[1] -= 2.0 * betas.sigma(b) * cors[b1*pp*(pp+1)/2 + (a + b*(b+1)/2)] * betaUpdate;
        }
        else{
            S[1] -= 2.0 * betas.sigma(b) * cors[b1*pp*(pp+1)/2 + (b + a*(a+1)/2)] * betaUpdate;
        }

        // Of course the penalty decomposes like the loss as well, so we can just add
        //   the contribution of pen(beta_ab)
        S[1] += pen.p(fabs(betaUpdate), ajb * lambda) - pen.p(0.0, ajb * lambda); // subtract off the value at zero in case p(0) != 0
    }

    // If we modified betas, return it to it's original value
    if(oldIdx_ab >= 0) betas.setValue(b, oldIdx_ab, oldBeta_ab); // 7/31/14: used to be oldIdx_ab > 0; pretty sure this was a bug since this value can be zero and still valid

}

//
// checkCycleSparse
//
//   Determine whether or not adding an edge between a and b (i.e. a -> b) induces a cycle in the current structure
//     of the DAG represented by betas. This is a modification of the Fei's original code, modified to leverage
//     sparsity and speed up computations.
//
//   Output: 0 = no cycle induced, 1 = cycle is induced
//
//   UPDATE 05/14/14: Arrays are indeed faster than vectors, but we do not know the size of the arrays at compile time.
//                     This is a big issue; moreover, the old code did not properly initialize these arrays which is
//                     forbidden anyway when the size is not known at compile time. For now, we have reverted back to
//                     using vectors, but this should be revisited in the future.
//
//   NOTES:
//     -see Fei's paper for original source for algorithm and his original code for the original implementation
//     -the vectors that have been commented out have been left as a reminder to use ARRAYS in this function instead
//       of vectors: the allocation cost of vectors turns out to be highly nontrivial and costly here and MUST be
//       avoided
//
bool checkCycleSparse(const int node,
                      const SparseBlockMatrix& betas,
                      int a,
                      int b
                      ){

    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG3) << "Function call: checkCycleSparse(" << a << ", " << b << ")";
        ccs_calls++;
    #endif

    a++; b++; // modification to adjust for re-indexing from zero

    if(a==b)  return true;

    int i, j, k, lo, nBot = 0, nTop = 0, SLeng = 1;
    bool bCycle = false;

//    std::vector<int> color(node, 0);
//    std::vector<int> S(node, 0);
    int color[_MAX_CCS_ARRAY_SIZE_] = {};
    int S[_MAX_CCS_ARRAY_SIZE_] = {};
    color[a-1] = 1;
    S[0] = a;

    while(SLeng > 0){

        i = S[nBot] - 1;
        SLeng--;
        nBot++;

        for(k = 0; k < betas.rowsizes(i); ++k){
            j = betas.row(i, k);

            if(fabs(betas.value(i, k)) > ZERO_THRESH){
                if((j+1) == b){
                    bCycle = true;
                    break;
                } else if(color[j]==0){
                    nTop++;
                    S[nTop] = j + 1;
                    SLeng++;
                    color[j] = 1;
                }
            }
        }

        if(bCycle) break;
    }

    return bCycle;
}

#endif
