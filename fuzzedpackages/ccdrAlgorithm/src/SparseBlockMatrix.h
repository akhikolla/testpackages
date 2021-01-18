//
//  SparseBlockMatrix.h
//  ccdr_proj
//
//  Created by Bryon Aragam (local) on 3/19/14.
//  Copyright (c) 2014-2015 Bryon Aragam (local). All rights reserved.
//

#ifndef SparseBlockMatrix_h
#define SparseBlockMatrix_h

#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

#include "defines.h"
//#include "log.h"

extern double ZERO_THRESH; // defined in algorithm.h

//------------------------------------------------------------------------------/
//   SPARSE BLOCK MATRIX CLASS
//------------------------------------------------------------------------------/

//
// Stores a matrix in "Sparse-Block" format, originally intended as a memory-efficient data structure for updating
//   and accessing the entries of the adjacency matrix of a directed acyclic graph (DAG).
//
// If an adjacency matrix represents a DAG, then it can be represented by p(p-1) / 2 blocks instead of a naive
//   representation with p(p-1) entries. Since a_ij nonzero => a_ji = 0, we can store the matrix in blocks
//   {a_ij, a_ji}, where one of these values is always zero. If both are zero, it means there is no edge between
//   node i and j in the DAG. Of course, a_ii = 0 for all i.
//
//   The Sparse-Block format consists of three essential components. The first two give the standard column-major
//     implementation of a sparse matrix with a pair of (row, val) vectors for EACH column 1,...,p. The third is
//     unique to the block structure discussed above.
//
//     1) rows: This is a vector of vectors, where the jth vector stores the indices of the nonzero ROWS in column
//               j; so that rows[j][k] = i represents the edge a_ij.
//     2) vals: Also a vector of vectors, where the jth vector stores the VALUES of the nonzero entries in column
//               j; so that vals[j][k] = val represents the value of edge a_ij.
//     3) blocks: Also a vector of vectors, this makes it easy to access the "sibling" a_ji of an entry a_ij; since
//                 we are using a sparse representation, if vals[j][k] = a_ij, we cannot find a_ji by simply requesting
//                 vals[k][j] since k is the SPARSE row index, not the actual row index (i.e. k != i).
//                To overcome this difficulty, blocks[j][k] stores the SPARSE row index of a_ji given that rows[j][k] = i;
//                 i.e. blocks[j][k] = the index k' such that rows[i][k'] = j = location of a_ji.
//
//   For sufficiently sparse DAGs, the extra memory cost of the blocks vector is more than made up for by allowing O(1)
//     access to the blocks {a_ij, a_ji}, which makes the single parameter updates in CCDr very efficient.
//
//   Finally, the vector of residual variances (sigmas / rhos) in the CCDr algorithm are also stored inside this data
//     structure since it makes writing the code for CCDr easier. This vector can safely be ignored if desired, and if
//     this becomes important to users in the future, we can easily add some methods for clearing out this vector so it
//     does not take up unnecessay space.
//

//
// TRANSLATION: rows[j][k] -> location of (i,j)
//              vals[j][k] -> value of a_ij
//              blocks[j][k] -> location of (j,i)
//
// NOTE: In the ensuing code, the following convention is used
//        i, row = ACTUAL row index
//        j, col = ACTUAL column index
//        k      = SPARSE row index
//
//       Finally, note that columns are always accessed "as-is", so there should be no confusion regarding column indices.
//

//
// nonzero
//
//   Returns true if the absolute value of z is larger than the threshold et for nonzero values (see ZERO_THRESH)
//
bool nonzero(double z){
    return (fabs(z) > ZERO_THRESH);
}

class SparseBlockMatrix{

public:
    //
    // Constructors
    //
    SparseBlockMatrix(int sizeOfMatrix);                                  // Default Constructor

    SparseBlockMatrix(const std::vector< std::vector<int> >& rows_in,     //
                      const std::vector< std::vector<double> >& vals_in,  // Explicit Constructor
                      const std::vector< std::vector<int> >& blocks_in);  //

    SparseBlockMatrix(const std::vector< std::vector<int> >& rows_in,     //
                      const std::vector< std::vector<double> >& vals_in,  // Explicit Constructor
                      const std::vector< std::vector<int> >& blocks_in,   //  (Initialize sigmas too)
                      const std::vector<double>& sigmas_in);              //

    //
    // Accessor functions
    //
    int row(int j, int k) const;                            // get row index
    double value(int j, int k) const;                       // get row value
    int block(int j, int k) const;                          // get sibling row index
    double sigma(int j) const;                              // get sigma value
    int find(int row, int col) const;                       // find the sparse row in rows[col] that holds the (row, col) element
    double findValue(int row, int col) const;               // find the edge weight that correspond to the (row, col) element
    double getSiblingValue(int j, int k) const;             // user-friendly getter for accessing sibling value
    bool isEmpty(int j) const;                              // returns 1 if column j has no parents, 0 otherwise
    int rowsizes(int j) const;                              // return the size of the jth sparse row (essentially the same as neighbourhoodSizes?)
    int neighbourhoodSize(int j) const;                     // return the number of parents at node j
    int recomputeNeighbourhoodSize(int j) const;            // manually recompute the number of parents at node j
    int activeSetSize() const;                              // return the number of blocks currently in the model (activeSetLength)
    int recomputeActiveSetSize(bool reset = false);         // manually recompute the number of nonzero values in the edge set and return a warning if warn = TRUE

    //
    // Mutator functions
    //
    void setValue(int j, int k, double v);                                          // set the value of an _existing_ edge in the model
    void setSigma(int j, double s);                                                 // set the value of a residual parameter (sigma)
    std::vector<double> addBlock(int row, int col, double valij, double valji);     // add a new block (i.e. an edge) to the model with values 'valij', 'valji'
    std::vector<double> updateBlock(int row, int col, double valij, double valji);  // update the value of an _existing_ block to the model with values 'valij', 'valji'
    void clearBlocks();       // zeroes out and frees memory associated with blocks vector (which is not needed for storage and access)

    //
    // Auxiliary member functions
    //
    int dim() const;            // dimension (i.e. # of nodes) in the model
    void print() const;         // print out the _full_ beta matrix
    void print(int r) const;    // print out the upper rxr principal submatrix of betas (for suppressing large output)

#ifdef _COMPILE_FOR_RCPP_
    //
    // Constructors used by Rcpp / R
    //
    void init(Rcpp::List rows_in,
              Rcpp::List vals_in,
              Rcpp::List blocks_in,
              Rcpp::NumericVector sigmas_in);
    SparseBlockMatrix(Rcpp::List sbm);          // Explicit Constructor

    //
    // Conversion to R List
    //
    Rcpp::List get_R(double lambda_R = -1);
#endif

private:
    //
    // The main components of the data structure
    //
    std::vector< std::vector<int> > rows;       // store the sparse row indices for each column
    std::vector< std::vector<double> > vals;    // store the sparse row values for each column
    std::vector< std::vector<int> > blocks;     // store the row index for block siblings
    std::vector<double> sigmas;                 // store the residual values (sigmas) from the CCDr algorithm

    //
    // Auxiliary variables
    //
    int pp;                                     // dimension of the model
    int activeSetLength;                        // total number of nonzero edges in model (the "active set")
    std::vector<int> neighbourhoodSizes;        // store the number of parents for each node (the "neighbourhood")

    //
    // Initialization method
    //
    void init(const std::vector< std::vector<int> >& rows_in,
              const std::vector< std::vector<double> >& vals_in,
              const std::vector< std::vector<int> >& blocks_in,
              const std::vector<double>& sigmas_in);

};

//
// Initialization method
//   Separate method for initializing main values of a SparseBlockMatrix object; used since initializer lists
//   are only allowed in C++11, which R does not support on CRAN yet.
//
void SparseBlockMatrix::init(const std::vector< std::vector<int> >& rows_in,
                             const std::vector< std::vector<double> >& vals_in,
                             const std::vector< std::vector<int> >& blocks_in,
                             const std::vector<double>& sigmas_in){
    if(rows_in.size() != vals_in.size() || vals_in.size() != blocks_in.size() || blocks_in.size() != rows_in.size()){
        ERROR_OUTPUT << "Dimension mismatch in input lists: Input dimensions do not match." << std::endl;
    }

    activeSetLength = 0;                    // initialize this value zero, it will be updated as we update the data vectors
    pp = static_cast<int>(rows_in.size());  // the dimension should be equal to the number of vectors (e.g. at the first level) in any of rows / vals / blocks
    sigmas.resize(pp, 0);                   // reserve necessary memory for sigmas vector and initialize all values to zero

    if(sigmas_in.size() != pp){
        ERROR_OUTPUT << "Dimension mismatch in sigmas input: Length of sigmas must match length of rows, vals, blocks." << std::endl;
    }

    // Populate the data structure using the supplied data
    for(int j = 0; j < pp; ++j){
        rows.push_back(rows_in[j]);
        vals.push_back(vals_in[j]);
        blocks.push_back(blocks_in[j]);
        sigmas[j] = sigmas_in[j];

        // Update neighbourhood and active set sizes

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  NOTE: The way we have set up the data structure, for every nonzero edge in the model,
        //         there will be a zero edge (the 'sibling' a_ji of the nonzero edge a_ij). Thus,
        //         rows[i].size() does NOT represent the number of nonzero values in column i.
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        neighbourhoodSizes.push_back(recomputeNeighbourhoodSize(j));
        activeSetLength += neighbourhoodSizes[j];
    }
}

// Default constructor
//   Creates an empty matrix of dimension 'sizeOfMatrix' --- this constructor essentially initializes all the vectors
//   to have zero components.
//
SparseBlockMatrix::SparseBlockMatrix(int sizeOfMatrix){
    activeSetLength = 0;    // since the matrix has no nonzero edges, its active set is empty
    pp = sizeOfMatrix;      // set the dimension appropriately
    sigmas.resize(pp, 0);   // reserve necessary memory for sigmas vector and initialize all values to zero

    // Create empty vectors in each slot for rows / vals / blocks
    for(int i = 0; i < pp; ++i){
        rows.push_back(std::vector<int>(0));
        vals.push_back(std::vector<double>(0));
        blocks.push_back(std::vector<int>(0));

        neighbourhoodSizes.push_back(0);
    }
}

// Explicit constructor
//   Takes in three separate STL 'vectors of vectors' that represent (rows, vals, blocks) and sets them appropriately
//
SparseBlockMatrix::SparseBlockMatrix(const std::vector< std::vector<int> >& rows_in,
                                     const std::vector< std::vector<double> >& vals_in,
                                     const std::vector< std::vector<int> >& blocks_in){

    // If sigmas is not passed to the constructor, assume all zeroes to begin with
    std::vector<double> sigmas_in(static_cast<int>(rows_in.size()), 0);

    // Initialize the data structure
    init(rows_in, vals_in, blocks_in, sigmas_in);
}

// Explicit constructor
//   Takes in three separate STL 'vectors of vectors' that represent (rows, vals, blocks) AND a conventional
//   STL vector that represents sigmas, and sets all of them appropriately
//
SparseBlockMatrix::SparseBlockMatrix(const std::vector< std::vector<int> >& rows_in,
                                     const std::vector< std::vector<double> >& vals_in,
                                     const std::vector< std::vector<int> >& blocks_in,
                                     const std::vector<double>& sigmas_in){
    // Initialize the data structure
    init(rows_in, vals_in, blocks_in, sigmas_in);
}

//
// The getters and setters for this class are currently pretty dumb, but they are included
//  in case complications arise, and/or a decision is made to switch containers for the
//  SparseBlockMatrix class.
//
//

// j = (true) column index
// k = sparse row index
int SparseBlockMatrix::row(int j, int k) const{
    return rows[j][k];
}

double SparseBlockMatrix::value(int j, int k) const{
    return vals[j][k];
}

int SparseBlockMatrix::block(int j, int k) const{
    return blocks[j][k];
}

double SparseBlockMatrix::sigma(int j) const{
    return sigmas[j];
}

// Returns the sparse row index (i.e. in rows[col]) for the edge (row, col)
//  If the requested edge is in the model, returns the sparse row
//  Otherwise, returns -1
//      e.g. i = rows[j][find(i, j)]
//
// NOTE: If a block has been "zeroed-out" (both edges are zero), this will still find the edge
//
int SparseBlockMatrix::find(int row, int col) const{

    int found = -1; // if found < 0, then the index was not found
    for(int k = 0; k < rowsizes(col); ++k){
        if(rows[col][k] == row){
            found = k;
            break;
        }
    }

    return found;
}

// Returns the value for the edge (row, col)
//  If the requested edge is in the model, returns the value
//  Otherwise, returns 0
//      e.g. a_ij = vals[j][find(i, j)]
//
// ***Currently only used in DEBUG mode***
//
double SparseBlockMatrix::findValue(int row, int col) const{
    //
    // In DEBUG mode, return 0 if find(row, col) < 0, otherwise just return the value
    //
    #ifdef _DEBUG_ON_
        int sparse_row = find(row, col);
        if(sparse_row < 0){
            FILE_LOG(logWARNING) << "findValue called on edge which does not exist in model: " << "row = " << row << " col = " << col;
            return 0;
        } else{
            return vals[col][sparse_row];
        }
    #else
        return vals[col][find(row, col)];
    #endif
}

// For an edge represented by column j and sparse row k, get the value of the corresponding sibling
double SparseBlockMatrix::getSiblingValue(int j, int k) const{
    //
    // In DEBUG mode, check if k is not in rows[j]
    //
    #ifdef _DEBUG_ON_
        if(k > rows[j].size()){
            FILE_LOG(logERROR) << "getSiblingValue called on edge which does not exist in model:";
            FILE_LOG(logERROR) << "(Sparse row) k = " << k << " (Column) j = " << j << std::endl;
            return 0;
        }
    #endif

    return vals[rows[j][k]][blocks[j][k]];
}

// Check to see if node j has no parents in the model
bool SparseBlockMatrix::isEmpty(int j) const{
    if(rowsizes(j) > 0)
        return false;
    else
        return true;
}

// Return the number of parents of node j in the model
int SparseBlockMatrix::rowsizes(int j) const{
    return static_cast<int>(rows[j].size());
}

// Return the number of parents at node j
int SparseBlockMatrix::neighbourhoodSize(int j) const{
    return neighbourhoodSizes[j];
}

// Manually recompute the number of parents at node j
int SparseBlockMatrix::recomputeNeighbourhoodSize(int j) const{
    // Two ways to compute this:
    //  1) Compute number of zeroes, then subtract from total
    //  2) Directly compute number of nonzeroes using count_if
    //
    // Both are included below for testing purposes (which is faster???)
    //
    int numZeroes = static_cast<int>(std::count(vals[j].begin(), vals[j].end(), 0));
    int numNonZeroes = static_cast<int>(vals[j].size()) - numZeroes;

//    int numNonZeroes = static_cast<int>(std::count_if(vals[j].begin(), vals[j].end(), nonzero));

    #ifdef _DEBUG_ON_
        if(numNonZeroes < 0){
            FILE_LOG(logERROR) << "Negative result while computing neighbourhoodSize at node " << j << "!";
        }

        std::ostringstream nhbd_out;


        nhbd_out << "Neighbourhood at X" << j << ": ";
        for(int i = 0; i < vals[j].size(); ++i){
            nhbd_out << vals[j][i] << " ";
        }
        nhbd_out << " | " << numNonZeroes;
        FILE_LOG(logDEBUG4) << nhbd_out.str();
    #endif

    return numNonZeroes;
}

// Return the current active set size
int SparseBlockMatrix::activeSetSize() const{
    return activeSetLength;
}

// Recompute activeSetLength from scratch (i.e. by checking each value individually to see if > 0)
//  Should only be used for debugging purposes
//
//  If reset = true, the result of this computation will overwrite the current value of activeSetLength
//
int SparseBlockMatrix::recomputeActiveSetSize(bool reset){
    int re_activeSetLength = 0;

    for(int j = 0; j < pp; ++j){
        re_activeSetLength += recomputeNeighbourhoodSize(j);
//        for(int k = 0; k < rowsizes(j); ++k){
//            if(fabs(value(j, k)) > ZERO_THRESH)
//                re_activeSetLength++;
//        }
    }

    #ifdef _DEBUG_ON_
        if(reset && (re_activeSetLength != activeSetLength)){
            OUTPUT << "\n\n!!!!!!!!!!!!!!!!\n";
            OUTPUT << "recomputeActiveSetLength: Number of nonzero edges ( = " << re_activeSetLength << ") not equal to activeSetLength ( = " << activeSetLength << ")!!!" << std::endl;
            OUTPUT << "Updating value of activeSetLength..." << std::endl;
            OUTPUT << "!!!!!!!!!!!!!!!!\n\n";

            FILE_LOG(logWARNING) << "recomputeActiveSetLength: Number of nonzero edges ( = " << re_activeSetLength << ") not equal to activeSetLength ( = " << activeSetLength << ")!!!";
        }
    #endif

    if(reset) activeSetLength = re_activeSetLength;

    return re_activeSetLength;
}

// Update / set the value of vals[j][k]
void SparseBlockMatrix::setValue(int j, int k, double v){

    #ifdef _DEBUG_ON_
        // in debug mode, check for existence of edge first
        if(k >= rows[j].size()){
            FILE_LOG(logERROR) << "Warning: setValue called on edge that does not exist in model!" << std::endl;
            FILE_LOG(logERROR) << ">>>>>>>> Killing function." << std::endl;
            return;
        }
    #endif

    vals[j][k] = v;
}

// Update / set the jth sigma parameter
void SparseBlockMatrix::setSigma(int j, double s){
    sigmas[j] = s;
}

// Add a NEW block to the sparse-block structure
//  Returns a vector containing the difference between the old values and the updated values; since
//  we are adding a block the "old values" are always zero and this is reflected in the calculations
//
// NOTE: row and col here refer to TRUE indices, not indices in the sparse structure
std::vector<double> SparseBlockMatrix::addBlock(int row, int col, double valij, double valji){

    #ifdef _DEBUG_ON_
        if(find(row, col) >= 0){
            FILE_LOG(logERROR) << "Warning: addBlock called on edge that already exists in model! @(" << row << ", " << col << ")";
            updateBlock(col, find(row, col), valij, valji);
        } else{
            FILE_LOG(logDEBUG1) << "Adding new block at (" << row << ", " << col << "):  " << valij << " / " << valji;
        }
    #endif

    rows[col].push_back(row); // add edge (row, col)
    rows[row].push_back(col); // add edge (col, row)

    vals[col].push_back(valij);
    vals[row].push_back(valji);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Check these calculations
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    blocks[col].push_back(static_cast<int>(rows[row].size()) - 1);
    blocks[row].push_back(static_cast<int>(rows[col].size()) - 1);
    activeSetLength++;   // don't forget to update the activeSet size

    // NOTE: These values may be negative; it is up to the getError() function to implement the desired error function
    //        (e.g. L1, L2, etc)
    std::vector<double> err(2, 0);
    err[0] = valij;     // Since we are adding a new block, the old coefficient values must be zero,
    err[1] = valji;     //  so err[0] = valij - 0 = valij (similarly for err[1])

    return err;
}

// Update an EXISTING block in the sparse-block structure
//  Returns a vector containing the difference between the old values and the updated values
//
// NOTE: j and k here do not refer to true indices, but indices in the sparse structure
//       j = COLUMN index
//       k = sparse ROW index
std::vector<double> SparseBlockMatrix::updateBlock(int j, int k, double valij, double valji){

    #ifdef _DEBUG_ON_
        if(k >= rows[j].size()){
            FILE_LOG(logERROR) << "Warning: updateBlock called on edge that does not exist in model!" << std::endl;
            FILE_LOG(logERROR) << ">>>>>>>> Killing function." << std::endl;
            std::vector<double> err(2, 0);
            return err;
        }
    #endif

    // Need to store the current values before updating so we can compute the error below
    //  NOTE: This is not actually necessary, since we can compute the error BEFORE updating. If
    //        efficiency becomes a concern here we can make this change.
    double oldij = vals[j][k];
    double oldji = getSiblingValue(j, k);

    #ifdef _DEBUG_ON_
        FILE_LOG(logDEBUG2) << "Updating block at (" << rows[j][k] << ", " << j << "):  " << oldij << " / " << oldji << " --> " << valij << " / " << valji;
    #endif

    // Update the values
    setValue(j, k, valij);
    setValue(rows[j][k], blocks[j][k], valji);

    // NOTE: These values may be negative; it is up to the getError() function to implement the desired error function
    //        (e.g. L1, L2, etc)
    std::vector<double> err(2, 0);
    err[0] = valij - oldij;
    err[1] = valji - oldji;

    return err;
}

// Clear out / free the memory associated with the blocks vector
//  This is useful when passing data back to R: Once the C++ code is finished running, the blocks vector is
//  pretty much useless, and just takes up space. We free this memory before passing it back to R to keep
//  the memory footprint down while the algorithm runs.
void SparseBlockMatrix::clearBlocks(){
    #ifdef _DEBUG_ON_
        OUTPUT << "Clearing all data associated with blocks vector for this matrix.";
    #endif

    // Clear all the internal vectors associated with blocks
    for(int i = 0; i < pp; ++i) blocks[i].clear();

    // Clear the entire blocks vector itself
    blocks.clear();
}

// Returns the dimension of the model (e.g. number of nodes)
int SparseBlockMatrix::dim() const{
    return pp;
}

// print out the full betas matrix
void SparseBlockMatrix::print() const{
    for(int i = 0; i < pp; ++i){
        for(int j = 0; j < pp; ++j){
            int found = -1;
            for(int k = 0; k < rows[j].size(); ++k){
                if(rows[j][k] == i){
                    found = k;
                    break;
                }
            }

            if(found >= 0)
        #ifdef _COMPILE_FOR_RCPP_
                Rprintf("%8.2f",  vals[j][found]);
        #else
                printf("%8.2f",  vals[j][found]);
        #endif
            else
        #ifdef _COMPILE_FOR_RCPP_
                Rprintf("%8d", 0);
        #else
                printf("%8d", 0);
        #endif

        }

        OUTPUT << std::endl << std::endl;
    }
}

// print only upper rxr principal submatrix
void SparseBlockMatrix::print(int r) const{
    r = std::min(pp, r); // if r > dimension then just print the whole thing

    for(int i = 0; i < r; ++i){
        for(int j = 0; j < r; ++j){
            int found = -1;
            for(int k = 0; k < rows[j].size(); ++k){
                if(rows[j][k] == i){
                    found = k;
                    break;
                }
            }

            if(found >= 0)
        #ifdef _COMPILE_FOR_RCPP_
                Rprintf("%8.2f",  vals[j][found]);
        #else
                printf("%8.2f",  vals[j][found]);
        #endif
            else
        #ifdef _COMPILE_FOR_RCPP_
                Rprintf("%8d", 0);
        #else
                printf("%8d", 0);
        #endif
        }

        OUTPUT << std::endl << std::endl;
    }
}

#ifdef _COMPILE_FOR_RCPP_
    // Takes in three separate R lists and an explicit NumericVector (compared to e.g., below)
    void SparseBlockMatrix::init(Rcpp::List rows_in,
                                 Rcpp::List vals_in,
                                 Rcpp::List blocks_in,
                                 Rcpp::NumericVector sigmas_in){

        if(rows_in.size() != vals_in.size() || vals_in.size() != blocks_in.size() || blocks_in.size() != rows_in.size()){
            // Use Rcerr to redirect error messages to R output (still need to figure out exception handling)
            ERROR_OUTPUT << "Dimension mismatch in input lists: Input dimensions do not match." << std::endl;
        }

        activeSetLength = 0; // initialize to zero
        pp = rows_in.size();
        sigmas.resize(pp, 0); // reserve necessary memory for sigmas vector and initialize all values to zero

        if(sigmas_in.size() != pp){
            ERROR_OUTPUT << "Dimension mismatch in sigmas input: Length of sigmas must match length of rows, vals, blocks." << std::endl;
        }

        // Collect the internal vectors inside each R list to populate the data structure
        for(int j = 0; j < pp; ++j){
            // Us as<> to convert R vectors to C++ STL vectors
            rows.push_back(Rcpp::as< std::vector<int> >(rows_in[j]));
            vals.push_back(Rcpp::as< std::vector<double> >(vals_in[j]));
            blocks.push_back(Rcpp::as< std::vector<int> >(blocks_in[j]));
            sigmas[j] = sigmas_in[j];

            // Update neighbourhood and active set sizes
            neighbourhoodSizes.push_back(recomputeNeighbourhoodSize(j));
            activeSetLength += neighbourhoodSizes[j];
        }
    }

    // Takes in an R list containing the components of an R SparseBlockMatrix object:
    //   list(rows, vals, blocks, sigmas)
    SparseBlockMatrix::SparseBlockMatrix(Rcpp::List sbm){

        Rcpp::List rows_in = Rcpp::as<Rcpp::List>(sbm["rows"]);
        Rcpp::List vals_in = Rcpp::as<Rcpp::List>(sbm["vals"]);
        Rcpp::List blocks_in = Rcpp::as<Rcpp::List>(sbm["blocks"]);
        Rcpp::NumericVector sigmas_in = Rcpp::as<Rcpp::NumericVector>(sbm["sigmas"]);

        init(rows_in, vals_in, blocks_in, sigmas_in);
    }
#endif

//---------------------------------------------------------------------------------------------------//
// ***MOVED TO rcpp_wrap.h***

//
// REMOVE THIS ifdef WHEN COMPILING IN RCPP
//  This ifdef is needed to compile in Xcode, but it causes an error when compiling in Rcpp
//
//#ifdef _COMPILE_TO_RCPP_
//    // We use this namespace here to avoid overly (/ridiculously) verbose code below since just about everything sits in the Rcpp namespace
//    using namespace Rcpp;
//
//    // Returns the SparseBlockMatrix as an R list (using Rcpp); for passing data back to R
//    //  Two cases:
//    //  1) Include lambda in list (lambda_R >= 0)
//    //  2) Ignore lambda (lambda_R < 0)
//    List SparseBlockMatrix::get_R(double lambda_R){
//        if(lambda_R < 0)
//            return List::create(_["rows"] = wrap(rows), _["vals"] = wrap(vals), _["sigmas"] = wrap(sigmas), _["blocks"] = wrap(blocks), _["length"] = wrap(activeSetLength));
//        else
//            return List::create(_["rows"] = wrap(rows), _["vals"] = wrap(vals), _["sigmas"] = wrap(sigmas), _["blocks"] = wrap(blocks), _["length"] = wrap(activeSetLength), _["lambda"] = wrap(lambda_R));
//    }
//#endif
//---------------------------------------------------------------------------------------------------//

#endif
