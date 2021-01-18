#include <R.h>
//#include <Rinternals.h>
//#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL

//#include <Rcpp.h> //fixing: In file included from calc_anclikes_sp.cpp:2:
// /Library/Frameworks/R.framework/Versions/2.15/Resources/library/RcppArmadillo/include/RcppArmadillo.h:26:6: 
// error: #error "The file 'Rcpp.h' should not be included. Please correct to include only 'RcppArmadillo.h'."


// min_element/max_element
#include <iostream>
#include <string>	// for e.g. std::string
#include <algorithm>
#include <cmath> 	// for e.g. abs (absolute value: http://www.cplusplus.com/reference/cmath/abs/ ) 
					// apparently, when you reference <cmath>, the functions can be accessed via e.g.
					// std::abs.  Source: http://stackoverflow.com/questions/8734230/math-interface-vs-cmath-in-c
#include "basics.h"

// Try this to expose FORTRAN functions (you have to make sure .f files are included in the compilation)
//#include "expokit_wrappers.cpp"		// for use of EXPOKIT exponentiation code NOT FASTER UGH!
using namespace std;


// Try this to expose FORTRAN functions
/*extern"C" {
	void wrapdgpadm_(int * ideg,int * m,double * t,double * H,int * ldh, double * wsp,int * lwsp,int * ipiv,int * iexph,int *ns,int *iflag );
	}
*/



extern "C" 
{
SEXP cpp_combn_zerostart(SEXP R_n, SEXP R_m, SEXP R_maxval)
	{
	using namespace std;
	
	// Convert to plain C++
	int cpp_nval = Rcpp::as<int>(R_n);
	int cpp_mval = Rcpp::as<int>(R_m);
	int cpp_maxval = Rcpp::as<int>(R_maxval);
	
	// Error checks
	// both values must be greater than 0, and n must be >= m
	// (you can't have 5 choose 10)
	if (cpp_nval < cpp_mval)
		{
		return Rcpp::wrap(0);
		}
	if (cpp_nval < 1)
		{
		return Rcpp::wrap(0);
		}
	if (cpp_nval < 1)
		{
		return Rcpp::wrap(0);
		}
	
	
	// Create pointer variables to hold the addresses to each
	int* n = &cpp_nval;
	int* m = &cpp_mval;
	
	// Choose n by m; calculate from the values stored at addresses n and m
	int Cnm;
	Cnm = nChoosek(*n, *m);
	
	// Error check
	if ((Cnm > cpp_maxval) || (Cnm < 1))
		{
		//cout << "\nERROR: n=" << cpp_nval << ", k=" << cpp_mval << ", n choose k=" << Cnm << " > maxval=", cpp_maxval;
		//cout << "\nCalculating something this big may crash your computer!  Returning 0.";
		return Rcpp::wrap(0);
		}
	
	
	// Declare and populate empty array of addresses to hold the combn results
	// Addresses for a 10x3 array; 30 total
	int* combmat = new int[(int)Cnm**(m+0)];
	
	// Run moncombn_zerostart; this will update the stuff in the addresses
	// stored in combmat
	moncombn_zerostart(combmat,n,m);
	
	// Write the contents of each of the Cnm times *m addresses to cout
	// convert double (float) Cnm to int
	int nrows;
	int ncols;
	nrows = *m;
	ncols = (int)Cnm;
	
	//vecsize now longer used, remove
	//int vecsize;
	//vecsize = nrows * ncols;
	//int combmat_vals[nrows * ncols];
	
	//Rcpp::NumericVector combmat_vals(vecsize);	// vector of size vecsize filled with 0s
	Rcpp::IntegerMatrix combmat_vals(nrows,ncols);	// vector of size vecsize filled with 0s

	// initialize row & column numbers, and the temporary number	
	int rownum = 0;
	int colnum = 0;
	int tmpnum = 0;
	
	//cout << nrows << " rows, " << ncols << "cols...\n";
	
	for (int j = 1; j <= Cnm**(m+0); j++)
		{
		//cout << "\n";
		//cout << *(combmat+j-1) << " ";
		//combmat_vals[rownum][colnum] = *(combmat+j-1);
		//tmpnum = Rcpp::as<int>(*(combmat+j-1));
		tmpnum = *(combmat+j-1);
		//combmat_vals[j-1] = tmpnum;

		//cout << "\n" << rownum << "," << colnum << ", ncols=" << ncols << ": " << tmpnum;
		combmat_vals(rownum, colnum) = tmpnum;
		
		// Increment column
		rownum++;
		
		// Reset column when you reach the end; increment the row.
		// Reset of rows is not necessary
		if (rownum >= (nrows))
			{
			rownum = 0;
			colnum++;
			}
		} // end forloop

	// Convert to an int
	//Rcpp::Matrix outcombs(combmat_vals);
	
	// Delete the object which combmat points to, which was opened by new
	delete[] combmat;
	// Set the pointer itself to NULL
	combmat = NULL;
	
	// example use Armadillo matrix
	// http://dirk.eddelbuettel.com/blog/2011/04/23/
	//arma::mat outcombs = Rcpp::as<arma::mat>(combmat_vals);
	
	//return Rcpp::wrap(combmat_vals);
	return combmat_vals;
	}
} // end extern "C"

	
// areas_list here, is just a vector of indices (0, 1, 2, etc...)
// areas_list here, is just a vector of indices (0, 1, 2, etc...)
extern "C" 
{

SEXP cpp_areas_list_to_states_list(SEXP R_areas_indices, SEXP R_maxareas, SEXP R_include_null_range)
	{
	using namespace std;

	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::IntegerVector areas_indices(R_areas_indices);
	int maxareas = Rcpp::as<int>(R_maxareas);
	bool include_null_range = Rcpp::as<bool>(R_include_null_range);
	
	// calculate the number of states, based on the number of areas
	int total_numareas = (int)areas_indices.size();

	// Declare numstates and initialize count with 1 or 0, depending on if null range is included
	// (default should be TRUE)
	int numstates;
	if (include_null_range == true)
		{
		numstates = 1;
		}
	else
		{
		numstates = 0;		
		}
	// This adds up the areas
	for (int i=1; i<=maxareas; i++)
		{
		int numstates_for_this_numareas = nChoosek(total_numareas, i);
		numstates = numstates + numstates_for_this_numareas;
		}
	
	// Create a list; make the null range NA if needed
	Rcpp::List states_list(numstates);
	int states_list_pos = 0;
	
	
	// I guess 1 == TRUE, here
	//cout << total_numareas << " " << numstates << " " << include_null_range << TRUE << FALSE << endl;
	if (include_null_range == 1)
		{
		states_list_pos = 0;
		
		// insert NA
		// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2011-February/001837.html
		states_list[states_list_pos] = -1;	// NA_REAL is for float(?)
		states_list_pos = 1;
		}
	else
		{
		// Just start listing states at 0!
		states_list_pos = 0;
		}


	// 
	int tmp_R_maxval = 1e+07;
	SEXP R_maxval = Rcpp::wrap(tmp_R_maxval);
	
	// Go through the states for this range size
	for (int tmp_numareas=1; tmp_numareas<=maxareas; tmp_numareas++)
		{
		// declare tmp_matrix_combinations as an SEXP
		//SEXP tmp_matrix_combinations;
		SEXP R_n = Rcpp::wrap(total_numareas);
		SEXP R_m = Rcpp::wrap(tmp_numareas);
		
		// Calculate the combn matrix of areas for this rangesize
		// Convert to NumericMatrix
		//tmp_matrix_combinations = cpp_combn_zerostart(R_n, R_m, R_maxval);
		//Rcpp::IntegerMatrix combmat_vals(tmp_matrix_combinations);
		Rcpp::IntegerMatrix combmat_vals = cpp_combn_zerostart(R_n, R_m, R_maxval); // screwing with this leads to AddressSanitizer faults
		
		// Put the states into the states_list
		
		// get the number of combinations, i.e. number of columns
		// Choose n by m; calculate from the values stored at addresses n and m
		int Cnm;
		Cnm = nChoosek(total_numareas, tmp_numareas);
		int ncols = Cnm;
		int nrows = tmp_numareas;
		
		for (int i=0; i<ncols; i++)
			{
			// Go through the rows of the combn matrix
			vector<int> areas_in_this_state;
			for (int j=0; j<nrows; j++)
				{
				// Get the subitem
				int subitem;
				subitem = combmat_vals(j,i);
				areas_in_this_state.push_back(subitem);
				}
			
			states_list[states_list_pos] = areas_in_this_state;
			
			// Increment states_list position
			states_list_pos++;
			}
		} //end go through the states for this range size
	
	// states_list is an Rcpp::List object, we should be able to
	// just use it...
	return states_list;
	}
} // end extern "C"








/* Rcpp function: take a geographic states list, dmat, and elist, and produce a (dense) dispersal */
/* probability matrix. */
/* (Later step: normalize it.) */
extern "C" 
{
SEXP cpp_states_list_to_DEmat(SEXP R_areas_indices, SEXP R_states_indices, SEXP R_dmat, SEXP R_elist, SEXP R_amat, SEXP
R_normalize_TF) {

	using namespace std;
	float zero_14_digit_precision = 0.000000000000001;

	// Convert to plain C++
	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::List areas_list(R_areas_indices);
	Rcpp::List states_indices(R_states_indices);
	Rcpp::NumericMatrix dmat(R_dmat);
	Rcpp::NumericVector elist(R_elist);
	Rcpp::NumericMatrix amat(R_amat);
	
	// Instead of TF, use 1/0
	int normalize_TF = Rcpp::as<int>(R_normalize_TF);
	
	/*Rcpp::NumericVector xb(rightprobs);
	float xs = Rcpp::as<float>(s);
	float xv = Rcpp::as<float>(v);
	float xj = Rcpp::as<float>(j);
	float xy = Rcpp::as<float>(y);
	
	Rcpp::NumericVector dmatc(dmat);
	int max_min_rangesize_c = Rcpp::as<int>(max_min_rangesize);
	Rcpp::NumericMatrix maxent01sub(maxent01s);
	Rcpp::NumericMatrix maxent01vic(maxent01v);
	Rcpp::NumericMatrix maxent01jump(maxent01j);
	Rcpp::NumericMatrix maxent01symp(maxent01y);
	*/
	/* Get the sizes of the vectors */
	/*int n_xa = xa.size();
	int n_xb = xb.size();
	int numstates = xl.size();
	*/
	
	/* Store the states_indices as a vector of vectors */
	/* http://stackoverflow.com/questions/8159872/r-list-of-numeric-vectors-c-2d-array-with-rcpp */

	int numareas = areas_list.size();
	int numstates = states_indices.size();
	
	/* Empty vector for the states */
	vector< vector<int> > states_vecs;
	for(int i=0; i<numstates; i++)
		{
		/* Get the subvector */
		vector<int> tmp = Rcpp::as<vector<int> > (states_indices[i]);
		
		/* Push to the original list */
		states_vecs.push_back(tmp);
		}
	
	// Make a DEmat
	int nrows = numstates;
	int ncols = numstates;
	
	// COO vs. regular
	
	// NOTE: DON'T DECLARE THINGS INSIDE if/for loops...because:
	// http://www.cplusplus.com/forum/beginner/27245/
	// If you declare something inside a while, for, if, function, it survives only in that scope.

	// square matrix
	Rcpp::NumericMatrix DEmat(nrows,ncols);
	
	// Create DEmat
	// loop through starting states
	for(int i=0; i<numstates; i++)
		{
		// starting state is a vector of area indices
		vector<int> starting_state = states_vecs[i];


		// loop through ending states
		for(int j=0; j<numstates; j++)
			{
			// ending state is a vector of area indices
			vector<int> ending_state = states_vecs[j];
			
			// Starting value for the cell
			float tmpval = 0.0;

			// Check that there is 1 more ending state than starting
			// state, and that all of the starting states are inside the ending state
			//bool boolTF = false; // not used

			// Are the starting and ending range sizes both size 1?
			// Then, amat might be relevant
			// ALSO -- exclude start=NULL, end=NULL
			if (( (starting_state.size() == 1) && (ending_state.size() == 1) ) && (starting_state[0] >= 0) && (ending_state[0] >= 0) )
				{

				//std::cout << starting_state.size() << "\n";
				//std::cout << ending_state.size() << "\n";
				
				// We can do this since they are of size 1
				int first_pos = 0;
				//std::cout << first_pos << "\n";
				
				int ancestral_area_index = starting_state[first_pos];
				int descendent_area_index = ending_state[first_pos];
				
				//std::cout << ancestral_area_index << "\n";
				//std::cout << descendent_area_index << "\n";
				
				// Skip, of course, if you are on the diagonal
				// Otherwise do it
				if (ancestral_area_index != descendent_area_index)
					{
					float amat_prob = amat(ancestral_area_index, descendent_area_index);
					//std::cout << amat_prob << "\n";
					DEmat(i,j) = amat_prob;
					//std::cout << DEmat(i,j) << "\n";
					//continue;
					}
				}
				


			// Is the ending range size 1 bigger than the starting range size?
			if ( (starting_state.size()+1) == ending_state.size() )
				{
				// are all starting states found in the ending states?
				if (all_ints_found(starting_state, ending_state) == true)
					{
					// print some stuff
					/* 
					int tmp = get_missing_int (ending_state, starting_state);
					//std::cout << starting_state.size() << "\n";
					//std::cout << ending_state.size() << "\n";
					printvec(starting_state);
					printvec(ending_state);
					//std::cout << "missing state=" << tmp << "\n";
					//std::cout << "\n";
					//std::cout << "\n";
					*/

					
					// Find the missing range (-1 hopefully throws an error as an index)
					int missing_end_area = -1;
					missing_end_area = get_missing_int(ending_state, starting_state);
					// The probability of reaching this state is the sum of ways to get there
					// from dmat
					
					// Go through the starting areas, calculate probability
					for (unsigned int k=0; k<starting_state.size(); k++)
						{
						int ancestral_area_index = starting_state[k];
						int descendent_area_index = missing_end_area;
						float dmat_prob = dmat(ancestral_area_index, descendent_area_index);
						
						// Add to the cell's tmpval
						tmpval = tmpval + dmat_prob;
						}
					// go to next iteration of the loop
					//std::cout << tmpval << " " ;
					DEmat(i,j) = tmpval;
					continue;
					//std::cout << "\nmissing_end_area=" << missing_end_area << ": " << tmpval << "\n";
					}
				}

			// Is the ending range NULL (i.e., -1) and the starting range is size 1 and is NOT -1
			//cout << ending_state[0] << " " << starting_state[0] << "\n";
			if ( (ending_state[0] < 0) && (starting_state[0] >= 0) )
				{
				if ( (starting_state.size() == 1) )
					{
					// Find the missing range (-1 hopefully throws an error as an index)
					int missing_start_area = -1;
					missing_start_area = starting_state[0];
					// The probability of reaching this state is the sum of ways to get there
					// from dmat

					// Get the extinction probability, and add it
					int descendent_area_index = missing_start_area;
					float e_prob = elist[descendent_area_index];

					tmpval = tmpval + e_prob;
					//std::cout << tmpval << " " ;
					DEmat(i,j) = tmpval;
					continue;
					}
				}


			// Is the ending range size 1 *smaller* than the starting range size?
			if ( starting_state.size() == (ending_state.size()+1) )
				{
				// are all ending states found in the starting states?
				if (all_ints_found(ending_state, starting_state) == true)
					{
					// print some stuff
					/* 
					int tmp = get_missing_int (ending_state, starting_state);
					//std::cout << starting_state.size() << "\n";
					//std::cout << ending_state.size() << "\n";
					printvec(starting_state);
					printvec(ending_state);
					//std::cout << "missing state=" << tmp << "\n";
					//std::cout << "\n";
					//std::cout << "\n";
					*/

					
					// Find the missing range (-1 hopefully throws an error as an index)
					int missing_start_area = -1;
					missing_start_area = get_missing_int(starting_state, ending_state);
					// The probability of reaching this state is the sum of ways to get there
					// from dmat

					// Get the extinction probability, and add it
					int descendent_area_index = missing_start_area;
					float e_prob = elist[descendent_area_index];
					
					// Add to the cell's tmpval
					tmpval = tmpval + e_prob;
					
					// Go to next iteration of the loop
					//std::cout << tmpval << " " ;
					DEmat(i,j) = tmpval;
					continue;
					}
				}
			//std::cout << tmpval << " " ;
			} // end loop through starting/ending states
		//cout << "\n" ;
		}
	
	// Normalize DEmat by row
	if (normalize_TF == 1)
		{
		// Go through the matrix by row
		for(int i=0; i<numstates; i++)
			{
			// sum the non-diagonal values in the row
			float rowsum=0.0;
			
			for(int j=0; j<numstates; j++)
				{
				// Skip, if it's on the diagonal
				if (i == j)
					{
					//Causes a warning in R CMD check:
					//"calc_anclikes_sp.cpp:505:13: warning: explicitly assigning a variable of type 'float' to itself [-Wself-assign]
                    //                  rowsum = rowsum;
                    //                    ~~~~~~ ^ ~~~~~~
					//	1 warning generated.
					//
					// So, commenting this out:
					//rowsum = rowsum;
					float junkval = 1.0;
					junkval = junkval + 1;
					}
				else
					{
					rowsum = rowsum + DEmat(i,j);
					}
				}

			// check precision issue
			// don't divide by rowsum unless you are normalizing diag to 1!!
			/*
			if (rowsum >= zero_14_digit_precision)
				{
				// Now, divide everything in this row by the rowsum
				for(int j=0; j<numstates; j++)
					{
					DEmat(i,j) = DEmat(i,j) / rowsum;
					}
				}
			*/
			
			// Put the negative rowsum on the diagonal
			// Make the diagonal 0.0 if the rowsum is super-duper low...
			int j = i;
			if (rowsum < zero_14_digit_precision)
				{
				DEmat(i,j) = 0.0;
				} else {
				DEmat(i,j) = -1 * rowsum;			
				}
			} // end processing this row
		} // end if normalize_TF


	// Print DEmat
	for(int i=0; i<numstates; i++)
		{
		for(int j=0; j<numstates; j++)
			{
			//cout << DEmat(i,j);
			}
		//cout << "\n" ;
		}

	// Print dmat
	for(int i=0; i<numareas; i++)
		{
		for(int j=0; j<numareas; j++)
			{
			//cout << dmat(i,j) << " ";
			}
		//cout << "\n" ;
		}

	//cout << "\n" ;
	//cout << DEmat.nrow() ;
	//cout << DEmat.ncol() ;
	
	return DEmat;
	}
} // end extern "C"








/* Rcpp function: take a geographic states list, dmat, and elist, and produce a (COO) dispersal */
/* probability matrix. */
/* (Later step: normalize it.) */
extern "C" 
{
SEXP cpp_states_list_to_DEmat_COO(SEXP R_areas_indices, SEXP R_states_indices, SEXP R_dmat, SEXP R_elist, SEXP R_amat, 
SEXP R_normalize_TF, SEXP R_min_precision) {

	using namespace std;
	
	// Hard-code default; values below this are treated as 0.0 (not included in sparse matrix)
	//float zero_14_digit_precision = 0.000000000000001;
	// Soft-code default
	float zero_14_digit_precision = Rcpp::as<float>(R_min_precision);
	
	// Convert to plain C++
	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::List areas_list(R_areas_indices);
	Rcpp::List states_indices(R_states_indices);
	Rcpp::NumericMatrix dmat(R_dmat);
	Rcpp::NumericVector elist(R_elist);
	Rcpp::NumericMatrix amat(R_amat);
	
	// Instead of TF, use 1/0
	int normalize_TF = Rcpp::as<int>(R_normalize_TF);
	
	
	/*Rcpp::NumericVector xb(rightprobs);
	float xs = Rcpp::as<float>(s);
	float xv = Rcpp::as<float>(v);
	float xj = Rcpp::as<float>(j);
	float xy = Rcpp::as<float>(y);
	
	Rcpp::NumericVector dmatc(dmat);
	int max_min_rangesize_c = Rcpp::as<int>(max_min_rangesize);
	Rcpp::NumericMatrix maxent01sub(maxent01s);
	Rcpp::NumericMatrix maxent01vic(maxent01v);
	Rcpp::NumericMatrix maxent01jump(maxent01j);
	Rcpp::NumericMatrix maxent01symp(maxent01y);
	*/
	/* Get the sizes of the vectors */
	/*int n_xa = xa.size();
	int n_xb = xb.size();
	int numstates = xl.size();
	*/
	
	/* Store the states_indices as a vector of vectors */
	/* http://stackoverflow.com/questions/8159872/r-list-of-numeric-vectors-c-2d-array-with-rcpp */

	int numareas = areas_list.size();
	int numstates = states_indices.size();
	
	/* Empty vector for the states */
	vector< vector<int> > states_vecs;
	for(int i=0; i<numstates; i++)
		{
		/* Get the subvector */
		vector<int> tmp = Rcpp::as<vector<int> > (states_indices[i]);
		
		/* Push to the original list */
		states_vecs.push_back(tmp);
		}
	
	// Make a DEmat
	int nrows = numstates;
	int ncols = numstates;
	
	// COO vs. regular
	
	// NOTE: DON'T DECLARE THINGS INSIDE if/for loops...because:
	// http://www.cplusplus.com/forum/beginner/27245/
	// If you declare something inside a while, for, if, function, it survives only in that scope.

	// square matrix
	Rcpp::NumericMatrix DEmat(nrows,ncols);


	// COO-formatted setup
	
	// three columns, pre-transposed
	int num_nonzeros;
	num_nonzeros = 0;	// number of nonzeros
	// initial allocation of the array size
	//int numcells = 100;	// not used
	
	/*
	// Declare pointer variables; 0 means no address allocated (null pointer)
	// nothrow argument: http://www.cplusplus.com/reference/std/new/nothrow/
	int * ia;				// row position (after tranposition, which we are forcing)
	int * ja;				// col position (after tranposition, which we are forcing)
	float * a;					// cell values
	
	ia = new (nothrow) int[numcells];
	if (ia == 0) {
		// error assigning memory. Take measures.
		//cout << "\nError: memory for ia (" << numcells << " cells) could not be allocated.\n";
		};
	ja = new (nothrow) int[numcells];
	if (ja == 0) {
		// error assigning memory. Take measures.
		//cout << "\nError: memory for ja (" << numcells << " cells) could not be allocated.\n";
		};
	//a = new (nothrow) float [numcells];
	a = new (nothrow) float[numcells];
	if (a == 0) {
		// error assigning memory. Take measures.
		//cout << "\nError: memory for a (" << numcells << " cells) could not be allocated.\n";
		};
	//cout << ia << " " << ja << " " << a << " " << endl;
	*/
	
	vector<int> ia;
	vector<int> ja;
	vector<float> a;
	
	
	// Create DEmat
	
	// loop through starting states
	for(int i=0; i<numstates; i++)
		{
		// starting state is a vector of area indices
		vector<int> starting_state = states_vecs[i];


		// loop through ending states
		for(int j=0; j<numstates; j++)
			{
			// ending state is a vector of area indices
			vector<int> ending_state = states_vecs[j];
			
			// Starting value for the cell
			float tmpval = 0.0;

			// Check that there is 1 more ending state than starting
			// state, and that all of the starting states are inside the ending state
			//bool boolTF = false;	// not used


			// Are the starting and ending range sizes both size 1?
			// Then, amat might be relevant
			// EXCLUDE NULL ranges
			if (( (starting_state.size() == 1) && (ending_state.size() == 1) ) && (starting_state[0] >= 0) && (ending_state[0] >= 0) )
				{
				// We can do this since they are of size 1
				int ancestral_area_index = starting_state[0];
				int descendent_area_index = ending_state[0];
				
				// Skip, of course, if you are on the diagonal
				// Otherwise do it
				if (ancestral_area_index != descendent_area_index)
					{
					tmpval = amat(ancestral_area_index, descendent_area_index);
					float amat_prob = amat(ancestral_area_index, descendent_area_index);
					tmpval = amat_prob;

					// If tmpval is >0, push to output list
					if (tmpval >= zero_14_digit_precision)
						{
						num_nonzeros++;
						ia.push_back(j+1);
						ja.push_back(i+1);
						a.push_back(tmpval);
						}
					continue;
					}
				}

			// Is the ending range size 1 bigger than the starting range size?
			if ( (starting_state.size()+1) == ending_state.size() )
				{
				// are all starting states found in the ending states?
				if (all_ints_found(starting_state, ending_state) == true)
					{
					// print some stuff
					/* 
					int tmp = get_missing_int (ending_state, starting_state);
					//std::cout << starting_state.size() << "\n";
					//std::cout << ending_state.size() << "\n";
					printvec(starting_state);
					printvec(ending_state);
					//std::cout << "missing state=" << tmp << "\n";
					//std::cout << "\n";
					//std::cout << "\n";
					*/

					
					// Find the missing range (-1 hopefully throws an error as an index)
					int missing_end_area = -1;
					missing_end_area = get_missing_int(ending_state, starting_state);
					// The probability of reaching this state is the sum of ways to get there
					// from dmat
					
					// Go through the starting areas, calculate probability
					for (unsigned int k=0; k<starting_state.size(); k++)
						{
						int ancestral_area_index = starting_state[k];
						int descendent_area_index = missing_end_area;
						float dmat_prob = dmat(ancestral_area_index, descendent_area_index);
						
						// Add to the cell's tmpval
						tmpval = tmpval + dmat_prob;
						}
					// go to next iteration of the loop
					//std::cout << tmpval << " " ;
					
					// If tmpval is >0, push to output list
					if (tmpval >= zero_14_digit_precision)
						{
						num_nonzeros++;
						ia.push_back(j+1);
						ja.push_back(i+1);
						a.push_back(tmpval);
						}
					continue;
					//std::cout << "\nmissing_end_area=" << missing_end_area << ": " << tmpval << "\n";
					}
				}

			// Is the ending range NULL (i.e., -1) and the starting range is size 1 and is NOT -1
			//cout << ending_state[0] << " " << starting_state[0] << "\n";
			if ( (ending_state[0] < 0) && (starting_state[0] >= 0) )
				{
				if ( (starting_state.size() == 1) )
					{
					// Find the missing range (-1 hopefully throws an error as an index)
					int missing_start_area = -1;
					missing_start_area = starting_state[0];
					// The probability of reaching this state is the sum of ways to get there
					// from dmat

					// Get the extinction probability, and add it
					int descendent_area_index = missing_start_area;
					float e_prob = elist[descendent_area_index];

					tmpval = tmpval + e_prob;
					//std::cout << tmpval << " " ;

					// If tmpval is >0, push to output list
					if (tmpval >= zero_14_digit_precision)
						{
						num_nonzeros++;
						ia.push_back(j+1);
						ja.push_back(i+1);
						a.push_back(tmpval);
						}
					continue;
					}
				}


			// Is the ending range size 1 *smaller* than the starting range size?
			if ( starting_state.size() == (ending_state.size()+1) )
				{
				// are all ending states found in the starting states?
				if (all_ints_found(ending_state, starting_state) == true)
					{
					// print some stuff
					/* 
					int tmp = get_missing_int (ending_state, starting_state);
					//std::cout << starting_state.size() << "\n";
					//std::cout << ending_state.size() << "\n";
					printvec(starting_state);
					printvec(ending_state);
					//std::cout << "missing state=" << tmp << "\n";
					//std::cout << "\n";
					//std::cout << "\n";
					*/

					
					// Find the missing range (-1 hopefully throws an error as an index)
					int missing_start_area = -1;
					missing_start_area = get_missing_int(starting_state, ending_state);
					// The probability of reaching this state is the sum of ways to get there
					// from dmat

					// Get the extinction probability, and add it
					int descendent_area_index = missing_start_area;
					float e_prob = elist[descendent_area_index];
					
					// Add to the cell's tmpval
					tmpval = tmpval + e_prob;
					
					// Go to next iteration of the loop
					//std::cout << tmpval << " " ;

					// If tmpval is >0, push to output list
					if (tmpval >= zero_14_digit_precision)
						{
						num_nonzeros++;
						ia.push_back(j+1);
						ja.push_back(i+1);
						a.push_back(tmpval);
						}
					continue;
					}
				}
			//std::cout << tmpval << " " ;
			} // end loop through starting/ending states
		//cout << "\n" ;
		}
	
	// Subtract 1 to take out the last ++;
	num_nonzeros = num_nonzeros-1;

	// Normalize DEmat by row
	if (normalize_TF == 1)
		{
		// Go through the COO matrix by row

		for(int i=0; i<numstates; i++)
			{
			// sum the non-diagonal values in the row (here, rows are ja)
			float rowsum=0.0;
			
			// go through all the indices of the COO matrix, add items that match (except diagonal)
			for (int tmp_COO_index=0; tmp_COO_index<(num_nonzeros); tmp_COO_index++)
				{
				// if ja (row) equals the row number (i), but isn't the diagonal, then add
				//if ( (ja[tmp_COO_index] == i) && ( ia[tmp_COO_index] != i) )
				if ( ( (ja[tmp_COO_index]-1) == i) && (  (ia[tmp_COO_index]-1) != i) )
					{
					rowsum = rowsum + a[tmp_COO_index];
					}
				}
			
			// Add the negative rowsum to the COO matrix at the diagonal position
			// but if rowsum is VERY small, DON'T add...
			if (rowsum >= zero_14_digit_precision)
				{
				num_nonzeros++;
				ia.push_back(i+1);
				ja.push_back(i+1);
				a.push_back(-1.0 * rowsum);
				} // end if rowsum < zero_14_digit_precision check
			} // end processing this row
		} // end if normalize_TF

	// Print DEmat
	for(int i=0; i<numstates; i++)
		{
		for(int j=0; j<numstates; j++)
			{
			//cout << DEmat(i,j);
			}
		//cout << "\n" ;
		}

	// Print dmat
	for(int i=0; i<numareas; i++)
		{
		for(int j=0; j<numareas; j++)
			{
			//cout << dmat(i,j) << " ";
			}
		//cout << "\n" ;
		}

	//cout << "num_nonzeros: " << num_nonzeros << "\n" ;
	//cout << DEmat.nrow() ;
	//cout << DEmat.ncol() ;
	
	/*
	// Copy the stored values to output flat arrays
	vector<int> ia_output[num_nonzeros];
	vector<int> ja_output[num_nonzeros];
	vector<float> a_output[num_nonzeros];
	for (int i=0; i<num_nonzeros; i++)
		{
		ia_output[i] = *(ia+i);
		ja_output[i] = *(ja+i);
		a_output[i] = *(a+i);
		}
	
	// Delete the dynamic pointer arrays
	delete [] ia;
	delete [] ja;
	delete [] a;
	*/
	// Create a list of vectors in Rcpp
	// http://stackoverflow.com/questions/3088650/how-do-i-create-a-list-of-vectors-in-rcpp
	Rcpp::List DEmat_COO = Rcpp::List::create(Rcpp::Named("ia") = ia,
                           Rcpp::Named("ja") = ja,
                           Rcpp::Named("a") = a);
	
	
	
	return DEmat_COO;
	}
} // end extern "C"






/* combine leftprobs and rightprobs through speciation model */
extern "C" 
{
SEXP cpp_calc_anclikes_sp(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP states_indices, SEXP s,
SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP maxent01y, SEXP
max_minsize_as_function_of_ancsize, SEXP Rsp_rowsums) {

	using namespace std;
	
	/* print if Rprintmat == 1 */
	//int printmat = Rcpp::as<int>(Rprintmat);
	//if (printmat == 1) {cout << "BEGIN probs_to_use:" << "\n";};
	
	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::NumericVector xa(leftprobs);
	Rcpp::NumericVector xb(rightprobs);
	Rcpp::List xl(states_indices);
	float xs = Rcpp::as<float>(s);
	float xv = Rcpp::as<float>(v);
	float xj = Rcpp::as<float>(j);
	float xy = Rcpp::as<float>(y);

	Rcpp::NumericMatrix dmatc(dmat);
	//int max_min_rangesize_c = Rcpp::as<int>(max_min_rangesize);
	Rcpp::IntegerVector max_min_rangesize_c(max_minsize_as_function_of_ancsize);
	
	// NOTE: THESE SHOULD BE a function of ancestral range size, f(ancsize), i.e. a matrix	
	Rcpp::NumericMatrix maxent01sub(maxent01s);
	Rcpp::NumericMatrix maxent01vic(maxent01v);
	Rcpp::NumericMatrix maxent01jump(maxent01j);
	Rcpp::NumericMatrix maxent01symp(maxent01y);
	Rcpp::NumericVector sp_rowsums(Rsp_rowsums);
	
	/* Get the sizes of the vectors */
	int n_xa = xa.size();
	int n_xb = xb.size();
	int numstates = xl.size();
	
	
	/* Store the states_indices as a vector of vectors */
	/* http://stackoverflow.com/questions/8159872/r-list-of-numeric-vectors-c-2d-array-with-rcpp */
	
	/* Empty vector for the states */
	vector< vector<int> > states_vecs;
	for(int i=0; i<numstates; i++)
		{
		/* Get the subvector */
		vector<int> tmp = Rcpp::as<vector<int> > (xl[i]);
		
		/* Push to the original list */
		states_vecs.push_back(tmp);
		}
	
	/* Make an int array of the range sizes of each state */
	
	/* This makes an vector of pointers, or something... (can't get returned to R by return) */
	/* int rangesizes [numstates]; */
	
	/* Define a vector of size rangesizes */
	Rcpp::NumericVector rangesizes (xa.size());
	/* Rcpp::NumericVector rangesizes( Rcpp::Vector<RTYPE>::Vector(const int&) tmpvec ); */
	for (int r=0; r<numstates; r++)
		{
		/* Rcpp::NumericVector tmp_areas( Rcpp::as<SEXP>(xl[r]) ); */
		/* Rcpp::NumericVector tmp_area = Rcpp::as<Rcpp::NumericVector> (xx[r]); */
		rangesizes[r] = states_vecs[r].size();
		/* rangesizes[r] = 1; */
		}
	
	/* The length of the output probabilities */
	Rcpp::NumericVector outprobs(n_xa);
	
	// The final result
	//int max_min_rangesize = *it5;
	
	// This prints the memory address
	//std::cout << max_min_rangesize << endl;
	// This prints the value(s) transferred
	//std::cout << max_min_rangesize_c << endl;
	
	
	//float testval = maxval(maxent01sub);
	
	/* l: go through each row of the ancestral states */
	for (int l = 0; l < numstates; l++)
		{
		/* Here, we don't need to keep track of all the cells, just add them up */
		float row_probval = 0;

		/* This is the ancestral state (the areas) for this row */
		vector<int> ancestral_areas = states_vecs[l];
		
		/* Get the ancestral range size */
		int ancsize = 0;
		ancsize = states_vecs[l].size();
		
		/* Go through the combinations of left and right descendant branches */
		// right branches = j
		for (int j = 0; j < n_xb; j++)
			{
			// left branches = i
			for (int i = 0; i < n_xa; i++)
				{
				float tmp_probval = 0;
				
				/* Get the lists of areas */
				/* Read each list item as a new SEXP (S expression) */
			
				/* areas on the left branch; and size */
				vector<int> lstate_areas( states_vecs[i] );
				int lsize = lstate_areas.size();
				
				/* areas on the right branch; and size */
				vector<int> rstate_areas( states_vecs[j] );
				int rsize = rstate_areas.size();

				// Shortcut
				// If the min rangesize is larger than max_min_rangesize_c, then
				// this cell has a value of 0 and you can continue
				int minsize_of_leftright = 0;
				minsize_of_leftright = std::min(lsize, rsize);
				// max_min_rangesize_c is a list of maximum sizes of the minimum (smaller) descendant
				// it's a function of the ancestral range size (ancsize)
				if (minsize_of_leftright > max_min_rangesize_c[(ancsize-1)])
					{
					// Just add zero
					row_probval += 0;
					//if (printmat == 1) {cout << "0" << " ";};
					continue; // go to next cell of the speciation matrix
					}
				
				// Perfect sympatry if ancestral_areas == lstate_areas && ancestral_areas == rstate_areas
				// check the sizes
				if ( (ancsize==lsize) && (ancsize==rsize) )
					{
					// check if they are all equal
					if ( all_ints_equal(lstate_areas, rstate_areas) == true )
						{
						// check if the DESCENDANT equals the ANCESTOR
						if ( all_ints_equal(lstate_areas, ancestral_areas) == true )
							{
							// sympatric prob (xy) * prob as a function of size * relprob leftright combination
							tmp_probval = xy * maxent01symp((ancsize-1), (ancsize-1)) * xa[i] * xb[j];
							
							row_probval += tmp_probval ;
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue; // go to next cell of the speciation matrix
							}
						}
					}
				
				
				//////////////////////////////////////////////////////////////////////
				// If one of the descendants is identical to the ancestor, 
				// (after we've excluded sympatry)
				// we can have jump dispersal or subset speciation
				//////////////////////////////////////////////////////////////////////
				//bool is_there_an_identical_daughter = false;
				
				// Does the left branch equal the ancestor?
				if ( (lsize == ancsize) && (all_ints_equal(lstate_areas, ancestral_areas)==true) )
					{
					// Is the left descendant smaller?
					if ( (rsize < ancsize) || (ancsize==1) )
						{
						// If so, then subset or jump speciation on the right branch
						
						// If the rstate_areas are ALL outside the ancestral range, it's jump/founder dispersal
						if (any_ints_equal(rstate_areas, ancestral_areas) == false)
							{
							// jump/founder event speciation on right branch
							// jump speciation on right branch
	
							// jump/founder event speciation on right branch
							// jump prob (xj) * prob as a function of size * relprob leftright combination
							// * dmatc((ancsize-1), (rsize-1)) = jump dispersal as a function of distance
							int try_jump_dispersal_based_on_dist = 1;
							float jprob_for_cell_based_on_distances = 0.0;
							if (try_jump_dispersal_based_on_dist == 1)
								{
								for (unsigned int ancarea_i=0; ancarea_i<ancestral_areas.size(); ancarea_i++)
									{
									for (unsigned int decarea_i=0; decarea_i<rstate_areas.size(); decarea_i++)
										{
										jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances +
	dmatc(ancestral_areas[ancarea_i], rstate_areas[decarea_i]);
										}
									}
								// Normalize (divide by the number of possible jump events) 
								// so that j parameter can be bigger
								float flt_aasize(ancestral_areas.size());
								float flt_rsize(rstate_areas.size());
								jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances / ( flt_aasize * flt_rsize);
								}
							else
								{
								jprob_for_cell_based_on_distances = 1.0;
								}
	
	
							tmp_probval = xj * maxent01jump((ancsize-1), (rsize-1)) * xa[i] * xb[j] *
	jprob_for_cell_based_on_distances;

							row_probval += tmp_probval;
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue;
							}
						else // otherwise it's partial overlap, or subset dispersal
							{
							// exclude partial overlap: ALL descendant states must match an ancestral state
							if (all_ints_found(lstate_areas, ancestral_areas) == true)
								{
								// subset speciation on right branch
								// subset prob (xs) * prob as a function of size * relprob leftright combination
								tmp_probval = xs * maxent01sub((ancsize-1), (rsize-1)) * xa[i] * xb[j];
								row_probval += tmp_probval;
								//if (printmat == 1) {cout << tmp_probval << " ";};
								continue;								
								} else {
								// partial overlap; probability 0 (for now)
								row_probval += 0;
								//if (printmat == 1) {cout << "0" << " ";};
								continue;								
								}						
							}
						}
					}
				// Does the right branch equal the ancestor?
				if ( (rsize == ancsize) && (all_ints_equal(rstate_areas, ancestral_areas)==true) )
					{
					// Is the right descendant smaller?
					if (( lsize < ancsize) || (ancsize==1) )
						{
						// If so, then subset or jump speciation on the left branch
						
						// If the lstate_areas are ALL outside the ancestral range, it's jump/founder dispersal
						if (any_ints_equal(lstate_areas, ancestral_areas) == false)
							{
							// jump/founder event speciation on left branch
							// jump speciation on left branch
	
							// jump/founder event speciation on left branch
							// jump prob (xj) * prob as a function of size * relprob leftleft combination
							// * dmatc((ancsize-1), (lsize-1)) = jump dispersal as a function of distance
							int try_jump_dispersal_based_on_dist = 1;
							float jprob_for_cell_based_on_distances = 0.0;
							if (try_jump_dispersal_based_on_dist == 1)
								{
								for (unsigned int ancarea_i=0; ancarea_i<ancestral_areas.size(); ancarea_i++)
									{
									for (unsigned int decarea_i=0; decarea_i<lstate_areas.size(); decarea_i++)
										{
										jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances +
	dmatc(ancestral_areas[ancarea_i], lstate_areas[decarea_i]);
										}
									}
								// Normalize (divide by the number of possible jump events) 
								// so that j parameter can be bigger
								jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances / (ancestral_areas.size() * lstate_areas.size());
								}
							else
								{
								jprob_for_cell_based_on_distances = 1.0;
								}
	
	
							tmp_probval = xj * maxent01jump((ancsize-1), (lsize-1)) * xa[i] * xb[j] *
	jprob_for_cell_based_on_distances;
								
							row_probval += tmp_probval;
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue;								
							}
						else // otherwise it's partial overlap, or subset dispersal
							{
							// exclude partial overlap: ALL descendant states must match an ancestral state
							if (all_ints_found(rstate_areas, ancestral_areas) == true)
								{
								// subset speciation on left branch
								// subset prob (xs) * prob as a function of size * relprob leftright combination
								tmp_probval = xs * maxent01sub((ancsize-1), (lsize-1)) * xa[i] * xb[j];
								row_probval += tmp_probval;
								//if (printmat == 1) {cout << tmp_probval << " ";};
								continue;								
								} else {
								// partial overlap; probability 0 (for now)
								row_probval += 0;
								//if (printmat == 1) {cout << tmp_probval << " ";};
								continue;								
								}						
							}
						}
					}
				
				// Vicariance can only happen when daughter states cover all of the ancestral states
				std::vector<int> combined_vector = merge_int_vectors(lstate_areas, rstate_areas);
				//std::cout << combined_vector.size() << endl;
				
				// The sizes must add up, and they must all match
				if ( (combined_vector.size() == ancestral_areas.size()) && (all_ints_equal(combined_vector,
ancestral_areas)==true)  )
					{
					int smaller_range_size = min(lsize, rsize);
					
					// vicariance speciation on right branch
					tmp_probval = xv * maxent01vic((ancsize-1), (smaller_range_size-1)) * xa[i] * xb[j];
					row_probval += tmp_probval;
					//if (printmat == 1) {cout << tmp_probval << " ";};
					continue;
					}
				
				
				// Remaining cases (hopefully none) are left blank
				row_probval += 0;
				//if (printmat == 1) {cout << "0" << " ";};

				/* int lstate_areas = xl[0].size(); */
				//float spmodel_prob = 1;
								
				/* Add up the probability of each combination */
				//row_probval += xa[i] * xb[j] * spmodel_prob * rangesizes[l];
				//row_probval = rangesize;
				}
			}
		// Add the probs for this row
		outprobs[l] = row_probval / sp_rowsums[l];
		//if (printmat == 1) {cout << "\n";};

		}
	
	/* Rcpp::NumericVector output(Rcpp::as<SEXP>(rangesizes)); */
	
	/* return the probability vector */
	//if (printmat == 1) {cout << "END probs to use\n";};
	//if (printmat == 1) {cout << "\n";};
	return outprobs;
	
	/* Return a single integer */
	/* return Rcpp::wrap(numstates); */
	//return Rcpp::wrap(*it);
	//return Rcpp::wrap(testval);
	}
} // end extern "C"







/* combine leftprobs and rightprobs through speciation model */
/* (sadly, we have to go through the sp matrix once just to get sums; still, we don't want to 
store the whole n*n^2 thing...) */
extern "C" 
{
SEXP cpp_calc_anclikes_sp_rowsums(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP states_indices,
SEXP s, SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP maxent01y, SEXP
max_minsize_as_function_of_ancsize) {

	using namespace std;

	/* print if Rprintmat == 1 */
	//int printmat = Rcpp::as<int>(Rprintmat);
	//if (printmat == 1) {cout << "BEGIN raw_probs_for_rowSums:" << "\n";};

	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::NumericVector xa(leftprobs);
	Rcpp::NumericVector xb(rightprobs);
	Rcpp::List xl(states_indices);
	float xs = Rcpp::as<float>(s);
	float xv = Rcpp::as<float>(v);
	float xj = Rcpp::as<float>(j);
	float xy = Rcpp::as<float>(y);

	Rcpp::NumericMatrix dmatc(dmat);
	//int max_min_rangesize_c = Rcpp::as<int>(max_min_rangesize);
	Rcpp::IntegerVector max_min_rangesize_c(max_minsize_as_function_of_ancsize);
	Rcpp::NumericMatrix maxent01sub(maxent01s);
	Rcpp::NumericMatrix maxent01vic(maxent01v);
	Rcpp::NumericMatrix maxent01jump(maxent01j);
	Rcpp::NumericMatrix maxent01symp(maxent01y);
	
	/* Get the sizes of the vectors */
	int n_xa = xa.size();
	int n_xb = xb.size();
	int numstates = xl.size();
	
	
	/* Store the states_indices as a vector of vectors */
	/* http://stackoverflow.com/questions/8159872/r-list-of-numeric-vectors-c-2d-array-with-rcpp */
	
	/* Empty vector for the states */
	vector< vector<int> > states_vecs;
	for(int i=0; i<numstates; i++)
		{
		/* Get the subvector */
		vector<int> tmp = Rcpp::as<vector<int> > (xl[i]);
		
		/* Push to the original list */
		states_vecs.push_back(tmp);
		}
	
	/* Make an int array of the range sizes of each state */
	
	/* This makes an vector of pointers, or something... (can't get returned to R by return) */
	/* int rangesizes [numstates]; */
	
	/* Define a vector of size rangesizes */
	Rcpp::NumericVector rangesizes (xa.size());
	/* Rcpp::NumericVector rangesizes( Rcpp::Vector<RTYPE>::Vector(const int&) tmpvec ); */
	for (int r=0; r<numstates; r++)
		{
		/* Rcpp::NumericVector tmp_areas( Rcpp::as<SEXP>(xl[r]) ); */
		/* Rcpp::NumericVector tmp_area = Rcpp::as<Rcpp::NumericVector> (xx[r]); */
		rangesizes[r] = states_vecs[r].size();
		/* rangesizes[r] = 1; */
		}
	
	/* The length of the output probabilities */
	Rcpp::NumericVector sp_rowsums(n_xa);
	
	// The final result
	//int max_min_rangesize = *it5;
	
	// This prints the memory address
	//std::cout << max_min_rangesize << endl;
	// This prints the value(s) transferred
	//std::cout << max_min_rangesize_c << endl;
	
	
	//float testval = maxval(maxent01sub);
	
	/* l: go through each row of the ancestral states */
	for (int l = 0; l < numstates; l++)
		{
		/* Here, we don't need to keep track of all the cells, just add them up */
		float row_probval = 0;

		/* This is the ancestral state (the areas) for this row */
		vector<int> ancestral_areas = states_vecs[l];
		
		/* Get the ancestral range size */
		int ancsize = 0;
		ancsize = states_vecs[l].size();
		
		/* Go through the combinations of left and right descendant branches */
		// right branches = j
		for (int j = 0; j < n_xb; j++)
			{
			// left branches = i
			for (int i = 0; i < n_xa; i++)
				{
				float tmp_probval = 0;
				
				/* Get the lists of areas */
				/* Read each list item as a new SEXP (S expression) */
			
				/* areas on the left branch; and size */
				vector<int> lstate_areas( states_vecs[i] );
				int lsize = lstate_areas.size();
				
				/* areas on the right branch; and size */
				vector<int> rstate_areas( states_vecs[j] );
				int rsize = rstate_areas.size();

				// Shortcut
				// If the min rangesize is larger than max_min_rangesize_c, then
				// this cell has a value of 0 and you can continue
				int minsize_of_leftright = 0;
				minsize_of_leftright = std::min(lsize, rsize);
				
				// max_min_rangesize_c is a list of maximum sizes of the minimum (smaller) descendant
				// it's a function of the ancestral range size (ancsize)
				if (minsize_of_leftright > max_min_rangesize_c[(ancsize-1)])
					{
					// Just add zero
					row_probval += 0;
					//if (printmat == 1) {cout << "0" << " ";};
					continue; // go to next cell of the speciation matrix
					}
				
				// Perfect sympatry if ancestral_areas == lstate_areas && ancestral_areas == rstate_areas
				// check the sizes
				if ( (ancsize==lsize) && (ancsize==rsize) )
					{
					// check if they are all equal
					if ( all_ints_equal(lstate_areas, rstate_areas) == true )
						{
						// check if the DESCENDANT equals the ANCESTOR
						if ( all_ints_equal(lstate_areas, ancestral_areas) == true )
							{
							// sympatric prob (xy) * prob as a function of size * relprob leftright combination
							tmp_probval = xy * maxent01symp((ancsize-1), (ancsize-1)) * 1 * 1;
							
							row_probval += tmp_probval ;
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue; // go to next cell of the speciation matrix
							}
						}
					}
				
				
				//////////////////////////////////////////////////////////////////////
				// If one of the descendants is identical to the ancestor, 
				// (after we've excluded sympatry)
				// we can have jump dispersal or subset speciation
				//////////////////////////////////////////////////////////////////////
				//bool is_there_an_identical_daughter = false;
				
				// Does the left branch equal the ancestor?
				// l = row of speciation matrix
				// i = left branch (change quickly)
				// j = right branch (change slowly)
				int printthis = 0;
				if (printthis == 1) {
				if ( l==0 ) {
				if ( (i==1) && (j==0) )
					{
					//cout << "\nThis should be a jump event" << endl;
					

					//cout << "ancestral_areas:" << endl;
					for (unsigned int m=0; m<ancestral_areas.size(); m++)
						{
						//cout << ancestral_areas[m] << " ";
						}
					//cout << endl;

					//cout << "lstate_areas:" << endl;
					for (unsigned int m=0; m<lstate_areas.size(); m++)
						{
						//cout << lstate_areas[m] << " ";
						}
					//cout << endl;

					//cout << "rstate_areas:" << endl;
					for (unsigned int m=0; m<rstate_areas.size(); m++)
						{
						//cout << rstate_areas[m] << " ";
						}
					//cout << endl;

					//cout << "ancsize	lsize	rsize:" << endl;
					//cout << ancsize << " " << lsize << " " << rsize << endl;
					//cout << "\nleft side:" << endl;
					//cout << (lsize == ancsize) << endl;
					//cout << all_ints_equal(lstate_areas, ancestral_areas) << endl;
					//cout << ((lsize == ancsize) && (all_ints_equal(lstate_areas, ancestral_areas)==true)) << endl;
					//cout << ((lsize < ancsize) || (ancsize==1)) << endl;
					//cout << any_ints_equal(rstate_areas, ancestral_areas) << endl;
					//cout << all_ints_found(lstate_areas, ancestral_areas) << endl;
					
					//cout << "right side:" << endl;
					//cout << (rsize == ancsize) << endl;
					//cout << all_ints_equal(rstate_areas, ancestral_areas) << endl;
					//cout << ( (rsize == ancsize) && (all_ints_equal(rstate_areas, ancestral_areas)==true) ) << endl;
					//cout << ((rsize < ancsize) || (ancsize==1)) << endl;
					//cout << any_ints_equal(lstate_areas, ancestral_areas) << endl;
					//cout << all_ints_found(rstate_areas, ancestral_areas) << endl;
					//cout << " " << endl;					
					} //end if
					} //end if
					} //end printthis
				if ( (lsize == ancsize) && (all_ints_equal(lstate_areas, ancestral_areas)==true) )
					{
					// Is the left descendant smaller?  (Or, ancsize=1)
					if ((rsize < ancsize) || (ancsize==1))
						{
						// If so, then subset or jump speciation on the right branch
						
						// If the rstate_areas are ALL outside the ancestral range, it's jump/founder dispersal
						if (any_ints_equal(rstate_areas, ancestral_areas) == false)
							{
							// jump/founder event speciation on right branch
							// jump prob (xj) * prob as a function of size * relprob leftright combination
							
							// * dmatc((ancsize-1), (rsize-1)) = jump dispersal as a function of distance
							int try_jump_dispersal_based_on_dist = 1;
							float jprob_for_cell_based_on_distances = 0.0;
							if (try_jump_dispersal_based_on_dist == 1)
								{
								for (unsigned int ancarea_i=0; ancarea_i<ancestral_areas.size(); ancarea_i++)
									{
									for (unsigned int decarea_i=0; decarea_i<rstate_areas.size(); decarea_i++)
										{
										jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances +
	dmatc(ancestral_areas[ancarea_i], rstate_areas[decarea_i]);
										}
									}
								// Normalize (divide by the number of possible jump events) 
								// so that j parameter can be bigger
								float flt_aasize(ancestral_areas.size());
								float flt_rsize(rstate_areas.size());
								jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances / ( flt_aasize * flt_rsize);
								}
							else
								{
								jprob_for_cell_based_on_distances = 1.0;
								}
	
	
							tmp_probval = xj * maxent01jump((ancsize-1), (rsize-1)) * 1 * 1 *
	jprob_for_cell_based_on_distances;
							
							row_probval += tmp_probval;
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue;
							}
						else // otherwise it's partial overlap, or subset dispersal
							{
							// exclude partial overlap: ALL descendant states must match an ancestral state
							if (all_ints_found(lstate_areas, ancestral_areas) == true)
								{
								// subset speciation on right branch
								// subset prob (xs) * prob as a function of size * relprob leftright combination
								tmp_probval = xs * maxent01sub((ancsize-1), (rsize-1)) * 1 * 1;
								row_probval += tmp_probval;
								//if (printmat == 1) {cout << tmp_probval << " ";};
								continue;								
								} else {
								// partial overlap; probability 0 (for now)
								row_probval += 0;
								//if (printmat == 1) {cout << "0" << " ";};
								continue;								
								}						
							}
						}
					}
				// Does the right branch equal the ancestor?
				if ( (rsize == ancsize) && (all_ints_equal(rstate_areas, ancestral_areas)==true) )
					{
					// Is the right descendant smaller?  (Or, ancsize=1)
					if ((lsize < ancsize) || (ancsize==1))
						{
						// If so, then subset or jump speciation on the left branch
						
						// If the lstate_areas are ALL outside the ancestral range, it's jump/founder dispersal
						if (any_ints_equal(lstate_areas, ancestral_areas) == false)
							{
							// jump speciation on left branch
	
							// jump/founder event speciation on left branch
							// jump prob (xj) * prob as a function of size * relprob leftleft combination
							// * dmatc((ancsize-1), (lsize-1)) = jump dispersal as a function of distance
							int try_jump_dispersal_based_on_dist = 1;
							float jprob_for_cell_based_on_distances = 0.0;
							if (try_jump_dispersal_based_on_dist == 1)
								{
								for (unsigned int ancarea_i=0; ancarea_i<ancestral_areas.size(); ancarea_i++)
									{
									for (unsigned int decarea_i=0; decarea_i<lstate_areas.size(); decarea_i++)
										{
										jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances +
	dmatc(ancestral_areas[ancarea_i], lstate_areas[decarea_i]);
										}
									}
								// Normalize (divide by the number of possible jump events) 
								// so that j parameter can be bigger
								jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances / (ancestral_areas.size() * lstate_areas.size());
								}
							else
								{
								jprob_for_cell_based_on_distances = 1.0;
								}
	
	
							tmp_probval = xj * maxent01jump((ancsize-1), (lsize-1)) * 1 * 1 *
	jprob_for_cell_based_on_distances;

							row_probval += tmp_probval;
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue;								
							}
						else // otherwise it's partial overlap, or subset dispersal
							{
							// exclude partial overlap: ALL descendant states must match an ancestral state
							if (all_ints_found(rstate_areas, ancestral_areas) == true)
								{
								// subset speciation on left branch
								// subset prob (xs) * prob as a function of size * relprob leftright combination
								tmp_probval = xs * maxent01sub((ancsize-1), (lsize-1)) * 1 * 1;
								row_probval += tmp_probval;
								//if (printmat == 1) {cout << tmp_probval << " ";};
								continue;								
								} else {
								// partial overlap; probability 0 (for now)
								row_probval += 0;
								//if (printmat == 1) {cout << "0" << " ";};
								continue;								
								}						
							}
						}
					}
				
				// Vicariance can only happen when daughter states cover all of the ancestral states
				std::vector<int> combined_vector = merge_int_vectors(lstate_areas, rstate_areas);
				//std::cout << combined_vector.size() << endl;
				
				// The sizes must add up, and they must all match
				if ( (combined_vector.size() == ancestral_areas.size()) && (all_ints_equal(combined_vector,
ancestral_areas)==true)  )
					{
					int smaller_range_size = min(lsize, rsize);
					
					// vicariance speciation on right branch
					tmp_probval = xv * maxent01vic((ancsize-1), (smaller_range_size-1)) * 1 * 1;
					row_probval += tmp_probval;
					//if (printmat == 1) {cout << tmp_probval << " ";};
					continue;
					}
				
				
				// Remaining cases (hopefully none) are left blank
				row_probval += 0;
				//if (printmat == 1) {cout << "0" << " ";};
				
				/* int lstate_areas = xl[0].size(); */
				//float spmodel_prob = 1;
								
				/* Add up the probability of each combination */
				//row_probval += 1 * 1 * spmodel_prob * rangesizes[l];
				//row_probval = rangesize;
				}
			}
		// Add the probs for this row
		sp_rowsums[l] = row_probval;

		//if (printmat == 1) {cout << "\n";};
		}
	
	/* Rcpp::NumericVector output(Rcpp::as<SEXP>(rangesizes)); */
	
	/* return the probability vector */
	//if (printmat == 1) {cout << "END raw_probs_for_rowSums\n";};
	//if (printmat == 1) {cout << "\n";};
	
	/*
	if (printmat == 1) {
		//cout << "BEGIN sp_rowsums\n";
		for (int i=0; i<sp_rowsums.size(); i++)
			{
			//cout << sp_rowsums[i] << "\n";		
			}
		//cout << "END sp_rowsums\n";
		};
	*/
	//if (printmat == 1) {cout << "\n";};
	
	
	return sp_rowsums;
	
	/* Return a single integer */
	/* return Rcpp::wrap(numstates); */
	//return Rcpp::wrap(*it);
	//return Rcpp::wrap(testval);
	}

} // end extern "C"










/* Originally, I just combined leftprobs and rightprobs through speciation model */
/* 						*/
/* E.g.: 				*/
/* 		A,A		A,B		A,AB	B,A		B,B		B,AB	AB,A	AB,B	AB,AB				*/
/* A	0.1																0.9					*/
/* B	0.2						0.8															*/
/* AB	0.3												0.7									*/
/* 																							*/
/* It became faster to at least calculate the rowsums once(sadly, we have to go through the */
/* sp matrix once just to get sums; still, we don't want to re-do this.						*/
/* 																							*/
/* It would be faster to calculate the conditional probabilities ONCE and store them, but	*/
/* we don't want to store the whole n*n^2 thing.  So, store a COO-formatted sp matrix of	*/
/* of non-zero values, and their coordinates in i (left) and j (right) states.				*/
/* */
/* */
/* */
extern "C" 
{
SEXP cpp_calc_anclikes_sp_COOprobs(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP states_indices,
SEXP s, SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP maxent01y, SEXP
max_minsize_as_function_of_ancsize) {

	using namespace std;

	float min_precision = 1e-30;

	/* print if Rprintmat == 1 */
	//int printmat = Rcpp::as<int>(Rprintmat);
	//if (printmat == 1) {cout << "BEGIN cpp_calc_anclikes_sp_COOprobs()" << "\n";};

	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::NumericVector xa(leftprobs);
	Rcpp::NumericVector xb(rightprobs);
	Rcpp::List xl(states_indices);
	float xs = Rcpp::as<float>(s);
	float xv = Rcpp::as<float>(v);
	float xj = Rcpp::as<float>(j);
	float xy = Rcpp::as<float>(y);

	Rcpp::NumericMatrix dmatc(dmat);
	//int max_min_rangesize_c = Rcpp::as<int>(max_min_rangesize);
	Rcpp::IntegerVector max_min_rangesize_c(max_minsize_as_function_of_ancsize);
	Rcpp::NumericMatrix maxent01sub(maxent01s);
	Rcpp::NumericMatrix maxent01vic(maxent01v);
	Rcpp::NumericMatrix maxent01jump(maxent01j);
	Rcpp::NumericMatrix maxent01symp(maxent01y);
	
	/* Get the sizes of the vectors */
	int n_xa = xa.size();
	int n_xb = xb.size();
	int numstates = xl.size();
	
	/* Setup the vectors for row, i (left side state) and j (right side state) */
	vector< vector<int> > nonzero_i;
	vector< vector<int> > nonzero_j;
	
	/* Setup the vectors for row, column probs */
	vector< vector<float> > nonzero_probs;
	
	/* Store the states_indices as a vector of vectors */
	/* http://stackoverflow.com/questions/8159872/r-list-of-numeric-vectors-c-2d-array-with-rcpp */
	
	/* Empty vector for the states */
	vector< vector<int> > states_vecs;
	for(int i=0; i<numstates; i++)
		{
		/* Get the subvector */
		vector<int> tmp = Rcpp::as<vector<int> > (xl[i]);
		
		/* Push to the original list */
		states_vecs.push_back(tmp);
		}
	
	/* Make an int array of the range sizes of each state */
	
	/* This makes an vector of pointers, or something... (can't get returned to R by return) */
	/* int rangesizes [numstates]; */
	
	/* Define a vector of size rangesizes */
	Rcpp::NumericVector rangesizes (xa.size());
	/* Rcpp::NumericVector rangesizes( Rcpp::Vector<RTYPE>::Vector(const int&) tmpvec ); */
	for (int r=0; r<numstates; r++)
		{
		/* Rcpp::NumericVector tmp_areas( Rcpp::as<SEXP>(xl[r]) ); */
		/* Rcpp::NumericVector tmp_area = Rcpp::as<Rcpp::NumericVector> (xx[r]); */
		rangesizes[r] = states_vecs[r].size();
		/* rangesizes[r] = 1; */
		}
	
	/* The length of the output probabilities */
	Rcpp::NumericVector sp_rowsums(n_xa);
	
	// The final result
	//int max_min_rangesize = *it5;
	
	// This prints the memory address
	//std::cout << max_min_rangesize << endl;
	// This prints the value(s) transferred
	//std::cout << max_min_rangesize_c << endl;
	
	
	//float testval = maxval(maxent01sub);
	
	/* l: go through each row of the ancestral states */
	for (int l = 0; l < numstates; l++)
		{
		/* Get the subvector of nonzero descendant for left and right branches (i and j) */
		vector<int> tmp_inums;
		vector<int> tmp_jnums;

		/* Get the subvector of nonzero column probabilities */
		vector<float> tmp_colprobs;

		/* Here, we don't need to keep track of all the cells, just add them up */
		float row_probval = 0;

		/* This is the ancestral state (the areas) for this row */
		vector<int> ancestral_areas = states_vecs[l];
		
		/* Get the ancestral range size */
		int ancsize = 0;
		ancsize = states_vecs[l].size();
		
		/* Go through the combinations of left and right descendant branches */
		// right branches = j
		for (int j = 0; j < n_xb; j++)
			{
			// left branches = i
			for (int i = 0; i < n_xa; i++)
				{
				float tmp_probval = 0;
				
				/* Get the lists of areas */
				/* Read each list item as a new SEXP (S expression) */
			
				/* areas on the left branch; and size */
				vector<int> lstate_areas( states_vecs[i] );
				int lsize = lstate_areas.size();
				
				/* areas on the right branch; and size */
				vector<int> rstate_areas( states_vecs[j] );
				int rsize = rstate_areas.size();

				// Shortcut
				// If the min rangesize is larger than max_min_rangesize_c, then
				// this cell has a value of 0 and you can continue
				int minsize_of_leftright = 0;
				minsize_of_leftright = std::min(lsize, rsize);
				
				// max_min_rangesize_c is a list of maximum sizes of the minimum (smaller) descendant
				// it's a function of the ancestral range size (ancsize)
				if (minsize_of_leftright > max_min_rangesize_c[(ancsize-1)])
					{
					// Just add zero
					row_probval += 0;
					//if (printmat == 1) {cout << "0" << " ";};
					continue; // go to next cell of the speciation matrix
					}
				
				// Perfect sympatry if ancestral_areas == lstate_areas && ancestral_areas == rstate_areas
				// check the sizes
				if ( (ancsize==lsize) && (ancsize==rsize) )
					{
					// check if they are all equal
					if ( all_ints_equal(lstate_areas, rstate_areas) == true )
						{
						// check if the DESCENDANT equals the ANCESTOR
						if ( all_ints_equal(lstate_areas, ancestral_areas) == true )
							{
							// sympatric prob (xy) * prob as a function of size * relprob leftright combination
							tmp_probval = xy * maxent01symp((ancsize-1), (ancsize-1)) * 1 * 1;
							row_probval += tmp_probval ;
							
							// Add the i index (left side state) and j index (right side state) to the list for
							// this ancestral row
							if (std::abs(tmp_probval) > min_precision)
								{
								tmp_inums.push_back(i);
								tmp_jnums.push_back(j);
								tmp_colprobs.push_back(tmp_probval);
								}
							
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue; // go to next cell of the speciation matrix
							}
						}
					}
				
				
				//////////////////////////////////////////////////////////////////////
				// If one of the descendants is identical to the ancestor, 
				// (after we've excluded sympatry)
				// we can have jump dispersal or subset speciation
				//////////////////////////////////////////////////////////////////////
				//bool is_there_an_identical_daughter = false;
				
				// Does the left branch equal the ancestor?
				// l = row of speciation matrix
				// i = left branch (change quickly)
				// j = right branch (change slowly)
				int printthis = 0;
				if (printthis == 1) {
				if ( l==0 ) {
				if ( (i==1) && (j==0) )
					{
					//cout << "\nThis should be a jump event" << endl;
					

					//cout << "ancestral_areas:" << endl;
					for (unsigned int m=0; m<ancestral_areas.size(); m++)
						{
						//cout << ancestral_areas[m] << " ";
						}
					//cout << endl;

					//cout << "lstate_areas:" << endl;
					for (unsigned int m=0; m<lstate_areas.size(); m++)
						{
						//cout << lstate_areas[m] << " ";
						}
					//cout << endl;

					//cout << "rstate_areas:" << endl;
					for (unsigned int m=0; m<rstate_areas.size(); m++)
						{
						//cout << rstate_areas[m] << " ";
						}
					//cout << endl;

					//cout << "ancsize	lsize	rsize:" << endl;
					//cout << ancsize << " " << lsize << " " << rsize << endl;
					//cout << "\nleft side:" << endl;
					//cout << (lsize == ancsize) << endl;
					//cout << all_ints_equal(lstate_areas, ancestral_areas) << endl;
					//cout << ((lsize == ancsize) && (all_ints_equal(lstate_areas, ancestral_areas)==true)) << endl;
					//cout << ((lsize < ancsize) || (ancsize==1)) << endl;
					//cout << any_ints_equal(rstate_areas, ancestral_areas) << endl;
					//cout << all_ints_found(lstate_areas, ancestral_areas) << endl;
					
					//cout << "right side:" << endl;
					//cout << (rsize == ancsize) << endl;
					//cout << all_ints_equal(rstate_areas, ancestral_areas) << endl;
					//cout << ( (rsize == ancsize) && (all_ints_equal(rstate_areas, ancestral_areas)==true) ) << endl;
					//cout << ((rsize < ancsize) || (ancsize==1)) << endl;
					//cout << any_ints_equal(lstate_areas, ancestral_areas) << endl;
					//cout << all_ints_found(rstate_areas, ancestral_areas) << endl;
					//cout << " " << endl;					
					} //end if
					} //end if
					} //end printthis
				if ( (lsize == ancsize) && (all_ints_equal(lstate_areas, ancestral_areas)==true) )
					{
					// Is the left descendant smaller?  (Or, ancsize=1)
					if ((rsize < ancsize) || (ancsize==1))
						{
						// If so, then subset or jump speciation on the right branch
						
						// If the rstate_areas are ALL outside the ancestral range, it's jump/founder dispersal
						if (any_ints_equal(rstate_areas, ancestral_areas) == false)
							{
							// jump speciation on right branch
							// jump/founder event speciation on right branch
							// jump prob (xj) * prob as a function of size * relprob leftright combination
							// * dmatc((ancsize-1), (rsize-1)) = jump dispersal as a function of distance
							int try_jump_dispersal_based_on_dist = 1;
							float jprob_for_cell_based_on_distances = 0.0;
							if (try_jump_dispersal_based_on_dist == 1)
								{
								for (unsigned int ancarea_i=0; ancarea_i<ancestral_areas.size(); ancarea_i++)
									{
									for (unsigned int decarea_i=0; decarea_i<rstate_areas.size(); decarea_i++)
										{
										jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances +
	dmatc(ancestral_areas[ancarea_i], rstate_areas[decarea_i]);
										}
									}
								// Normalize (divide by the number of possible jump events) 
								// so that j parameter can be bigger
								float flt_aasize(ancestral_areas.size());
								float flt_rsize(rstate_areas.size());
								jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances / ( flt_aasize * flt_rsize);
								}
							else
								{
								jprob_for_cell_based_on_distances = 1.0;
								}
	
	
							tmp_probval = xj * maxent01jump((ancsize-1), (rsize-1)) * 1 * 1 *
	jprob_for_cell_based_on_distances;

							row_probval += tmp_probval;

							// Add the i index (left side state) and j index (right side state) to the list
							// for this ancestral row
							if (std::abs(tmp_probval) > min_precision)
								{
								tmp_inums.push_back(i);
								tmp_jnums.push_back(j);
								tmp_colprobs.push_back(tmp_probval);
								}

							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue;
							}
						else // otherwise it's partial overlap, or subset dispersal
							{
							// exclude partial overlap: ALL descendant states must match an ancestral state
							if (all_ints_found(lstate_areas, ancestral_areas) == true)
								{
								// subset speciation on right branch
								// subset prob (xs) * prob as a function of size * relprob leftright combination
								tmp_probval = xs * maxent01sub((ancsize-1), (rsize-1)) * 1 * 1;
								row_probval += tmp_probval;

								// Add the i index (left side state) and j index (right side state) to the list for
								// this ancestral row
								if (std::abs(tmp_probval) > min_precision)
									{
									tmp_inums.push_back(i);
									tmp_jnums.push_back(j);
									tmp_colprobs.push_back(tmp_probval);
									}

								//if (printmat == 1) {cout << tmp_probval << " ";};
								continue;								
								} else {
								// partial overlap; probability 0 (for now)
								row_probval += 0;
								//if (printmat == 1) {cout << "0" << " ";};
								continue;								
								}						
							}
						}
					}
				// Does the right branch equal the ancestor?
				if ( (rsize == ancsize) && (all_ints_equal(rstate_areas, ancestral_areas)==true) )
					{
					// Is the right descendant smaller?  (Or, ancsize=1)
					if ((lsize < ancsize) || (ancsize==1))
						{
						// If so, then subset or jump speciation on the left branch
						
						// If the lstate_areas are ALL outside the ancestral range, it's jump/founder dispersal
						if (any_ints_equal(lstate_areas, ancestral_areas) == false)
							{
							// jump speciation on left branch
							// jump/founder event speciation on left branch
							// jump prob (xj) * prob as a function of size * relprob leftleft combination
							// * dmatc((ancsize-1), (lsize-1)) = jump dispersal as a function of distance
							int try_jump_dispersal_based_on_dist = 1;
							float jprob_for_cell_based_on_distances = 0.0;
							if (try_jump_dispersal_based_on_dist == 1)
								{
								for (unsigned int ancarea_i=0; ancarea_i<ancestral_areas.size(); ancarea_i++)
									{
									for (unsigned int decarea_i=0; decarea_i<lstate_areas.size(); decarea_i++)
										{
										jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances +
	dmatc(ancestral_areas[ancarea_i], lstate_areas[decarea_i]);
										}
									}
								// Normalize (divide by the number of possible jump events) 
								// so that j parameter can be bigger
								jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances / (ancestral_areas.size() * lstate_areas.size());
								}
							else
								{
								jprob_for_cell_based_on_distances = 1.0;
								}
	
	
							tmp_probval = xj * maxent01jump((ancsize-1), (lsize-1)) * 1 * 1 *
	jprob_for_cell_based_on_distances;

							row_probval += tmp_probval;

							// Add the i index (left side state) and j index (right side state) to the list
							// for this ancestral row
							if (std::abs(tmp_probval) > min_precision)
								{
								tmp_inums.push_back(i);
								tmp_jnums.push_back(j);
								tmp_colprobs.push_back(tmp_probval);
								}
							
							//if (printmat == 1) {cout << tmp_probval << " ";};
							continue;								
							}
						else // otherwise it's partial overlap, or subset dispersal
							{
							// exclude partial overlap: ALL descendant states must match an ancestral state
							if (all_ints_found(rstate_areas, ancestral_areas) == true)
								{
								// subset speciation on left branch
								// subset prob (xs) * prob as a function of size * relprob leftright combination
								tmp_probval = xs * maxent01sub((ancsize-1), (lsize-1)) * 1 * 1;
								row_probval += tmp_probval;

								// Add the i index (left side state) and j index (right side state) to the list
								// for this ancestral row
								if (std::abs(tmp_probval) > min_precision)
									{
									tmp_inums.push_back(i);
									tmp_jnums.push_back(j);
									tmp_colprobs.push_back(tmp_probval);
									}

								//if (printmat == 1) {cout << tmp_probval << " ";};
								continue;								
								} else {
								// partial overlap; probability 0 (for now)
								row_probval += 0;
								//if (printmat == 1) {cout << "0" << " ";};
								continue;								
								}						
							}
						}
					}
				
				// Vicariance can only happen when daughter states cover all of the ancestral states
				std::vector<int> combined_vector = merge_int_vectors(lstate_areas, rstate_areas);
				//std::cout << combined_vector.size() << endl;
				
				// The sizes must add up, and they must all match
				if ( (combined_vector.size() == ancestral_areas.size()) && (all_ints_equal(combined_vector,
ancestral_areas)==true)  )
					{
					int smaller_range_size = min(lsize, rsize);
					
					// vicariance speciation on right branch
					tmp_probval = xv * maxent01vic((ancsize-1), (smaller_range_size-1)) * 1 * 1;
					row_probval += tmp_probval;
					
					// Add the i index (left side state) and j index (right side state) to the list
					// for this ancestral row
					if (std::abs(tmp_probval) > min_precision)
						{
						tmp_inums.push_back(i);
						tmp_jnums.push_back(j);
						tmp_colprobs.push_back(tmp_probval);
						}
					
					//if (printmat == 1) {cout << tmp_probval << " ";};
					continue;
					}
				
				
				// Remaining cases (hopefully none) are left blank
				row_probval += 0;
				//if (printmat == 1) {cout << "0" << " ";};
				
				/* int lstate_areas = xl[0].size(); */
				//float spmodel_prob = 1;
								
				/* Add up the probability of each combination */
				//row_probval += 1 * 1 * spmodel_prob * rangesizes[l];
				//row_probval = rangesize;
				}
			}
		// Add the probs for this row
		sp_rowsums[l] = row_probval;
		
		// Add each row to its list of rows
		nonzero_i.push_back(tmp_inums);
		nonzero_j.push_back(tmp_jnums);
		nonzero_probs.push_back(tmp_colprobs);
		
		// Delete these three sublists
		//delete tmp_inums;
		//delete tmp_jnums;
		//delete tmp_colprobs;
		
		
		//if (printmat == 1) {cout << "\n";};
		}
	
	/* Rcpp::NumericVector output(Rcpp::as<SEXP>(rangesizes)); */
	
	/* return the probability vector */
	//if (printmat == 1) {cout << "END cpp_calc_anclikes_sp_COOprobs()\n";};
	//if (printmat == 1) {cout << "\n";};
	
	/*
	if (printmat == 1) {
		//cout << "BEGIN sp_rowsums\n";
		for (int i=0; i<sp_rowsums.size(); i++)
			{
			//cout << sp_rowsums[i] << "\n";		
			}
		//cout << "END sp_rowsums\n";
		};
	*/
	//if (printmat == 1) {cout << "\n";};
	
	
	//return sp_rowsums;
	return(Rcpp::List::create(nonzero_i, nonzero_j, nonzero_probs));
	
	
	/* Return a single integer */
	/* return Rcpp::wrap(numstates); */
	//return Rcpp::wrap(*it);
	//return Rcpp::wrap(testval);
	}

} // end extern "C"



/* Here, we will try to speed this up, by targeted formation of the COOprobs list. To wit:*/
/* */
/* */
/*
Possible faster speciation probability calculation

Take states
order by range sizes
make a list of start indices for size 1, size 2, size 3, etc.

Subset speciation:
1. iterate through ancestors
- for each ancestor, assume left branch = ancestor
- iterate through right branches of appropriate sizes

1. Iterate through ancestors=left branches
2. For each ancestor, the sisters are all subsets of size 1, size 2, etc.
3. flip list for right ancestors being widespread


Jump speciation:
1. iterate through ancestors
- for each ancestor, assume left branch = ancestor

2. For each ancestor, the sisters are all subsets of size 1, size 2, etc. OUTSIDE of ancestral range

3. flip list for right ancestors being widespread


Vicariance speciation:
1. iterate through ancestors
- for each ancestor, assume left branch = bigger
- for each ancestor, assume right branch = smaller

2. Go through each of the possible right branches of size 1, size 2, etc. INSIDE of ancestral range
2a. Left branch gets the opposite

3. flip list for right ancestors being more widespread
*/
/* */
/* */
/* Originally, I just combined leftprobs and rightprobs through speciation model */
/* 						*/
/* E.g.: 				*/
/* 		A,A		A,B		A,AB	B,A		B,B		B,AB	AB,A	AB,B	AB,AB				*/
/* A	0.1																0.9					*/
/* B	0.2						0.8															*/
/* AB	0.3												0.7									*/
/* 																							*/
/* It became faster to at least calculate the rowsums once(sadly, we have to go through the */
/* sp matrix once just to get sums; still, we don't want to re-do this.						*/
/* 																							*/
/* It would be faster to calculate the conditional probabilities ONCE and store them, but	*/
/* we don't want to store the whole n*n^2 thing.  So, store a COO-formatted sp matrix of	*/
/* of non-zero values, and their coordinates in i (left) and j (right) states.				*/
/* */
/* */
/* */
extern "C" 
{
SEXP cpp_calc_anclikes_sp_COOweights_faster(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP
states_indices, SEXP s, SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP
maxent01y, SEXP max_minsize_as_function_of_ancsize) {

	using namespace std;
	//std::cout << "A0" << endl;

	float min_precision = 1e-30;
	min_precision = min_precision + 0.0; 	// just so that it's used
	/* print if Rprintmat == 1 */
	
	//printmat not needed now that we've commented out all printmats
	//int printmat = Rcpp::as<int>(Rprintmat);
	//if (printmat == 1) {cout << "BEGIN cpp_calc_anclikes_sp_COOprobs_fast()" << "\n";};

	/* Define the numeric vectors and put in the data into C++ from R */
	/* These should be 1s when just doing the initial setup */
	Rcpp::NumericVector xa(leftprobs);
	Rcpp::NumericVector xb(rightprobs);
	/* The above should be 1s when just doing the initial setup */

	Rcpp::List xl(states_indices);
	
	
	
	float xs = Rcpp::as<float>(s);
	float xv = Rcpp::as<float>(v);
	float xj = Rcpp::as<float>(j);
	float xy = Rcpp::as<float>(y);

	Rcpp::NumericMatrix dmatc(dmat);
	//int max_min_rangesize_c = Rcpp::as<int>(max_min_rangesize);
	Rcpp::IntegerVector max_min_rangesize_c(max_minsize_as_function_of_ancsize);
	Rcpp::NumericMatrix maxent01sub(maxent01s);
	Rcpp::NumericMatrix maxent01vic(maxent01v);
	Rcpp::NumericMatrix maxent01jump(maxent01j);
	Rcpp::NumericMatrix maxent01symp(maxent01y);
	
	//std::cout << "A1" << endl;
	
	/* Get the sizes of the vectors */
	int n_xa = xa.size();
	//int n_xb = xb.size();	// not used
	int numstates = xl.size();

	
	/* Setup the vectors for row, i (left side state) and j (right side state) */
	vector<int> nonzero_i;
	vector<int> nonzero_j;
	vector<int> nonzero_l;	// ancestral state indexes
	
	/* Setup the vectors for row, column probs */
	vector<float> nonzero_probs;
	
	/* Store the states_indices as a vector of vectors */
	/* http://stackoverflow.com/questions/8159872/r-list-of-numeric-vectors-c-2d-array-with-rcpp */
	
	/* Empty vector of vectors for the states */
	/* May be the memory problem; fix by creating an internal int vector and inserting it */
	vector< vector<int> > states_vecs;
	//vector<int> tempvec;
	//tempvec.push_back(tmp);
	
	for(int i=0; i<numstates; i++)
		{
		/* Get the subvector */
		vector<int> tmp = Rcpp::as<vector<int> > (xl[i]);
		
		/* Push to the original list */
		states_vecs.push_back(tmp);
		}
	
	/* Make an int array of the range sizes of each state */
	
	/* This makes an vector of pointers, or something... (can't get returned to R by return) */
	/* int rangesizes [numstates]; */
	
	// Make a list of start indices for each range size
	vector<int> range_size_category_indexes;
	int oldsize=0;
	int newsize=0;

	/* Define a vector of size rangesizes */
	Rcpp::NumericVector rangesizes (numstates);
	/* Rcpp::NumericVector rangesizes( Rcpp::Vector<RTYPE>::Vector(const int&) tmpvec ); */
	for (int r=0; r<numstates; r++)
		{
		/* Rcpp::NumericVector tmp_areas( Rcpp::as<SEXP>(xl[r]) ); */
		/* Rcpp::NumericVector tmp_area = Rcpp::as<Rcpp::NumericVector> (xx[r]); */
		newsize = states_vecs[r].size();
		rangesizes[r] = newsize;
		
		// If the new size is bigger than the previous state range size, record this
		if (newsize > oldsize)
			{
			range_size_category_indexes.push_back(r);
			oldsize = newsize;
			}
		/* rangesizes[r] = 1; */
		}
	
	/* The length of the output probabilities */
	Rcpp::NumericVector sp_rowsums(n_xa);
	
	float tmp_probval=0.0;
	
	// The final result
	//int max_min_rangesize = *it5;
	
	// No longer works -- 2013-04-10
	// This prints the memory address
	// std::cout << max_min_rangesize << endl;
	// This prints the value(s) transferred
	// std::cout << max_min_rangesize_c << endl;
	
	//std::cout << "A2" << endl;
	
	//float testval = maxval(maxent01sub);

	/* Identify sympatry (range-copying) */
	/* Identify sympatry (range-copying) */
	/* Identify sympatry (range-copying) */
	//if (printmat >= 1) {cout << "Start sympatry (range-copying)" << "\n";};
	if (std::abs(xy) > min_precision)
		{
		/* l: go through each row of the ancestral states */

		/* Get the subvector of nonzero descendant for left and right branches (i and j) */
		/* For range-copying, we easily know the maximum size of the vector (= number of ancestral states) */
		vector<int> tmp_inums(numstates);
		vector<int> tmp_jnums(numstates);
		vector<int> tmp_lnums(numstates);

		/* Get the subvector of nonzero column probabilities */
		vector<float> tmp_colprobs(numstates);
		
		/* count the number of states we put in; this is the actual size of the output vector */
		int filled_counter = 0;
		filled_counter = filled_counter + 0;	// just so it gets used

		
		// Go through each ancestral state
		for (int l = 0; l < numstates; l++)
			{

			/* This is the ancestral state (the areas) for this row */
			vector<int> ancestral_areas = states_vecs[l];
			
			/* Get the ancestral range size */
			//int ancsize = 0;
			int ancsize = states_vecs[l].size();
			
			// If the size is allowed, then add sympatric range-copying!
			if (ancsize <= max_min_rangesize_c[(ancsize-1)])
				{
				// sympatric prob (xy) * prob as a function of size * relprob leftright combination
				tmp_probval = xy * maxent01symp((ancsize-1), (ancsize-1)) * 1 * 1;
				
				// Add the i index (left side state) and j index (right side state) to the list for this ancestral row
				if (std::abs(tmp_probval) > min_precision)
					{
					tmp_inums[filled_counter] = l;
					tmp_jnums[filled_counter] = l;
					tmp_lnums[filled_counter] = l;
					tmp_colprobs[filled_counter] = tmp_probval;
					filled_counter++;	// increment the filled counter
					}
				else
					{
					continue; // go to next ancestral state in the list
					}
				
				}
			else
				{
				continue; // go to next ancestral state in the list
				}
			} // End iteration through ancestral states
		
		// Now, add the sympatric speciation events to the main COO speciation lists
		// in efficient fashion
		// using vector1.insert():
		// http://stackoverflow.com/questions/201718/concatenating-two-stl-vectors
		nonzero_i.insert(nonzero_i.end(), tmp_inums.begin(), tmp_inums.begin()+filled_counter);
		nonzero_j.insert(nonzero_j.end(), tmp_jnums.begin(), tmp_jnums.begin()+filled_counter);
		nonzero_l.insert(nonzero_l.end(), tmp_lnums.begin(), tmp_lnums.begin()+filled_counter);
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs.begin(), tmp_colprobs.begin()+filled_counter);
		
		// In this case, we don't have to do the reverse
		}
	
	//std::cout << "A3" << endl;

	/* Identify sympatry (range-subset) */
	/* Identify sympatry (range-subset) */
	/* Identify sympatry (range-subset) */
	//if (printmat >= 1) {cout << "Start sympatry (range-subset)" << "\n";};
	if (std::abs(xs) > min_precision)
		{			
		/* l: go through each row of the ancestral states */

		/* Get the subvector of nonzero descendant for left and right branches (i and j) */
		/* For range-copying, we easily know the maximum size of the vector (= number of ancestral states) */
		//if (printmat >= 1) {cout << "checkpoint 0a" << "\n";};

		vector<int> tmp_inums2;
		//if (printmat >= 1) {cout << "checkpoint 0b" << "\n";};

		vector<int> tmp_jnums2;
		vector<int> tmp_lnums2;

		/* Get the subvector of nonzero column probabilities */
		vector<float> tmp_colprobs2;
		//if (printmat >= 1) {cout << "checkpoint 0c" << "\n";};

		/* count the number of states we put in; this is the actual size of the output vector */
		int filled_counter = 0;
		filled_counter = filled_counter + 0;	// just so it gets used
		
		//if (printmat >= 1) {cout << "checkpoint 1" << "\n";};

		// Go through each ancestral state of range size 2 or greater
		for (int l = 0; l < numstates; l++)
			{

			/* This is the ancestral state (the areas) for this row */
			vector<int> ancestral_areas = states_vecs[l];
			
			/* Get the ancestral range size */
			//int ancsize = 0;
			int ancsize = states_vecs[l].size();
			
			// Skip if the ancsize is 1
			if (ancsize < 2)
				{
				continue;
				}
			
			// For every range size between 2 and the maximum minrange, list every combination
			// WITHIN the anc areas--the tough thing here is that you don't have the indexes

			// Go through the states of the possible smaller descendants
			//if (printmat >= 1) {cout << "checkpoint 2" << "\n";};
			// NOTE: This for-loop needs "int", not "unsigned int", because only .size() functions return unsigned int
			for (int desc_size=1; desc_size<=max_min_rangesize_c[(ancsize-1)]; desc_size++)
				{
				//if (printmat >= 1) {cout << desc_size << "	" << max_min_rangesize_c[(ancsize-1)] << "	" << ancsize << "	" << ancsize-1 << "	" << range_size_category_indexes.size() << "\n";};
				
				// Error check so you don't overflow at e.g. 3 areas in anc and desc
				// (shouldn't happen, but still)
				int tmpsize = (int) range_size_category_indexes.size();		// recast the .size() unsigned int as int
				if (desc_size >= tmpsize)
					{
					continue;
					}
				
				int start_state_index = range_size_category_indexes[(desc_size-1)];
				int end_state_index = range_size_category_indexes[desc_size] - 1;
				//if (printmat >= 1) {cout << "checkpoint 3" << "\n";};
				
				//if (printmat >= 1) {cout << "checkpoint 3a	" << start_state_index << "	" << end_state_index << "\n";};

				// Go through the appropriate state indexes
				for (int state_index=start_state_index; state_index<=end_state_index; state_index++)
					{
					//if (printmat >= 1) {cout << "checkpoint 3a	" << start_state_index << "	" << end_state_index << "	" << state_index << "\n";};
					/* areas on the right branch; and size */
					vector<int> rstate_areas( states_vecs[state_index] );
					int rsize = rstate_areas.size();
					//if (printmat >= 1) {cout << "checkpoint 3z" << "\n";};

					// Include only if ALL the descendant right states are in the ancestor
					if (all_ints_found(rstate_areas, ancestral_areas) == true)
						{
						// Also, rule out exact sympatry (y, range-copying)
						if (all_ints_found(ancestral_areas, rstate_areas) == true)
							{
							continue;
							}

						// subset speciation on right branch
						// subset prob (xs) * prob as a function of size * relprob leftright combination
						tmp_probval = xs * maxent01sub((ancsize-1), (rsize-1)) * 1 * 1;
						//row_probval += tmp_probval;

						// Add the i index (left side state) and j index (right side state) to the list for
						// this ancestral row
						if (std::abs(tmp_probval) > min_precision)
							{
							//if (printmat >= 1) {cout << "checkpoint 3b" << "\n";};

							tmp_lnums2.push_back(l);			// ancestral state index
							tmp_inums2.push_back(l);			// left state index
							tmp_jnums2.push_back(state_index);	// right state index
							tmp_colprobs2.push_back(tmp_probval);// prob val
							}
						else
							{
							continue;
							}
						}
					else
						{
						continue;
						}
					}	// End iteration through states of a certain size
				}	// End iteration through desc sizes
			} // End iteration through ancestral states
		//if (printmat >= 1) {cout << "checkpoint 4" << "\n";};
		
		// Now, add the sympatric speciation events to the main COO speciation lists
		// in efficient fashion
		// using vector1.insert():
		// http://stackoverflow.com/questions/201718/concatenating-two-stl-vectors
		nonzero_i.insert(nonzero_i.end(), tmp_inums2.begin(), tmp_inums2.end());
		nonzero_j.insert(nonzero_j.end(), tmp_jnums2.begin(), tmp_jnums2.end());
		nonzero_l.insert(nonzero_l.end(), tmp_lnums2.begin(), tmp_lnums2.end());
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs2.begin(), tmp_colprobs2.end());
		
		// In this case, we DO have to do the reverse
		nonzero_i.insert(nonzero_i.end(), tmp_jnums2.begin(), tmp_jnums2.end());
		nonzero_j.insert(nonzero_j.end(), tmp_inums2.begin(), tmp_inums2.end());
		nonzero_l.insert(nonzero_l.end(), tmp_lnums2.begin(), tmp_lnums2.end());
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs2.begin(), tmp_colprobs2.end());
		//if (printmat >= 1) {cout << "checkpoint 5" << "\n";};
		
		}


	//std::cout << "A4" << endl;


	/* Identify vicariance (range-splitting) */
	/* Identify vicariance (range-splitting) */
	/* Identify vicariance (range-splitting) */
	//if (printmat >= 1) {cout << "Start vicariance (range-splitting)" << "\n";};
	if (std::abs(xv) > min_precision)
		{			
		/* l: go through each row of the ancestral states */

		/* Get the subvector of nonzero descendant for left and right branches (i and j) */
		/* For range-copying, we easily know the maximum size of the vector (= number of ancestral states) */
		vector<int> tmp_inums;
		vector<int> tmp_jnums;
		vector<int> tmp_lnums;

		/* Get the subvector of nonzero column probabilities */
		vector<float> tmp_colprobs;

		vector<int> tmp_inums_LRsame_size;
		vector<int> tmp_jnums_LRsame_size;
		vector<int> tmp_lnums_LRsame_size;

		/* Get the subvector of nonzero column probabilities */
		vector<float> tmp_colprobs_LRsame_size;

		
		/* count the number of states we put in; this is the actual size of the output vector */
		int filled_counter = 0;
		filled_counter = filled_counter + 0;	// just so it gets used

		//std::cout << "A4.1" << endl;
		
		// Go through each ancestral state of range size 2 or greater
		for (int l = 0; l < numstates; l++)
			{

			/* This is the ancestral state (the areas) for this row */
			vector<int> ancestral_areas = states_vecs[l];
			
			/* Get the ancestral range size */
			//int ancsize = 0;
			int ancsize = states_vecs[l].size();
			
			// Skip if the ancsize is 1
			if (ancsize < 2)
				{
				continue;
				}
			
			// For every range size between 2 and the maximum minrange, list every combination
			// WITHIN the anc areas--the tough thing here is that you don't have the indexes

			//std::cout << "A4.2" << endl;

			// Go through the states of the possible smaller descendants
			// NOTE: This for-loop needs "int", not "unsigned int", because only .size() functions return unsigned int
			for (int desc_size=1; desc_size<=max_min_rangesize_c[(ancsize-1)]; desc_size++)
				{
				// Error check so you don't overflow at e.g. 3 areas in anc and desc
				// (shouldn't happen, but still)
				int tmpsize = (int) range_size_category_indexes.size();		// recast the .size() unsigned int as int
				if (desc_size >= tmpsize)
					{
					continue;
					}
				
				int start_state_index = range_size_category_indexes[(desc_size-1)];
				int end_state_index = range_size_category_indexes[desc_size] - 1;
				
				
				//std::cout << "A4.3" << endl;

				
				// Go through the appropriate state indexes
				// NOTE: This for-loop needs "int", not "unsigned int", because only .size() functions return unsigned int
				for (int state_index=start_state_index; state_index <= end_state_index; state_index++)
					{
					/* areas on the right branch; and size */
					vector<int> rstate_areas( states_vecs[state_index] );
					int rsize = rstate_areas.size();
					//std::cout << "A4.4" << endl;

					// Include only if ALL the descendant right states are in the ancestor
					if (all_ints_found(rstate_areas, ancestral_areas) == true)
						{
						//std::cout << "A4.5" << endl;

						// The LEFT states are what's left from ancestral areas after rstate_areas are removed
						vector<int> lstate_areas( ancestral_areas );


						//std::cout << "A4.6" << endl;
						
						// OLD: remove rstates do in reverse order to avoid fuckups
						// find the states to remove
						vector<bool> needs_to_be_removed(lstate_areas.size());
						for (unsigned int larea_index=0; larea_index < lstate_areas.size(); larea_index++)
							{
							// Check if the ancestral area is in rstate_areas; if so, remove
							// vector contains ONE int with value lstate_areas[larea_index]
							vector <int> tmp_one_lstate_area_vector(1, lstate_areas[larea_index]);
							if (all_ints_found(tmp_one_lstate_area_vector, rstate_areas) == true)
								{
								// Edit the list of true/false
								needs_to_be_removed[larea_index] = true;
								} else {
								needs_to_be_removed[larea_index] = false;
								}
							}
							
							
						// Doesn't work:
						//lstate_areas.erase(larea_index);
						
						// works: 	
						// http://stackoverflow.com/questions/4115279/
						// most-efficient-way-of-erasing-deleting-multiple-stdvector-elements-while-retai
						/*
						int last = 0;
						for(int i=0; i<lstate_areas.size(); ++i, ++last)
							{
							while(needs_to_be_removed[i])
								{
								++i;
								}
							if (i >= lstate_areas.size())
								{
								break;
								}
							lstate_areas[last] = lstate_areas[i];   
							}
						lstate_areas.resize(last);
						// Hopefully worketh...
						*/
						//std::cout << "A4.7" << endl;
						//std::cout << "needs_to_be_removed" << endl;
						//printBoolVec(needs_to_be_removed);
						//std::cout << "needs_to_be_removed.size()" << endl;
						//std::cout << needs_to_be_removed.size() << endl;
						

						// Actually, just make a new fucking list
						vector <int> new_lstate_areas;
						
						// We get a memory error here, probably because the 
						// list rarely doesn't get initialized
						// 2013-04-10 -- FIXED & on CRAN
						bool has_new_lstate_areas_been_initialized = false;
						
						for (unsigned int i=0; i<needs_to_be_removed.size(); i++)
							{
							if (needs_to_be_removed[i] == false)
								{
								new_lstate_areas.push_back(lstate_areas[i]);
								has_new_lstate_areas_been_initialized = true;
								}
							}
						
						// If it never gets initialized, then there are NO new_lstate_areas, and we should just continue without adding anything
						// This does fix the bug
						if (has_new_lstate_areas_been_initialized == false)
							{
							continue;
							}
						
						// Now you have just the remaining, true, corresponding lstate_areas
						int lsize = new_lstate_areas.size();
						//std::cout << "lsize: " << lsize << endl;
						//std::cout << "A4.8" << endl;
						//std::cout << "printvec new_lstate_areas" << endl;
						//printvec(new_lstate_areas);
						//std::cout << "end new_lstate_areas" << endl;
						
						// Search for the appropriate index inside this size class
						int start_of_possible_lindexes_to_search = range_size_category_indexes[(lsize-1)];
						int end_of_possible_lindexes_to_search = range_size_category_indexes[lsize] - 1;
						for (int possible_lindex=start_of_possible_lindexes_to_search; possible_lindex <= end_of_possible_lindexes_to_search; possible_lindex++)
							{
							
							// Fixing vicariance calculation; the solution was don't double-count
							// vicariance events where leftsize = rightsize
							/*
							//cout << "\nleft: "; 
							for (int x=0; x<new_lstate_areas.size(); x++) {cout << new_lstate_areas[x] << " ";}
							//cout << "right: "; 
							for (int x=0; x<rstate_areas.size(); x++) {cout << rstate_areas[x] << " ";}
							//cout << "anc: "; 
							for (int x=0; x<ancestral_areas.size(); x++) {cout << ancestral_areas[x] << " ";}
							//cout << "poss.match: "; 
							for (int x=0; x<states_vecs[possible_lindex].size(); x++) {cout << states_vecs[possible_lindex][x] ;}
							//cout << " index:" << possible_lindex;
							*/
							//std::cout << "A4.9" << endl;
							
							//printvec(new_lstate_areas);
							//std::cout << possible_lindex << endl;
							//printvec(states_vecs[possible_lindex]);
							//std::cout << xv << endl;
							//std::cout << ancsize << endl;
							//std::cout << rsize << endl;
							//std::cout << maxent01vic((ancsize-1), (rsize-1)) << endl;
							//std::cout << all_ints_equal(new_lstate_areas, states_vecs[possible_lindex]) << endl;
							
							if (all_ints_equal(new_lstate_areas, states_vecs[possible_lindex]))
								{
								//cout << " MATCH\n" << possible_lindex;
								
								// vicariance speciation on right branch
								tmp_probval = xv * maxent01vic((ancsize-1), (rsize-1)) * 1 * 1;

								//std::cout << "A4.10" << endl;


								// Add the i index (left side state) and j index (right side state) to the list for
								// this ancestral row
								if (std::abs(tmp_probval) > min_precision)
									{
									
									// 1. Push back both the identified state AND the converse, WHEN the two
									// descendant range sizes are DIFFERENT
									if (new_lstate_areas.size() != rstate_areas.size())
										{
										tmp_lnums.push_back(l);			// ancestral state index
										tmp_inums.push_back(possible_lindex);	// left state index
										tmp_jnums.push_back(state_index);		// right state index
										tmp_colprobs.push_back(tmp_probval);// prob val
										} else {
										// 2. However, only do it ONCE when the range sizes are the SAME, since the 
										// algorithm above gets those independently
										tmp_lnums_LRsame_size.push_back(l);			// ancestral state index
										tmp_inums_LRsame_size.push_back(possible_lindex);	// left state index
										tmp_jnums_LRsame_size.push_back(state_index);		// right state index
										tmp_colprobs_LRsame_size.push_back(tmp_probval);// prob val
										}
									
									}
								// exit loop
								break;
								}
							else
								{
								continue;
								}
							} // End iteration of left indexes (larger range sizes)
						} // End if statement for an appropriate right desc under vicariance
					} // End iteration through right states of a certain size
				}	// End iteration through right desc sizes
			} // End iteration through ancestral states
		
		//std::cout << "A5" << endl;

		
		// 2. However, only do it ONCE when the range sizes are the SAME, since the 
		// algorithm above gets those independently
		nonzero_i.insert(nonzero_i.end(), tmp_inums_LRsame_size.begin(), tmp_inums_LRsame_size.end());
		nonzero_j.insert(nonzero_j.end(), tmp_jnums_LRsame_size.begin(), tmp_jnums_LRsame_size.end());
		nonzero_l.insert(nonzero_l.end(), tmp_lnums_LRsame_size.begin(), tmp_lnums_LRsame_size.end());
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs_LRsame_size.begin(), tmp_colprobs_LRsame_size.end());
		
		
		// 1. Push back both the identified state AND the converse, WHEN the two
		// descendant range sizes are DIFFERENT
		// Now, add the vicariance speciation events to the main COO speciation lists
		// in efficient fashion
		// using vector1.insert():
		// http://stackoverflow.com/questions/201718/concatenating-two-stl-vectors
		nonzero_i.insert(nonzero_i.end(), tmp_inums.begin(), tmp_inums.end());
		nonzero_j.insert(nonzero_j.end(), tmp_jnums.begin(), tmp_jnums.end());
		nonzero_l.insert(nonzero_l.end(), tmp_lnums.begin(), tmp_lnums.end());
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs.begin(), tmp_colprobs.end());
		
		// In this case, we DO have to do the reverse
		nonzero_i.insert(nonzero_i.end(), tmp_jnums.begin(), tmp_jnums.end());
		nonzero_j.insert(nonzero_j.end(), tmp_inums.begin(), tmp_inums.end());
		nonzero_l.insert(nonzero_l.end(), tmp_lnums.begin(), tmp_lnums.end());
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs.begin(), tmp_colprobs.end());
		
		}
	//std::cout << "A6" << endl;

	/* Identify jump (new ranges) */
	/* Identify jump (new ranges) */
	/* Identify jump (new ranges) */
	//if (printmat >= 1) {cout << "Start jump (new ranges)" << "\n";};
	if (std::abs(xj) > min_precision)
		{
		/* l: go through each row of the ancestral states */

		/* Get the subvector of nonzero descendant for left and right branches (i and j) */
		/* For range-copying, we easily know the maximum size of the vector (= number of ancestral states) */
		vector<int> tmp_inums;
		vector<int> tmp_jnums;
		vector<int> tmp_lnums;

		/* Get the subvector of nonzero column probabilities */
		vector<float> tmp_colprobs;
		
		/* count the number of states we put in; this is the actual size of the output vector */
		int filled_counter = 0;
		filled_counter = filled_counter + 0;	// just so it gets used
		
		// Go through each ancestral state of range size 2 or greater
		for (int l = 0; l < numstates; l++)
			{

			/* This is the ancestral state (the areas) for this row */
			vector<int> ancestral_areas = states_vecs[l];
			
			/* Get the ancestral range size */
			//int ancsize = 0;
			int ancsize = states_vecs[l].size();
			
			// Skip if the ancsize is ALL AREAS
			// CHECK
			if (ancsize > 20)
				{
				continue;
				}
			
			// For every range size between 1 and (1 less than all areas), list every combination
			// OUTSIDE the anc areas--the tough thing here is that you don't have the indexes

			// Go through the states of the possible smaller descendants
			// NOTE: This for-loop needs "int", not "unsigned int", because only .size() functions return unsigned int
			for (int desc_size=1; desc_size<=max_min_rangesize_c[(ancsize-1)]; desc_size++)
				{
				// Error check so you don't overflow at e.g. 3 areas in anc and desc
				// (shouldn't happen, but still)
				int tmpsize = (int) range_size_category_indexes.size();		// recast the .size() unsigned int as int
				if (desc_size >= tmpsize)
					{
					continue;
					}
				
				int start_state_index = range_size_category_indexes[(desc_size-1)];
				int end_state_index = range_size_category_indexes[desc_size] - 1;
				
				// Go through the appropriate state indexes
				for (int state_index=start_state_index; state_index<=end_state_index; state_index++)
					{
					/* areas on the right branch; and size */
					vector<int> rstate_areas( states_vecs[state_index] );
					int rsize = rstate_areas.size();

					// Include only if ALL the descendant right states are OUTSIDE the ancestor,
					// i.e., NONE inside the ancestor
					if (any_ints_equal(rstate_areas, ancestral_areas) == false)
						{
						// jump speciation on right branch

						// jump/founder event speciation on right branch
						// jump prob (xj) * prob as a function of size * relprob leftright combination
						// * dmatc((ancsize-1), (rsize-1)) = jump dispersal as a function of distance
						int try_jump_dispersal_based_on_dist = 1;
						float jprob_for_cell_based_on_distances = 0.0;
						if (try_jump_dispersal_based_on_dist == 1)
							{
							for (unsigned int ancarea_i=0; ancarea_i<ancestral_areas.size(); ancarea_i++)
								{
								for (unsigned int decarea_i=0; decarea_i<rstate_areas.size(); decarea_i++)
									{
									jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances +
dmatc(ancestral_areas[ancarea_i], rstate_areas[decarea_i]);
									}
								}
							// Normalize (divide by the number of possible jump events) 
							// so that j parameter can be bigger
							float flt_aasize(ancestral_areas.size());
							float flt_rsize(rstate_areas.size());
							jprob_for_cell_based_on_distances = jprob_for_cell_based_on_distances / ( flt_aasize * flt_rsize);
							}
						else
							{
							jprob_for_cell_based_on_distances = 1.0;
							}


						tmp_probval = xj * maxent01jump((ancsize-1), (rsize-1)) * 1 * 1 *
jprob_for_cell_based_on_distances;
						//row_probval += tmp_probval;

						// Add the i index (left side state) and j index (right side state) to the list for
						// this ancestral row
						if (std::abs(tmp_probval) > min_precision)
							{
							tmp_lnums.push_back(l);			// ancestral state index
							tmp_inums.push_back(l);			// left state index
							tmp_jnums.push_back(state_index);	// right state index
							tmp_colprobs.push_back(tmp_probval);// prob val
							}
						else
							{
							continue;
							}
						}
					else
						{
						continue;
						}
					}	// End iteration through states of a certain size
				}	// End iteration through desc sizes
			} // End iteration through ancestral states
	
		//std::cout << "A7" << endl;

	
		// Now, add the sympatric speciation events to the main COO speciation lists
		// in efficient fashion
		// using vector1.insert():
		// http://stackoverflow.com/questions/201718/concatenating-two-stl-vectors
		nonzero_i.insert(nonzero_i.end(), tmp_inums.begin(), tmp_inums.end());
		nonzero_j.insert(nonzero_j.end(), tmp_jnums.begin(), tmp_jnums.end());
		nonzero_l.insert(nonzero_l.end(), tmp_lnums.begin(), tmp_lnums.end());
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs.begin(), tmp_colprobs.end());
		
		// In this case, we DO have to do the reverse
		nonzero_i.insert(nonzero_i.end(), tmp_jnums.begin(), tmp_jnums.end());
		nonzero_j.insert(nonzero_j.end(), tmp_inums.begin(), tmp_inums.end());
		nonzero_l.insert(nonzero_l.end(), tmp_lnums.begin(), tmp_lnums.end());
		nonzero_probs.insert(nonzero_probs.end(), tmp_colprobs.begin(), tmp_colprobs.end());
			
		//std::cout << "A8" << endl;
		
		
		}


	/* Rcpp::NumericVector output(Rcpp::as<SEXP>(rangesizes)); */
	
	/* return the probability vector */
	//if (printmat == 1) {cout << "END cpp_calc_anclikes_sp_COOweights_faster()\n";};
	//if (printmat == 1) {cout << "\n";};
	//std::cout << "A9" << endl;


	
	//return sp_rowsums;
	return(Rcpp::List::create(nonzero_l, nonzero_i, nonzero_j, nonzero_probs));
	
	
	/* Return a single integer */
	/* return Rcpp::wrap(numstates); */
	//return Rcpp::wrap(*it);
	//return Rcpp::wrap(testval);
	}
} // end extern "C"





// Take COO_weights_columnar, which has:
// COO_weights_columnar[[1]] = ancestral indexes
// COO_weights_columnar[[2]] = left indexes
// COO_weights_columnar[[3]] = right indexes
// COO_weights_columnar[[4]] = probability of this split
//
// Return rowsums
extern "C" 
{
SEXP cpp_calc_rowsums_for_COOweights_columnar(SEXP RCOO_weights_columnar_anc_i_list, SEXP RCOO_probs_list, SEXP Rnumstates)
	{
	
	// Cast for Rcpp
	Rcpp::IntegerVector COO_weights_columnar_anc_i_list(RCOO_weights_columnar_anc_i_list);
	Rcpp::NumericVector COO_probs_list(RCOO_probs_list);


	// numstates is the total number of states
	// we will use numstates-1 (no null range)
	int numstates = Rcpp::as<int>(Rnumstates);

	// Initialize a vector to hold the rowsums
	// Fill with (numstates-0) zeros
	//vector <float> rowsums((numstates-0), 0.0); // can't be returned
	Rcpp::NumericVector rowsums((numstates-0), 0.0);
	
	for (int i=0; i<COO_weights_columnar_anc_i_list.size(); i++)
		{
		rowsums[COO_weights_columnar_anc_i_list[i]] = rowsums[COO_weights_columnar_anc_i_list[i]] + COO_probs_list[i];
		}
	
	return(rowsums);
	}
} // end extern "C"





// Take COO_weights_columnar, which has:
// COO_weights_columnar[[1]] = ancestral indexes
// COO_weights_columnar[[2]] = left indexes
// COO_weights_columnar[[3]] = right indexes
// COO_weights_columnar[[4]] = probability of this split
//
// Combine with rowsums and left/right indices to produce actual split likelihoods
extern "C" 
{
SEXP cpp_calc_splitlikes_using_COOweights_columnar(SEXP leftprobs, SEXP rightprobs, SEXP RCOO_weights_columnar_anc_i_list, SEXP RCOO_left_i_list, SEXP RCOO_right_j_list, SEXP RCOO_probs_list, SEXP Rsp_rowsums)
	{
	
	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::NumericVector xa(leftprobs);
	Rcpp::NumericVector xb(rightprobs);
	Rcpp::NumericVector sp_rowsums(Rsp_rowsums);
	
	// Cast for Rcpp
	Rcpp::IntegerVector COO_weights_columnar_anc_i_list(RCOO_weights_columnar_anc_i_list);
	Rcpp::IntegerVector COO_left_i_list(RCOO_left_i_list);
	Rcpp::IntegerVector COO_right_j_list(RCOO_right_j_list);
	Rcpp::NumericVector  COO_probs_list(RCOO_probs_list);

	// numstates is the total number of states
	// we will use numstates-1 (no null range)
	//int numstates = Rcpp::as<int>(Rnumstates);

	// Initialize a vector to hold the rowsums
	// Fill with (numstates-1) zeros
	// (or, same length as sp_rowsums)
	//vector <float> tmp_split_likes(sp_rowsums.size(), 0.0);
	//vector <float> split_likes(sp_rowsums.size(), 0.0);
	Rcpp::NumericVector tmp_split_likes(sp_rowsums.size(), 0.0);
	Rcpp::NumericVector split_likes(sp_rowsums.size(), 0.0);
	
	for (int i=0; i<COO_weights_columnar_anc_i_list.size(); i++)
		{
		float tmp_probval = xa[COO_left_i_list[i]] * xb[COO_right_j_list[i]] * COO_probs_list[i];
		
		tmp_split_likes[COO_weights_columnar_anc_i_list[i]] = tmp_split_likes[COO_weights_columnar_anc_i_list[i]] + tmp_probval;
		}
	
	
	// Now, divide by the rowsums
	for (int j=0; j<tmp_split_likes.size(); j++)
		{
		float tmp_probval = tmp_split_likes[j] / sp_rowsums[j];
		split_likes[j] = tmp_probval;
		}
	
	/* 
	
	# Behavior of:
	# rcpp_calc_splitlikes_using_COOweights_columnar

	tmpca_2 = c(1, 5e-25, 1e-12)
	tmpcb_2 = c(5e-25, 1, 1e-12)

	COO_weights_columnar = list(c(0L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 1L, 0L, 1L), c(0L, 
	1L, 2L, 2L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L), c(0L, 1L, 0L, 1L, 
	2L, 2L, 0L, 1L, 1L, 0L, 0L, 1L), c(3.33333332491748e-06, 3.33333332491748e-06, 
	3.33333332491748e-06, 3.33333332491748e-06, 3.33333332491748e-06, 
	3.33333332491748e-06, 3.33333332491748e-06, 3.33333332491748e-06, 
	2.99998998641968, 2.99998998641968, 2.99998998641968, 2.99998998641968
	))

	printmat = FALSE

	Rsp_rowsums = c(5.999983, 5.999983, 0.000020)
	#Rsp_rowsums = c(6, 6, 6)

	Rcpp_leftprobs = tmpca_2
	Rcpp_rightprobs = tmpcb_2


	# Print the matrix output to screen?	
	if (printmat == TRUE)
	{
	Rprintmat = 1
	} else {
	Rprintmat = 0
	}

	RCOO_weights_columnar_anc_i_list = COO_weights_columnar[[1]]
	RCOO_left_i_list = COO_weights_columnar[[2]]
	RCOO_right_j_list = COO_weights_columnar[[3]]
	RCOO_probs_list = COO_weights_columnar[[4]]

	# Call the fast C++ function
	splitlikes = .Call( "cpp_calc_splitlikes_using_COOweights_columnar", leftprobs=as.numeric(Rcpp_leftprobs), rightprobs=as.numeric(Rcpp_rightprobs), RCOO_weights_columnar_anc_i_list=as.integer(RCOO_weights_columnar_anc_i_list), RCOO_left_i_list=as.integer(RCOO_left_i_list), RCOO_right_j_list=as.integer(RCOO_right_j_list), RCOO_probs_list=as.numeric(RCOO_probs_list), Rsp_rowsums=as.numeric(Rsp_rowsums), PACKAGE = "cladoRcpp" )
	splitlikes


	# R version of the above calculations:
	# 
	num_splits_scenarios = length(COO_weights_columnar[[1]])
	num_ancestral_ranges = length(Rsp_rowsums)
	tmp_split_likes = rep(0, times=num_ancestral_ranges)
	cat("\n")
	for (i in 1:num_splits_scenarios)
	{
	tmp_probval = Rcpp_leftprobs[1+COO_weights_columnar[[2]][i]] * Rcpp_rightprobs[1+COO_weights_columnar[[3]][i]] * COO_weights_columnar[[4]][i]
	tmp_probval
	cat(tmp_probval,"\n")
	tmp_split_likes[1+COO_weights_columnar[[1]][i]] = tmp_split_likes[1+COO_weights_columnar[[1]][i]] + tmp_probval
	}

	split_likes = tmp_split_likes / Rsp_rowsums
	split_likes


	cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], COO_weights_columnar[[4]])
	Rsp_rowsums

	
	*/ 	
	
	return(split_likes);
	}
} // end extern "C"






// Taking the COO_probs_list as an input, and the leftprobs and rightprobs, combine them to produce ancprobs
extern "C" 
{
SEXP cpp_calc_anclikes_sp_using_COOprobs(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP
RCOO_left_i_list, SEXP RCOO_right_j_list, SEXP RCOO_probs_list, SEXP Rsp_rowsums)
	{
	using namespace std;

	//float min_precision = 1e-30; // not used

	/* print if Rprintmat == 1 */
	//int printmat = Rcpp::as<int>(Rprintmat);	// not used
	//if (printmat == 1) {cout << "BEGIN cpp_calc_anclikes_sp_using_COOprobs()" << "\n";};

	/* Define the numeric vectors and put in the data into C++ from R */
	Rcpp::NumericVector xa(leftprobs);
	Rcpp::NumericVector xb(rightprobs);
	Rcpp::NumericVector sp_rowsums(Rsp_rowsums);
	Rcpp::List COO_left_i_list(RCOO_left_i_list);
	Rcpp::List COO_right_j_list(RCOO_right_j_list);
	Rcpp::List COO_probs_list(RCOO_probs_list);

	/* Now, go through the ancestral states (speciation matrix rows), and do the floating-point calculation for each*/
	//std::vector<float> anc_relprobs (sp_rowsums.size(), 0.0);
	Rcpp::NumericVector anc_relprobs(sp_rowsums.size(), 0.0);

	for (int m=0; m<anc_relprobs.size(); m++)
		{
		/* This is the indexes of the left states (i) */
		vector<int> left_state_indices = COO_left_i_list[m];
		vector<int> right_state_indices = COO_right_j_list[m];
		vector<float> this_row_COO_probs_list = COO_probs_list[m];
		
		// Go through the nonzero cells for this particular row/ancestral state and add them up
		float tmp_probval = 0.0;
		for (unsigned int n=0; n<left_state_indices.size(); n++)
			{
			// Take the probability of the left descendant, times the right descendant, times the prob of that
			// split given the ancestor, divided by the sum of the probs for this row when prob of each left and right
			// state is 1 (this is in sp_rowsums)
			float cellval;
			cellval = xa[left_state_indices[n]] * xb[right_state_indices[n]] * this_row_COO_probs_list[n] /
sp_rowsums[m];
			
			// add to tmp_probval
			tmp_probval += cellval;
			}
		anc_relprobs[m] = tmp_probval;
		}

	// End message
	//if (printmat == 1) {cout << "END cpp_calc_anclikes_sp_COOprobs()\n";};
	//if (printmat == 1) {cout << "\n";};
	
	return(anc_relprobs);
	}
} // end extern "C"






