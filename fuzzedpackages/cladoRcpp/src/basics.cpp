#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <vector>

// Actually, you DO want this
// You don't want this is .cpp, it causes double-call errors, e.g.:
// basics.cpp: In function ‘SEXPREC* mult2probvect(SEXPREC*, SEXPREC*)’:
// basics.cpp:59: error: declaration of C function ‘SEXPREC* mult2probvect(SEXPREC*, SEXPREC*)’ conflicts with
// basics.h:26: error: previous declaration ‘SEXPREC* mult2probvect()’ here
// make: *** [basics.o] Error 1
//
// BUT, #include "basics.h" makes maxval & minval findable in cladoRcpp.so
//
#include "basics.h"

// for combn function
#include <math.h>

using namespace std;

// C++ combn function (slow in R)
// from:
// http://www.biostatisticien.eu/textes/rc0408.pdf
void moncombn(int* combmat, int* n, int* m)
	{
	int i, j, e, h, nmmp1, mp1;
	int* a;
	a=new int[*(m+0)];
	for (i=1;i<=*(m+0);i=i+1)
		{
		*(a+i-1)=i;
		}
	e=0;
	h=*(m+0);
	for (i=1;i<=*(m+0);i=i+1)
		{
		*(combmat+i-1)=i;
		}
	i=2;
	
	nmmp1=*(n+0) - *(m+0) + 1;
	mp1=*(m+0) + 1;
	
	while(*(a+0) != nmmp1)
		{
		if(e < *(n+0) - h)
			{
			h=1;
			e=*(a+*(m+0)-1);
			*(a+*(m+0) - h)=e + 1;
			for (j=1;j<=*(m+0);j=j+1)
				{
				*(combmat+(i-1)**(m+0)+j-1)=*(a+j-1);
				}
			i=i+1;
			}
		else
			{
			h=h + 1;
			e=*(a+mp1 - h-1);

			for (j=1;j<=h;j=j+1)
				{
				*(a+*(m+0) - h + j-1)=e + j;
				}
			for (j=1;j<=*(m+0);j=j+1)
				{
				*(combmat+(i-1)**(m+0)+j-1)=*(a+j-1);
				}
			i=i + 1;
			}
		}
	//On libere de la memoire
	delete[] a;
	}



void moncombn_zerostart(int* combmat, int* n, int* m)
	{
	int i, j, e, h, nmmp1, mp1;
	
	// a is a pointer to the beginning of an array of integers of length stored in m
	// declare
	int* a;
	// initialize address
	a=new int[*(m+0)];
	
	// Set the values referenced at address a to 1...m
	for (i=1; i<=*(m+0); i=i+1)
		{
		*(a+i-1)=i;
		}
		
	// h is the value of the 
	e=0;
	
	// h is the value stored at address m
	// i.e., number of rows
	h=*(m+0);
	
	// fill in the first column with e.g. 1,2,3
	for (i=1; i<=*(m+0); i=i+1)
		{
		// -1 at the end to make it zerostart
		*(combmat+i-1)=i-1;
		}
	i=2;
	
	// #cols - #rows + 1
	nmmp1=*(n+0) - *(m+0) + 1;
	mp1=*(m+0) + 1;
	
	while(*(a+0) != nmmp1)
		{
		if(e < *(n+0) - h)
			{
			h=1;
			e=*(a+*(m+0)-1);
			*(a+*(m+0) - h)=e + 1;
			for (j=1;j<=*(m+0);j=j+1)
				{
				// -1 at the end to make it zerostart
				*(combmat+(i-1)**(m+0)+j-1) = *(a+j-1)-1;
				}
			i=i+1;
			}
		else
			{
			h=h + 1;
			e=*(a+mp1 - h-1);

			for (j=1;j<=h;j=j+1)
				{
				*(a+*(m+0) - h + j-1)=e + j;
				}
			for (j=1;j<=*(m+0);j=j+1)
				{
				// -1 at the end to make it zerostart
				*(combmat+(i-1)**(m+0)+j-1)=*(a+j-1)-1;
				}
			i=i + 1;
			}
		}
	//On libere de la memoire
	delete[] a;
	}




// N choose K
// http://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
// # Setting all just to be int; the maxval is typically just 10 million
int nChoosek( int n, int k )
{
	// Error check
    if (k > n) return 0;
    
    // Speedup if k way bigger than n
    if (k * 2 > n) k = n-k;

	// Error checks
    if (k <= 0) return 1;
    if (n <= 0) return 0;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


// print an int vector via cout
void printvec (vector<int> myvector1)
	{
	for (unsigned int i=0; i<myvector1.size(); i++)
		{
		//std::cout << myvector1[i] << " ";
		}
	//std::cout << "\n";
	}

// print an bool vector via cout
void printBoolVec (vector<bool> myvector1)
	{
	for (unsigned int i=0; i<myvector1.size(); i++)
		{
		//std::cout << myvector1[i] << " ";
		}
	//std::cout << "\n";
	}


// Merge two vectors (e.g. to check vicariance)
// 
// http://www.cplusplus.com/reference/algorithm/sort/
vector<int> merge_int_vectors (vector<int> myvector1, vector<int> myvector2)
	{
	// Get total vector size
	int s1 = myvector1.size();
	int s2 = myvector2.size();
	int bigsize = s1 + s2;
	
	// define the output vector, outvec
	vector<int> outvec(bigsize);
	
	// sort the two input vectors
	sort(myvector1.begin(), myvector1.end());
	sort(myvector2.begin(), myvector2.end());
	
	// perform the merge
	merge(myvector1.begin(), myvector1.end(), myvector2.begin(), myvector2.end(), outvec.begin());
	
	return outvec;
	}


// equal algorithm example
// http://www.cplusplus.com/reference/algorithm/equal/
/* bool mypredicate (int i, int j)
	{
	return (i==j);
	}
 */
bool all_ints_equal (vector<int> myvector1, vector<int> myvector2)
	{
	//int myints[] = {20,40,60,80,100};          //   myints: 20 40 60 80 100
	//vector<int>myvector (myints,myints+5);     // myvector: 20 40 60 80 100
	
	// start the beginning of the myvector1 to the end of it,
	// compare to start of myvector2 (automatically goes to the end of my vector2)
	if (equal (myvector1.begin(), myvector1.end(), myvector2.begin()))
		{
		return true;
		}
	else
		{
		return false;
		}
	return 0;
	}


bool any_ints_equal (vector<int> myvector1, vector<int> myvector2)
	{
	//int myints[] = {20,40,60,80,100};          //   myints: 20 40 60 80 100
	//vector<int>myvector (myints,myints+5);     // myvector: 20 40 60 80 100
	
	
	// Just do for loops
	bool are_there_any_matches_TF = false;
	
	// j = loop through the longer (ancestral states) vector
	for (unsigned int j=0; j < myvector2.size(); j++)
		{
		for (unsigned int i=0; i < myvector1.size(); i++)
			{
			if (myvector2[j] == myvector1[i])
				{
				// If a match is found, return true
				are_there_any_matches_TF = true;
				return are_there_any_matches_TF;
				}
			}
		}
	// If none found, return false
	return are_there_any_matches_TF;
	}


bool all_ints_found (vector<int> myvector1, vector<int> myvector2)
	{
	//int myints[] = {20,40,60,80,100};          //   myints: 20 40 60 80 100
	//vector<int>myvector (myints,myints+5);     // myvector: 20 40 60 80 100
	
	
	// Just do for loops
	bool did_all_vec1_match_TF = true;

	// i = loop through the shorter (descendant states) vector
	for (unsigned int i=0; i < myvector1.size(); i++)
		{
		// j = loop through the longer (ancestral states) vector
		bool tmp_match_found_TF = false;
		
		for (unsigned int j=0; j < myvector2.size(); j++)
			{
			if (myvector1[i] == myvector2[j])
				{
				// If a match is found, return true
				tmp_match_found_TF = true;
				continue;
				}
			}
		
		// If this is still false, then you're done
		if (tmp_match_found_TF == false)
			{
			did_all_vec1_match_TF = false;
			return did_all_vec1_match_TF;
			}
		}
	// If none found, return false
	return did_all_vec1_match_TF;
	}


// find the *first* int in the first list that is missing in the secon
// should be combined with all_ints_found TF test to avoid bugs
int get_missing_int (vector<int> myvector1, vector<int> myvector2)
	{
	for (unsigned int i=0; i<myvector1.size(); i++)
		{
		int tmp_int_to_find = myvector1[i];
		
		// Set the flag to false; when you find it, it will be set to
		// true.  If this never happens, this int is missing.
		bool tmp_int_is_in_list = false;
		for (unsigned int j=0; j<myvector2.size(); j++)
			{
			int tmp_int_to_that_might_match = myvector2[j];
			
			// If the item is found, update the flag and 
			// skip the rest of this loop
			if (tmp_int_to_find == tmp_int_to_that_might_match)
				{
				tmp_int_is_in_list = true;
				continue;
				}
			}
		// Now check if is is still false;
		if (tmp_int_is_in_list == false)
			{
			// This is the first element that is not found
			return tmp_int_to_find;
			}
		// Otherwise, continue the search
		}
	
	// If you get to the end without finding anything missing, return -1
	return -1;
	}

// input an SEXP, return an int
float maxval(SEXP maxent01sub)
	{
	//using namespace Rcpp;
	using namespace std;
	
	Rcpp::NumericVector tmpvec(maxent01sub);
	
	
	list<float> li;
	for(int i=0; i<tmpvec.size(); i++)
		{
		li.push_back(tmpvec[i]);
		}
	list<float>::const_iterator it; // declare an iterator
	it = max_element(li.begin(), li.end());
	float tmpmaxval = *it;
	
	// free the memory
	// can't, it is a struct _List_const_iterator, not a pointer
	//delete it;
	
	return tmpmaxval;
	}


// input an SEXP, return an int
float minval(SEXP minent01sub)
	{
	//using namespace Rcpp;
	using namespace std;
	
	Rcpp::NumericVector tmpvec(minent01sub);
	
	
	list<float> li;
	for(int i=0; i<tmpvec.size(); i++)
		{
		li.push_back(tmpvec[i]);
		}
	list<float>::const_iterator it; // declare an iterator
	it = min_element(li.begin(), li.end());
	float tmpminval = *it;
	
	// free the memory
	// can't, it is a struct _List_const_iterator, not a pointer
	//delete it;
	
	return tmpminval;
	}

/* Multiply two probability vectors by each other */
SEXP mult2probvect(SEXP leftprobs, SEXP rightprobs) {
	
	/* Define the numeric vectors and put in the data from R */
	Rcpp::NumericVector xa(leftprobs);
	Rcpp::NumericVector xb(rightprobs);
	
	/* Get the sizes of the vectors */
	int n_xa = xa.size();
	int n_xb = xb.size();
	
	/* The length of the overlap */
	/* int nab = n_xa + n_xb - 1; */
	/* Initialize the output vector */
	/* Rcpp::NumericVector xab(nab); */

	/* The length of the output vector */
	int nab = n_xa * n_xb;
	int counter = 0;


	/* Initialize the output vector */
	Rcpp::NumericVector xab(nab);
	
	/* Go through leftprobs many times (i), repeat for each rightprobs (j) */
	for (int j = 0; j < n_xb; j++)
		{
		for (int i = 0; i < n_xa; i++)
			{
			xab[counter] += xa[i] * xb[j];
			
			/* increment the counter */
			counter++;
			}
		}
	return xab;
	}


/* Example Rcpp function: convolve 2 arrays */
/* Same result as, in R: */
/* convolve(a, b, conj=TRUE, type="open") */
SEXP convolve3cpp(SEXP a, SEXP b) {
	/* Define the numeric vectors and put in the data from R */
	Rcpp::NumericVector xa(a);
	Rcpp::NumericVector xb(b);
	
	/* Get the sizes of the vectors */
	int n_xa = xa.size(), n_xb = xb.size();
	
	/* The length of the overlap */
	int nab = n_xa + n_xb - 1;

	Rcpp::NumericVector xab(nab);
	for (int i = 0; i < n_xa; i++)
		for (int j = 0; j < n_xb; j++)
			xab[i + j] += xa[i] * xb[j];
	return xab;
	}

