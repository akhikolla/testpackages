//***********************************************************************

//
//
//                              COMMON.H
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC 
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : Out of class functions used by FISPRO, part of library fispro

//**********************************************************************

#ifndef __COMMON_H
#define __COMMON_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include <string>
#include <cstring>
#include <stdlib.h> 
#include <math.h>
#include <stdexcept>
#include <vector>
#include <list>
#include <cstring>
#include <algorithm>
#include <time.h>
using std::ifstream;
using std::vector;
using std::allocator;
using namespace std;

//! Config file parameter : string delimitor
#define STRING_SEP '\''  
//! Config file parameter : begin of a series of numbers relative to a same data item
#define START_NB '['  
//! Config file parameter : end delimitor of a series of numbers relative to a same data item
#define END_NB   ']'  
//! Config file parameter : field delimitor within a series of numbers relative to a same data item
#define SEPARE       ',' 
#define ESPACE       ' ' 

//! Default parameter for the SearchNb function.
#define DELIM_START  1     
//! Default parameter for the SearchNb function.
#define DELIM_END    0    

//!  Missing values handling:  Identification in ascii files
#define CODE_NAN    "NA"   
//!  Upper limit used for random generation (word size: 8 bits)
#define FIS_RAND_MAX ((unsigned) (1 << (sizeof(int)*8 - 1)) - 1) 

//! Format used by some output results files
#define FORMAT_DOUBLE "%12.3f "
//! Default value, used by FISOUT::InitPossibles, to consider that 2 rule conclusions are equal.
#define EQUAL_CONC EPSILON
//! Tolerance value on group center coordinates used by the Kmeans functions 
#define CVG_THRES 1e-10

#define DATAFILEMAXPATH 300

#define INFINI  1e6
//! Used by calculus
#define EPSILON 1e-6
//! Used by constructors
#define EPSILON2 1e-12 

//when kernel does not exist
#define EMPTYVALUE -1.00001010

#define DBL_EQUAL(x, y) (fabs((x) - (y)) < EPSILON)	// equality test
#define DBL_INF(x, y) (((y) - (x)) > EPSILON)		// inf test
#define DBL_SUP(x, y) (((x) - (y)) > EPSILON)		// sup test
#define DBL_INF_EQUAL(x, y) (((x) - (y)) < EPSILON)
#define DBL_SUP_EQUAL(x, y) (((y) - (x)) < EPSILON)

//! To handle error messages
extern char ErrorMsg[];

//! To handle variable names
extern char ** VarNameG;
extern int NbVarG;


//! To allocate temporary buffers
extern char DataFilePath[];

//******************************************************************
//
//
//                 OUT OF CLASS FUNCTIONS (General use)
//
//
//******************************************************************


//! Returns the maximum line length in stream 'f'
int MaxLineSize( ifstream &f );
//! Returns the maximum line length in file file_name
int MaxLineSize( const char *file_name );

//! Return the field separator in data file (space and tab are not allowed) 
//! 'hdr' is true when the names of vars are in the first line.
char ReadSeparator(const char *FileName, int & hdr);

//! Returns the number of numbres separated by 'separator' in the sting pointed by 'buf'. 
//! Used by the SampleFileSize function.
int CntNbs( char *buf, char separator, char begin=0, char end=0 );


//******************************************************************
//
//           Random generation (Park and Miller version)
//! Set the static seed used by the Park and Miller random generator 
int setseed ( long int );
//! The Park and Miller implementation of a random generator 
double randpm(void);  
//! Random generation of floating point numbers (Park and Miller implementation).  
//! Return:  between 0 and 1
double FisRand(void);                           

//! Ramdom generation from a gaussian distribution. 
//! calls uniform random FisRand n times (roll a dice n times)
//! 'sig' the standard deviation.
//! returns a random value from a gaussian distrib with zero mean and sigma sd
//! 12 is a good default value
double GaussDice( double sig, int n=12);
//
//******************************************************************
//
//           Missing value Handling
//
//******************************************************************
//! Generate a NaN by returning sqrt(-1.0)
double FisMknan(void);   
//! An ansi version of isnan:  return true is 'x' is NaN, else false.
int FisIsnan(double x);   
//******************************************************************




//******************************************************************
//
//           Sort functions
//
//******************************************************************

//! Function called by qsort to sort in a descending order
int CmpDbl(const void * a, const void * b);

//! Function called by qsort to sort in an ascending order. 
//! NaN's are considered as the smallest value and are all put at the end of array.
int CmpDblAsc(const void * a, const void * b);


//! Compute basic statistics of the distribution 'T'. 
//! 'Tsize' is the number of elements in the array, 
//! 'nb' is the number of extreme values not to take into account, 
//! 'median', 'mean', 'min' and 'max' are the corresponding computed values 
//! of the distribution, delta is the standard deviation. 
//! If 'estim' is not zero then 'delta' is estimated, denominator equal to 'Tsize' - 2* 'nb' - 1, 
//! else it is computed with denominator equal to 'Tsize' - 2* 'nb'.
void StatArray(double * T, int Tsize, int nb, double & median, double & mean, double & delta, double & min, double & max, int estim = 0);

// plus quartiles
void StatArrayQuart(double * T, int Tsize, int nb, double & firstq, double & median, double & thirdq, double & min, double & max);

//! Determines the number of unique values in array 'T' of size 'n',  
//! 'NPos' is filled with the number of unique values, 'ValPos' is a pointer on
//! the unique values themselves. The corresponding allocation is done within the function.
void InitUniq(double *T, int n, double * & ValPos, int & NPos);

//! Determines the number of unique values in array 'T' of size 'n',  
//! 'nu' is filled with the number of unique values, 'U' is a pointer on
//! the unique values themselves, 'p' a pointer on the number of occurrences 
//! of each unique value in 'T'. 
//! The corresponding allocation of 'p' and 'U' (arrays of size 'nu') 
//! is done within the function. The user must make sure these pointers are free 
//! of any allocation and initialized, for example to NULL.  
//! 'thres' is a tolerance on value equality. 
//! Return: 0 if everything is fine, -1 if the last unique value is a NaN.
int SortUniq(double * T, int Tsize, double ** U, int ** p, int &nu, double thres);


//******************************************************************
//
//           Search functions
//
//******************************************************************

//! Copy into 'chaine' the part of the string of 'source' located between  
//! the  first two occurrences of 'Delim'. 
//! 'Delim' is the string separator (default: STRING_SEP).  
//1 Return:  1 if no 'delim' are found else 0.
int SearchStr(const char * source, char * chaine, char Delim = STRING_SEP);

//! Searches in the 'source' string the numbers separated by the 'separator'  
//! whose default value is SEPARE. This argument allows the use of a different   
//! separator for the data file and the config file for instance.  
//! The search is achieved from 'DelimBeg' to 'DelimEnd'  
//! (default the whole string). 
//! 'val' is a pointer on the array to be filled by the function. It should be 
//! properly allocated. 
//! 'size' is the size of 'val'. It corresponds to the number of expected (or maximum)  values. 
//! Return:  the number of values read, and stored in 'val', or -1 if error.
int SearchNb(const char * source, double * val, int size,
	     char separator = SEPARE, int DelimBeg = DELIM_START, int DelimEnd = DELIM_END);

//! Search the variable names within string source, sep being the field separator. 
//! Stop when size names are found or when the whole string has been examined. 
//! Fill the VarNameG array and modify the NbVarG value. 
//! Return the number of names.
int SearchVarNames(const char * source, int size, char sep);

//! returns the file name of source excluding the path.
int FileNameIndex(const char * source);

//******************************************************************
//
//                      Sample file reading
//
//******************************************************************

//! Read the file whose name is pointed by 'N' and stored the data in a 2-dim array 
//! pointed out by the returned value. The allocation of the array is done by the function. 
//! The expected format of the file is one item by line, each value (of input and eventually output) 
//! separated by the interface separator. 
//! The number of columns and the number of rows are read from the file and returned 
//! through thh references 'NCol' and 'MRows'. 
//! The returned value is NULL if the array is too big to be held in memory.
double ** ReadSampleFile(const char *N,  int &NCol, int &NRow);

//! Function used by ReadSampleFile to determine the data sizes before allocating 
//! the pointer on the array. 
//! sep is the Field Separator and hdr indicates the presence of a header line. 
//! linelen is the maximum length of a line.
void SampleFileSize(const char *FileName, int &NbCol, int &NbRow, int &LineLen, char sep, int hdr);

//! Function used by ReadSampleFile
void ReadItems(const char *FileName, int NbVals, int NbInd, double ** vv, int bufsize, char sep, int hdr);

//! This function could be used to read only one item from the stream 'flot'. 
//! The file format is the one used by the function ReadSampleFile. 
//! 'bufsize' is the maximum size of a line to be read. 
//! 'sep' is the field separator. 
//! 'N' corresponds to the number of colums of the file.
//! 'Ind' is a pointer on the values which will be filled by the fonction. 
//! This pointer has to be previously allocated (size 'N').  
//! Return:  the number of values read or -1 if an error occurred. 
int ReadOneItem(ifstream & flot, int bufsize, char sep, double * Ind, int N);

//! Fills the array 'C' with the column number 'col' of 'T'. 
//! 'C' has to be previously allocated. 
//! 'n' is a number of items in one column. 
//! The argument 'col' cannot be checked within the function. 
void GetColumn(double ** T, int n, int col, double *C);

//******************************************************************
//
//           k-means clustering 
//
//******************************************************************

//! Mono dimensional unsupervised clustering. 
//! 'T' is the input array, 'ni' is the number of data, 
//! 'C' is the array of group centers to be filled by the function, 
//! 'nc' is the number of groups (or centers). 
//! if 'norm' is true the data will be scaled into the unit space. 
//! Return value:  the  number of iterations to reach the convergence 
//! on the position centers. The convergence threshold is define, 
//! within the function, by a constant, CVG_THRES set to 1e-14. 
//! The array of the centers should be previously allocated and properly
//! initialized by possible values. Otherwise, some of them will not be 
//! used by the algorithm. 
int Kmeans(double *T, int ni, double *C, int nc, int norm = false);

//! Function used by the Kmeans algorithm to assign a value 'v' to 
//! the closest group. 'cg' are the center coordinates and 'ng' 
//! the number of groups. The distance is the euclidean one. 
//! Returns the number of the closest group. 
int AssignClas(double v, double *cg, int ng);

//! Allocate and initialize the C array with n values between min and max
void InitCentres(double *&C, int n, double min, double max);

//! This function checks that clusters are not empty. 
//! It can be used afer center calculation with Kmeans. 
//! Eliminates empty clusters (NaN) and modifies the number of groups 'nc' 
//! Returns the number of empty (and removed) clusters. 
int KmeansNE(double *T, int ni, double * C, int & nc);

//! Multidimensional version of the Kmeans. 
//! 'dim' is the work space dimension, common to both 'T' and 'C'. 
int Kmeans(double **T, int ni, double **C, int nc, int dim, int norm = false);

//! Multidimensional version of AssignClas.  
//! 'dim' is the work space dimension, common to both 'v' and 'cg'. 
int AssignClas(double *v, int dim, double **cg, int ng);

//! Multidimensional version of KmeansNE. 
//! 'dim' is the work space dimension, common to both 'T' and 'C'. 
int KmeansNE(double **T, int ni, double **C, int & nc, int dim);

//  Global normalise to use whith sample data
void Normalize ( double ** SampleData, int col , int nbrow , double min, double max);
//  Global unnormalise to use whith sample data
void UnNormalize ( double ** SampleData, int col , int nbrow , double min, double max);
//  Reading template from file for fuzzy input
void ReadTemplate(char * file, double &KW, double &SW);
//  Writing template file for fuzzy input
void WriteTemplate(char * file, double KW, double SW);
#endif

//*************************  COMMON.H    **************************

