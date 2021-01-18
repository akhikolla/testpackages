#include "fastcmh_cpp.h"

/* CODE DEPENDENCIES */
//#include"time_keeping.c"
//#include "./chi2.h"



/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of buffer used to read files
#define NGRID 500 //Number of samples in the grid of tentative corrected significance thresholds, not counting 1 as it is a trivial threshold
#define LOG10_MIN_PVAL -30.0 //Minimum tentative corrected significance threshold = 10^{LOG10_MIN_PVAL}

/* -------------------------------------------- GLOBAL VARIABLES -----------------------------------------------------------*/

const int ASCII_0 = 48;
const int ASCII_1 = 49;
const int ASCII_9 = 57;
const int ASCII_NEWLINE = 10;
const int MAX_ASCII = 256;
const int ASCII_OTHER = 127;
const int ASCII_OTHER_NEWLINE = 126;

bool saveAllPvals = true;
bool showProcessing = false;
bool doFDR = false;
bool useDependenceFDR = false;

std::string timingString = "";
double timeExecution = 0.0;
double timeInitialisation = 0.0;
double timeFileIO= 0.0;
double timeComputeSigThreshold = 0.0;
double timeComputeSigIntervals = 0.0;
//at the moment, not saving memory, so converting to double
// size_t peakMemoryUsageInBytes = 0;
int peakMemoryUsageInBytes = 0;

//rather save the summary as a string
std::string summaryString = "";
// Summary file (now passed back as a parameter) File to output varied information about the process (number of interval processed and so on)



// max hisgoram file (now passed back as parameter) maxcmh_hist_file; // File to output histogram of maximum attainable CMH statistics (only reliable in the testable range)

// Number of observations, N
long long N;
// Number of observations in positive class
long long n;
// Sequence length
long long L;
long long L_max;
// Number of tables
long long K;
// Number of observations per table
long long *Nt;
// Number of observations in positive class per table
long long *nt;
// Number of observations per table
long long *cum_Nt; //Cumulative sum of Nt
// Now some precomputed quantities to save time when computing the CMH test statistic
long long *Nt_nt; // Ni-ni for each of the K tables
long long *hypercorner_bnd; // max(ni,Ni-ni) for each of the K tables
double *gammat; // ni/Ni for each of the K tables
double *gammabint; // (ni/Ni)*(1-ni/Ni) for each of the K tables
// And some others to save time when computing the maximum CMH test statistic in the "top right" hypercorner
double *f_vals, *g_vals, *betas;
double f_sum, g_sum, Tcmh_max_corner_l, Tcmh_max_corner_r, Tcmh_aux_corner;
long long *idx_betas_sorted;


// Current interval length
long long l;
// Number of testable intervals
long long m;
// Target FWER
double alpha;
// Final corrected significance threshold
double delta_opt;

// Grid of logarithmically spaced corrected significance thresholds. The sequence is ordered
// from larger thresholds to smaller thresholds.
double *pgrid;
// Current tentative corrected significance threshold and index of it in the grid
double pth;
int idx_th;
// Step size in the grid
double log10_p_step;


// Vector of class labels
char *Y_tr;
// The original dataset matrix stored as a LxN matrix, where L is the sequence length
// and N the number of observed sequences. This exploits the row-major storage of C matrices.
char **X_tr;
// Another LxN matrix, storing the values of nodes in the layer above.
char **X_par;

// A LxK-dimensional matrix storing the frequency of each interval in each table as they are processed
long long **freq_par;
// A (Ngrid+1)-dimensional vector such that freq_cnt[j] = #intervals with maximum attainable
// CMH statistic in the j-th bucket
long long *freq_cnt;

// Queue of testable intervals in the layer below
long long *testable_queue;
long long testable_queue_front;
long long testable_queue_length;

// Auxiliary variable to keep track of current layer
long long last_tau;

/* PROFILING VARIABLES */
long long n_intervals_processed;
long long n_pvalues_computed;
long long n_significant_intervals;

/* FUNCTION DECLARATIONS */


int get_N_n(char *);
int read_labels_file(char *);
int get_L(char *);
int read_dataset_file(char *, char *);
int get_K(char *);
int read_covariates_file(char *);
int qsort_cmp_betas(const void *, const void *);
inline int bucket_idx(double);
// A extern function related to the Chi2 distribution
extern double Chi2_sf(double,double); //Chi2 Survival function, from Dominik's EasyGWASCore


//variables for time_keeping

/* GLOBAL VARIABLES (TIME SPENT) */
clock_t t_init, t_end;
double time_initialisation=0;
double time_IO=0;
double time_comp_threshold=0;
double time_comp_significant_intervals=0;
clock_t tic,toc;
// timing_file (now passed back as parameter) File to output information about the runtime of the algorithm, defined in the significant_interval_search_* C files


//GLOBAL VECTOR
//save SIGNIFICANT tau,l,pval
std::vector<long long> sigTau;
std::vector<long long> sigL;
std::vector<double> sigPval;

//MASSIVE change - replace "all" with "testable" - otherwise confusing.
//save ALL tau,l,pval
std::vector<long long> allTestableTau;
std::vector<long long> allTestableL;
std::vector<double> allTestablePval;

//histogram - observation and frequency
// std::vector<long long> histObs;
// std::vector<long long> histFreq;
std::vector<int> histObs;
std::vector<int> histFreq;


//save FDR tau,l,pval
std::vector<long long> fdrTau;
std::vector<long long> fdrL;
std::vector<double> fdrPval;





// http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-g
template < typename T > std::string AnotherToString( const T& num )
{
    std::ostringstream ss;
    ss << num ;
    return ss.str() ;
}

// template <typename T>
// string ToString(T val)
// {
//     return std::to_string(val);
// }
//    //previous method
//     stringstream stream;
//     stream << val;
//     return stream.str();

//need to make my own to_string, since it is only part of C++11
//but now using C++11
// string my_to_string(int i){
//     return std::to_string(i);
// }
//    //prvious method
//     stringstream ss;
//     ss << i;
//     string str = ss.str();
//     return(str);





//-----------------------------gamma.h-----------------------------------//

double regularizedLowerIncompleteGamma(double x, double alpha) {
	if(x <= 0.0 || alpha <= 0.0) return 0.0;
	
    double gamma_f = exp(alpha*log(x) - lgamma(alpha) - x);
	if(x < alpha + 1.0) {
		double i = alpha;
		double tmp_sum = 1.0;
		double sum = tmp_sum;
		while(tmp_sum/sum > 1e-10) {
			i++;
			tmp_sum *= x/i;
			sum += tmp_sum;
		}
		return gamma_f*sum/alpha;
	} else { 
        //Solve via evaluting continued fractions
		//Evaluation of continued fraction
		double a=1.0-alpha;
		double b=1+x+a;
		double pa1 = 1.0;
		double pb1 = x;
		double pa2 = x + 1.0;
		double pb2 = b*x;
		double func = pa2/pb2;
		double pa,pb,ratio,tmp;
		double i = 0;

		while(1) {
			i++;
			a++;
			b += 2.0;
			pa = b*pa2-a*i*pa1;
			pb = b*pb2-a*i*pb1;
			if(pb) {
				ratio = pa/pb;
				tmp = fabs((func-ratio));
				if(tmp<=1e-10*ratio) break;
				func = ratio;
			} else tmp=1.0;
			pa1=pa2;
			pb1=pb2;
			pa2=pa;
			pb2=pb;
			if(i>100) break;//Maximum number if iterations
		}
		return 1.0-func*gamma_f;
	}
}


/*
Function to compute a complemented incomplete gamma function based on a continued fraction
*/
double complementedIncompleteGamma(double x, double alpha) {
	if((x <= 0) || ( alpha <= 0))
		return 1.0;
	
    if((x < 1.0) || (x < alpha))
		return 1.0 - regularizedLowerIncompleteGamma(x,alpha);
	
    double gamma_f = exp(alpha*log(x) - lgamma(alpha) - x);
	
	//continued fraction
	double y = 1.0 - alpha;
	double z = 1.0 + x + y;
	double c = 0.0;
	double pkm2 = 1.0;
	double qkm2 = x;
	double pkm1 = x + 1.0;
	double qkm1 = z * x;
	double func = pkm1/qkm1;
	double i = 0;
	double ratio,tmp,pk,qk;

	while(1) {
		i++;
		c += 1.0;
		y += 1.0;
		z += 2.0;
		pk = pkm1 * z  -  pkm2 * y*c;
		qk = qkm1 * z  -  qkm2 * y*c;
		if( qk != 0 ) {
			ratio = pk/qk;
			tmp = fabs( (func - ratio)/ratio );
			if(tmp<=1e-10*ratio) break;
			func = ratio;
		} else {
			tmp = 1.0;
		}
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;
		if( fabs(pk) > 1e32) {
			pkm2 *= 1e-32;
			pkm1 *= 1e-32;
			qkm2 *= 1e-32;
			qkm1 *= 1e-32;
		}
		if(i>100) break; //Max number of iterations
	}
	return func * gamma_f;
}

//-----------------------------chi2.h-----------------------------------//


/*
Computed the survival function of the chi2 distribution function for x and k
*/
inline double Chi2_sf(double x, double k) {
    return complementedIncompleteGamma(0.5*x,0.5*k);
}

/*
Computed the cumulative distribution function of the chi2 distribution function for x and k
*/
inline double Chi2_cdf(double x, double k) {
    if (k==2.0) {
        return 1.0 - exp(-0.5*x);
    } else {
        return regularizedLowerIncompleteGamma(0.5*x,0.5*k);
    }
}

/*
Probability density function of Chi2
*/
inline double Chi2_pdf(double x, double k) {
    if (x<0.0) return 0.0;
    return pow(x,0.5*k-1.0)*exp(-0.5*x)/(pow(2.0,0.5*k)*tgamma(0.5*k));
}


//-----------------------------time_keeping.c-----------------------------------//


// Measure running time
clock_t measureClocks(){
  return clock();
}

// Measure peak memory usage
// size_t measurePeakMemory(){
//   struct rusage t;
//   getrusage(RUSAGE_SELF, &t);
//   return (size_t)t.ru_maxrss;
// }

// Display execution time and memory consumption
void profileCode(){
// 	size_t peak_memory;
    //clear the string
    timingString.clear();
    //timing string starts
    timingString.append("CODE PROFILING\n");

    timeExecution = (t_end - t_init)/ CLOCKS_PER_SEC;
    timingString.append("Total execution time: " +
            AnotherToString(timeExecution) + " (s).\n");

    timeInitialisation = time_initialisation;
    timingString.append("\tInitialisation time: " +
            AnotherToString(timeInitialisation) + " (s).\n");
    
    timeFileIO = time_IO;
    timingString.append("\tFile I/O time: " + 
            AnotherToString(timeFileIO) + " (s).\n");

    timeComputeSigThreshold = time_comp_threshold;
    timingString.append("\tTime to compute corrected significance threshold: " +
            AnotherToString(timeComputeSigThreshold) + " (s).\n");

    timeComputeSigIntervals = time_comp_significant_intervals;
    timingString.append("\tTime to find significant intervals: " +
            AnotherToString(timeComputeSigIntervals) + " (s).\n");


// 	peak_memory = measurePeakMemory();
//     peakMemoryUsageInBytes  = measurePeakMemory();
    timingString.append("\tPeak memory usage: " +
            AnotherToString(peakMemoryUsageInBytes) + " (bytes).\n");
    
}



/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the main variables, call I/O routines and allocate memory
 * Input arguments are self-explanatory
 * */
// 	sis_init(*X_file_arg, *Y_file_arg, *C_file_arg, alpha_arg, L_max_arg);
int sis_init(char *X_filename, char *Y_filename, char *C_filename, double target_fwer, long long l_max){
	long long j; //Loop variable
// 	clock_t tic,toc;//Internal ones, do not overwrite the ones in time_keeping.c
	double log10_p;

	// Compute total number of observations and number of observations in minority class
// 	tic = measureClocks();
    int get_N_n_SUCCESS = EXIT_FAILURE;
    try {
        get_N_n_SUCCESS = get_N_n(Y_filename);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

    int getL_SUCCESS = EXIT_FAILURE;
    try {
        getL_SUCCESS = get_L(X_filename);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


    int getK_SUCCESS = EXIT_FAILURE;
    try {
        getK_SUCCESS = get_K(C_filename);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

// 	toc = measureClocks();
// 	time_IO += toc-tic;

	// Store core constants
	alpha = target_fwer;
	L_max = l_max;

	// Initialise grid of candidate corrected significance thresholds
	pgrid = (double *)malloc((NGRID+1)*sizeof(double));
	for(log10_p=0,log10_p_step=-LOG10_MIN_PVAL/NGRID,j=0; j<=NGRID; log10_p-=log10_p_step, j++) pgrid[j] = pow(10,log10_p);
	// Initialise threshold values
	idx_th = 1; pth = pgrid[idx_th];


	// Allocate space for per table number of observations and number of observations in positive class
	Nt = (long long *)calloc(K, sizeof(long long));
    try{
        if(!Nt){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array Nt\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	nt = (long long *)calloc(K, sizeof(long long));
    try{
        if(!nt){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array nt\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }



	cum_Nt = (long long *)calloc(K+1, sizeof(long long));
    try {
        if(!cum_Nt){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array cum_Nt\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	// And read covariates file, filling in the array Nt
// 	tic = measureClocks();
    int read_covariates_file_SUCCESS = EXIT_FAILURE;
    try {
        read_covariates_file_SUCCESS = read_covariates_file(C_filename);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

// 	toc = measureClocks();
// 	time_IO += toc-tic;

	// Allocate space for class labels
	Y_tr = (char *)malloc(N*sizeof(char));
    try {
        if(!Y_tr){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array Y_tr\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	// And store them in memory from file, also computing nt along the way
// 	tic = measureClocks();
    int read_labels_file_SUCCESS = EXIT_FAILURE;
    try {
        read_labels_file_SUCCESS = read_labels_file(Y_filename);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

// 	toc = measureClocks();
// 	time_IO += toc-tic;

	// Initialise dataset matrix
	X_tr = (char **)malloc(L*sizeof(char *));
    try {
        if(!X_tr){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array X_tr\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	X_tr[0] = (char *)calloc(L*N, sizeof(char));
    try {
        if(!X_tr[0]){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array X_tr[0]\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	for(j=1; j<L; j++) X_tr[j] = X_tr[0] + j*N;
	// Same for parents
	X_par = (char **)malloc(L*sizeof(char *));
    try {
        if(!X_par){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array X_par\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	X_par[0] = (char *)calloc(L*N, sizeof(char));
    try {
        if(!X_par[0]){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array X_par[0]\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	for(j=1; j<L; j++) X_par[j] = X_par[0] + j*N;

// 	tic = measureClocks();
    int read_dataset_file_SUCCESS = EXIT_FAILURE;
    try {
        read_dataset_file_SUCCESS = read_dataset_file(X_filename, X_tr[0]);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }
// 	toc = measureClocks();
// 	time_IO += toc-tic;


	// Allocate memory for several vectors

	// First some small vectors to precompute some common magnitudes used to evaluate the test statistics
	Nt_nt = (long long *)calloc(K, sizeof(long long));
    try {
        if(!Nt_nt){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array Nt_nt\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	hypercorner_bnd = (long long *)calloc(K, sizeof(long long));
    try {
        if(!hypercorner_bnd){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array hypercorner_bnd\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	gammat = (double *)calloc(K, sizeof(double));
    try {
        if(!gammat){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array gammat\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	gammabint = (double *)calloc(K, sizeof(double));
    try {
        if(!gammabint){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array gammabint\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	// Fill them in
	for(j=0; j<K; j++){
		Nt_nt[j] = Nt[j]-nt[j];
		hypercorner_bnd[j] = (nt[j] > Nt_nt[j]) ? nt[j] : Nt_nt[j];
		gammat[j] = ((double)nt[j])/Nt[j];
		gammabint[j] = gammat[j]*(1-gammat[j]);
	}

	// Some other small vectors which will be used in the function isprunable to avoid recomputing things
	f_vals = (double *)calloc(K, sizeof(double));
    try {
        if(!f_vals){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array f_vals\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	g_vals = (double *)calloc(K, sizeof(double));
    try {
        if(!g_vals){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array g_vals\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	betas = (double *)calloc(K, sizeof(double));
    try {
        if(!betas){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array betas\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	idx_betas_sorted = (long long *)calloc(K, sizeof(long long));
    try {
        if(!idx_betas_sorted){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array idx_betas_sorted\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	// Now some larger data structures used during the enumeration procedure
	testable_queue = (long long *)calloc(L, sizeof(long long));
    try {
        if(!testable_queue){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array testable_queue\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	freq_par = (long long **)calloc(L, sizeof(long long *));
    try {
        if(!freq_par){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array freq_par\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	freq_par[0] = (long long *)calloc(L*K, sizeof(long long));
    try {
        if(!freq_par[0]){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array freq_par[0]\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	for(j=1; j<L; j++) freq_par[j] = freq_par[0] + j*K;
	freq_cnt = (long long *)calloc((NGRID+1), sizeof(long long));

    try {
        if(!freq_cnt){
            throw std::runtime_error("Error in function sis_init: couldn't allocate memory for array freq_cnt\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


    if ( (get_N_n_SUCCESS==EXIT_FAILURE) or 
            (getL_SUCCESS==EXIT_FAILURE) or 
            (getK_SUCCESS==EXIT_FAILURE) or 
            (read_labels_file_SUCCESS==EXIT_FAILURE) or 
            (read_dataset_file_SUCCESS==EXIT_FAILURE) or
            (read_covariates_file_SUCCESS==EXIT_FAILURE)){
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;
}

void output_maxcmh_histogram(){
    long long i;
    for(i=0; i<=NGRID; i++) {
        //save the histogram observation and frequency count
        //actually casting long long to int in vector...
        histObs.push_back(i);
        histFreq.push_back(freq_cnt[i]);
    }
}

/* Free all allocated memory and give some output for debugging purposes */
void sis_end(){
	// Execution time and peak memory consumption
	profileCode();

    // Output histogram of max attainable CMH statistics (only reliable in the testable range)
    output_maxcmh_histogram();    

	// Free allocated memory
	free(Nt); free(nt); free(cum_Nt); free(hypercorner_bnd);
	free(Nt_nt); free(gammat); free(gammabint);
	free(f_vals); free(g_vals); free(betas); free(idx_betas_sorted);
	free(pgrid);
	free(Y_tr);
	free(X_tr[0]); free(X_par[0]);
	free(X_tr); free(X_par);
	free(freq_par[0]); free(freq_par); free(freq_cnt);
	free(testable_queue);
}



/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

/* Computes the CMH p-value as a function of the margins x, n and N and the cell counts a for the K tables */
double compute_pval(long long a, long long *x){
	long long k;
	double num = a, den = 0;
	for(k=0; k<K; k++){
		num -= x[k]*gammat[k];
		den += x[k]*(1-((double)x[k])/Nt[k])*gammabint[k];
	}
	num *= num;
	if(den==0) return 1;
	else return Chi2_sf(num/den,1);

}

/* Computes the minimum attainable CMH p-value depending on the margins x, n and N for the K tables */
double compute_minpval(long long *x){
	long long k;
	double left_tail_num = 0, right_tail_num = 0, den = 0;
	double aux1, aux2;
	for(k=0; k<K; k++){
		aux1 = x[k]-Nt_nt[k]; aux2 = x[k]*gammat[k];
		left_tail_num += ((aux1 > 0) ? aux1 : 0) - aux2;
		right_tail_num += ((x[k] > nt[k]) ? nt[k] : x[k]) - aux2;
		den += x[k]*(1-((double)x[k])/Nt[k])*gammabint[k];
	}
	left_tail_num *= left_tail_num; right_tail_num *= right_tail_num;
	if(den==0) return 1;
	else return Chi2_sf(((left_tail_num > right_tail_num) ? left_tail_num : right_tail_num)/den,1);
}

/* Given the margins of the K tables and the minimum attainable CMH p-value, checks whether the interval
 * and all other intervals containing it can be pruned from the search space
 * */
int isprunable(long long *x){
	long long j,k;
	// If for any of the K tables, its margin x is smaller than the maximum of n and N-n, then we cannot prune
	// the interval (we are not in the "top-right" hypercorner)
	for(k=0; k<K; k++) if(x[k] < hypercorner_bnd[k]) return 0;

	// Compute the maximum value of left handside function
	for(j=0,k=0; k<K; k++){
		// Discard all dimensions for which x[k]==Nt[k], as they don't contribute to the function neither in the
		// numerator nor the denominator
		if(x[k] < Nt[k]){
			f_vals[j] = gammat[k]*(Nt[k]-x[k]);
			g_vals[j] = gammabint[k]*x[k]*(1-((double)x[k])/Nt[k]);
			betas[j] = g_vals[j]/f_vals[j];
			idx_betas_sorted[j] = j;
			j++;
		}
	}
	qsort(idx_betas_sorted,j,sizeof(long long),qsort_cmp_betas); //Equivalent to argsort(betas[0:j])
	f_sum = 0; g_sum = 0; Tcmh_max_corner_l = 0;
	for(k=0; k<j; k++){
		f_sum += f_vals[idx_betas_sorted[k]];
		g_sum += g_vals[idx_betas_sorted[k]];
		Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
		Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? Tcmh_max_corner_l : Tcmh_aux_corner; //Tcmh_max_corner_l=max(Tcmh_max_corner_l,Tcmh_aux_corner)
	}
	// Compute the maximum value of right handside function
	for(j=0,k=0; k<K; k++){
		// Discard all dimensions for which x[k]==Nt[k], as they don't contribute to the function neither in the
		// numerator nor the denominator
		if(x[k] < Nt[k]){
			f_vals[j] = (1-gammat[k])*(Nt[k]-x[k]);
			//g_vals doesn't change, hence it does not need to be recomputed
			betas[j] = g_vals[j]/f_vals[j];
			idx_betas_sorted[j] = j;
			j++;
		}
	}
	qsort(idx_betas_sorted,j,sizeof(long long),qsort_cmp_betas); //Equivalent to argsort(betas[0:j])
	f_sum = 0; g_sum = 0; Tcmh_max_corner_r = 0;
	for(k=0; k<j; k++){
		f_sum += f_vals[idx_betas_sorted[k]];
		g_sum += g_vals[idx_betas_sorted[k]];
		Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
		Tcmh_max_corner_r = (Tcmh_max_corner_r >= Tcmh_aux_corner) ? Tcmh_max_corner_r : Tcmh_aux_corner; //Tcmh_max_corner_r=max(Tcmh_max_corner_r,Tcmh_aux_corner)
	}
	return Chi2_sf(((Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l),1) > pth;
}






//Create separate vectors for 
void process_first_layer_pvalues(){
	long long tau, j, k, queue_idx, a;
	char *X_tr_aux;
	long long *aux_ptr;
	double pval_min_int, pval;
	// Process each length 1 interval
	for(tau=0; tau<L; tau++){
		// Direct pointer to the relevant feature vector
		X_tr_aux = X_tr[tau];
		// Compute number of 1s in the interval for each of the K tables
		for(k=0; k<K; k++){
			aux_ptr = &freq_par[tau][k];
			for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) *aux_ptr += X_tr_aux[j];
		}
		// If the interval is testable...
		// Update frequency-buckets and number of testable intervals
		#ifndef NO_SINGLE_FEATURES
			pval_min_int = compute_minpval(freq_par[tau]);
            //TODO:
            //here we check if pval_min of test is smaller than threshold
			if(pval_min_int <= pth){//Definition of testability in terms of critical values
				// Compute global cell count considering all tables together
				a = 0;
				for(j=0; j<N; j++) if(X_tr_aux[j]) a += Y_tr[j];
				// Compute the p-value
				pval = compute_pval(a,freq_par[tau]); n_pvalues_computed++;
				// Check if the P-value is significant
				if(saveAllPvals) {
                    //TODO: do we really want to save all of these values?
                    //saving ALL values to vector
                    allTestableL.push_back(l+1);
                    //need to add 1 to tau
                    allTestableTau.push_back(tau+1);
                    allTestablePval.push_back(pval);
                }
				if(pval <= delta_opt) { 
                    //saving significant values to vector
                    sigL.push_back(l+1);
                    //need to add 1 to tau
                    sigTau.push_back(tau+1);
                    sigPval.push_back(pval);
                    
                    n_significant_intervals++; 
                }
			}
		#endif
		// If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || isprunable(freq_par[tau]) || isprunable(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}









//Instead of writing to file, saved to vector
void process_intervals_pvalues(){
	long long tau, j, k, queue_idx, a;
	char *X_tr_aux, *X_par_aux;
	long long *aux_ptr;
	double pval_min_int, pval;
	// While testable-interval queue is not empty, continue to process intervals
	while(testable_queue_length){
		// Pop a testable interval from the queue
		tau = testable_queue[testable_queue_front];
		testable_queue_front = (testable_queue_front<(L-1)) ? testable_queue_front + 1: 0;
		testable_queue_length--;
		// Check if we have started processing a new layer by detecting non-monotonicity in tau
		if(tau < last_tau) {
			l++;
			#ifndef NO_VERBOSE
            if (showProcessing){
                Rcpp::Rcout << "\tProcessing layer " << (l+1) << "...\n" << std::endl;
//                 printf("\tProcessing layer %lld...\n",l+1);
            }
			#endif
		}
		if((L_max>0) && ((l+1) > L_max)) {
			#ifndef NO_VERBOSE
            if (showProcessing){
                Rcpp::Rcout << "\tMaximum interval length achieved at l=" << (l+1) << ". Stopping enumeration...\n" << std::endl;
//                 printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n",l+1);
            }
			#endif
			break;
		}
		last_tau = tau;
		// In this case, the testable region does not change, so we don't need to check if the interval
		// has to be pruned now. If it was appended to the queue, then it has to be processed for sure
		// Compute OR and frequency of the interval for each of the K tables
		X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
		for(k=0; k<K; k++){
			aux_ptr = &freq_par[tau][k];
			for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; (*aux_ptr)++;}
		}
		// If the interval is testable, increase counter of testable items and frequency-buckets and
		// check if the corrected significance threshold must be reduced
		pval_min_int = compute_minpval(freq_par[tau]);
		if(pval_min_int <= pth){//Definition of testability in terms of critical values
			// Compute global cell count considering all tables together
			a = 0;
			for(j=0; j<N; j++) if(X_par_aux[j]) a += Y_tr[j];
			// Compute the P-value
			pval = compute_pval(a,freq_par[tau]); n_pvalues_computed++;
			// Check if the P-value is significant
            //
            // TODO - place where sigints are written to file...
			if(saveAllPvals) {
                //saving ALL values to vector
                allTestableL.push_back(l+1);
                //need to add 1 to tau
                allTestableTau.push_back(tau+1);
                allTestablePval.push_back(pval);
            }
			if(pval <= delta_opt) { 
                //saving significant values to vector
                sigL.push_back(l+1);
                //need to add 1 to tau
                sigTau.push_back(tau+1);
                sigPval.push_back(pval);

                n_significant_intervals++; 
            }
		}
		// If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || isprunable(freq_par[tau]) || isprunable(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}

/* Wrapper function that encapsulates the functionality required to find significant intervals */
void find_significant_intervals(){
	// Give feedback to user
	#ifndef NO_VERBOSE
    if (showProcessing){
        Rcpp::Rcout << "\n\nSCANNING DATASET FOR SIGNIFICANT INTERVALS...\n\n" << std::endl;
//         printf("\n\nSCANNING DATASET FOR SIGNIFICANT INTERVALS...\n\n");
    }
	#endif
	// Initialise the queue as empty
	testable_queue_front = 0; testable_queue_length = 0;
	// Initialise current layer index and current number of computed p-values to 0
	l = 0; n_pvalues_computed = 0; n_significant_intervals = 0;
	// Clear the current layer frequency counters
	memset(freq_par[0],0,L*K*sizeof(long long));
	// Initialise the value of the OR vectors of current layer to original dataset
	memcpy(X_par[0],X_tr[0],L*N*sizeof(char));
	// Process the upper-most layer (i.e. the layer composed of length 1 intervals)
	#ifndef NO_VERBOSE
    if (showProcessing){
        Rcpp::Rcout << "\tProcessing layer  " << (l+1) << "...\n" << std::endl;
//         printf("\tProcessing layer %lld...\n",l+1);
    }
	#endif
	process_first_layer_pvalues();
	// Artificially initialise last_tau to L-1 to ensure the first iteration of process_intervals()
	// increases the number of layers processed if the testable queue is non-empty
	last_tau = L-1;
	// Process the rest of layers (i.e. intervals of length > 1) until the pruning naturally stops the execution
	process_intervals_pvalues();
	// Report number of significant intervals found

    summaryString.append( "Number of significantly associated intervals found: " +
            AnotherToString(n_significant_intervals)  +
            "\n");
}

/* -------------------FUNCTIONS TO FIND THE CORRECTED SIGNIFICANCE THRESHOLD -------------------------------------- */

/* Decrease the minimum p-value threshold one level
 */
void decrease_threshold(){
	// Remove the intervals which become untestable after the change
	m -= freq_cnt[idx_th];
	// Change threshold
	idx_th++; pth = pgrid[idx_th];
}

void process_first_layer_threshold(){
	long long tau, j, k, queue_idx;
	char *X_tr_aux;
	long long *aux_ptr;
	double pmh_min_val;
	// Process each length 1 interval
	for(tau=0; tau<L; tau++){
		n_intervals_processed++;
		// Compute number of 1s in the interval for each of the K tables
		X_tr_aux = X_tr[tau];
		for(k=0; k<K; k++){
			aux_ptr = &freq_par[tau][k];
			for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) *aux_ptr += X_tr_aux[j];
		}
		#ifndef NO_SINGLE_FEATURES
			// If the interval is testable...
			// Update frequency-buckets and number of testable intervals
			pmh_min_val = compute_minpval(freq_par[tau]);
			if(pmh_min_val <= pth){//Definition of testability in terms of critical values
				// Increase counter of appropriate bucket and overall number of testable intervals
				freq_cnt[bucket_idx(pmh_min_val)]++; m++;
				// Update threshold until FWER condition is satisfied again
				while((m*pth) > alpha) decrease_threshold();
			}
		#endif
		// If either the current interval or the previous one are prunable
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || isprunable(freq_par[tau]) || isprunable(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}

void process_intervals_threshold(){
	long long tau, j, k, queue_idx;
	char *X_tr_aux, *X_par_aux;
	long long *aux_ptr;
	double pmh_min_val;
	// While testable-interval queue is not empty, continue to process intervals
	while(testable_queue_length){
		// Pop a testable interval from the queue
		tau = testable_queue[testable_queue_front];
		testable_queue_front = (testable_queue_front<(L-1)) ? testable_queue_front + 1: 0;
		testable_queue_length--;
		// Check if we have started processing a new layer by detecting non-monotonicity in tau
		if(tau < last_tau) {
			l++;
			#ifndef NO_VERBOSE
            if (showProcessing){
                Rcpp::Rcout << "\tProcessing layer  " << (l+1) << "...\n" << std::endl;
//                 printf("\tProcessing layer %lld...\n",l+1);
            }
			#endif
		}
		if((L_max>0) && ((l+1) > L_max)) {
			#ifndef NO_VERBOSE
            if (showProcessing){
                Rcpp::Rcout << "\tMaximum interval length achieved at l=" << (l+1) << "Stopping enumeration...\n" << std::endl;
//                 printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n",l+1);
            }
			#endif
			break;
		}
		last_tau = tau;
		// Check any of the two parents is prunable, stop processing. Notice that this check is necessary
		// even if the current interval was appended to the testable queue, because the threshold and
		// testability regions might have been modified between the time in which the current interval
		// was appended to the queue and the time in which it is being processed
		if(isprunable(freq_par[tau]) || isprunable(freq_par[tau+1])) continue;
		n_intervals_processed++;
		// Compute OR and frequency of the interval for each of the K tables
		X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
		for(k=0; k<K; k++){
			aux_ptr = &freq_par[tau][k];
			for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; (*aux_ptr)++;}
		}
		// If the interval is testable, increase counter of testable items and frequency-buckets and
		// check if the corrected significance threshold must be reduced
		pmh_min_val = compute_minpval(freq_par[tau]);
		if(pmh_min_val <= pth){//Definition of testability in terms of critical values
			// Increase counter of appropriate bucket and overall number of testable intervals
			freq_cnt[bucket_idx(pmh_min_val)]++; m++;
			// Update threshold until FWER condition is satisfied again
			while((m*pth) > alpha) decrease_threshold();
		}
		// If either the current interval or the previous one are prunable
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || isprunable(freq_par[tau]) || isprunable(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}

/* Function to give some feedback about the computation of the significance threshold */
void output_significance_threshold(){
	long long k;
    //clear the string
    summaryString.clear();
    //summary string starts
    summaryString.append("DATASET CHARACTERISTICS:\n" );

    summaryString.append("\tN = " + 
            AnotherToString(N) +   
            "\tn = " + 
            AnotherToString(n) +   
            "\tL = " + 
            AnotherToString(L) +
            "\n");

	for(k=0; k<K; k++) {
        summaryString.append("\tN[" +  AnotherToString(k) + "] = " + 
                AnotherToString(Nt[k]) +   
                ", \tn[" +  AnotherToString(k) + "] = " +
                AnotherToString(nt[k]) +   
                "\n");
    }
        
    summaryString.append("RESULTS: \n");

	// Number of intervals processed, proportion of intervals pruned
    double computedProcessed = ((double)(200*n_intervals_processed))/(L*(L+1));
    summaryString.append("Intervals processed: " +
            AnotherToString(n_intervals_processed) + 
            "(" + AnotherToString(computedProcessed) + " of total)\n");

    summaryString.append("Maximum testable interval length: " + 
            AnotherToString(l+1) );

	if(L_max==0){
        summaryString.append("Maximum interval length to be processed: unlimited\n");
    } else {
        summaryString.append("Maximum interval length to be processed: " + 
                AnotherToString(L_max) + "\n");
    }
        
    summaryString.append("Last testability threshold: " +
            AnotherToString(pth) + "\n");

    summaryString.append("Number of testable intervals: " +
            AnotherToString(m) + "\n");

    summaryString.append("Corrected significance threshold at level " +
            AnotherToString(alpha) + ":" + 
            AnotherToString(delta_opt) + "\n");
}

/* Wrapper function that encapsulates the functionality required to find the corrected significance threshold */
void compute_corrected_significance_threshold(){
	// Give feedback to user
	#ifndef NO_VERBOSE

    if (showProcessing){
        Rcpp::Rcout << "COMPUTING CORRECTED SIGNIFICANCE THRESHOLD...\n" << std::endl;
//         printf("COMPUTING CORRECTED SIGNIFICANCE THRESHOLD...\n");
    }
	#endif
	// Initialise the queue as empty
	testable_queue_front = 0; testable_queue_length = 0;
	// Initialise current layer index, current number of testable intervals and current number of intervals processed to 0
	l = 0; m = 0; n_intervals_processed = 0;
	// Initialise the value of the OR vectors of current layer to original dataset
	memcpy(X_par[0],X_tr[0],L*N*sizeof(char));
	// Process the upper-most layer (i.e. the layer composed of length 1 intervals)
	#ifndef NO_VERBOSE
    if (showProcessing){
        Rcpp::Rcout << "\tProcessing layer " << (l+1) << "...\n" << std::endl;
//         printf("\tProcessing layer %lld...\n",l+1);
    }
	#endif
	process_first_layer_threshold();
	// Artificially initialise last_tau to L-1 to ensure the first iteration of process_intervals()
	// increases the number of layers processed if the testable queue is non-empty
	last_tau = L-1;
	// Process the rest of layers (i.e. intervals of length > 1) until the pruning naturally stops the execution
	process_intervals_threshold();
	// Set final corrected significance threshold
	delta_opt = alpha/m;
	// Print results to stdout
	output_significance_threshold();
}

/*--------------------------------------------------------- FILE I/O --------------------------------------------------------*/
/* Do a first scan of the file containing the class labels to compute the total number of observations, N,
 * and the total number of observations in the positive class, n
 * */
int get_N_n(char *labels_file){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[MAX_ASCII];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
    int int_read_buf_aux;

	// Initialise both counters to 0 (the variables are defined as global variables in wy.c)
	N = 0; n = 0;

	//Try to open file, giving an error message if it fails
    try {
        if(!(f_labels = fopen(labels_file,"r"))){
            string message = "Error in function get_N_n when opening file ";
            message.append(labels_file);
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
    try {
        if(!read_buf){
            throw std::runtime_error("Error in function get_N_n: couldn't allocate memory for array read_buf\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Initialize the char to int converter
	for(i=0;i<MAX_ASCII;i++) char_to_int[i] = ASCII_OTHER;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int[ASCII_0] = 0; char_to_int[ASCII_1] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
        try {
            if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
                string message = "Error in function get_N_n while reading the file ";
                message.append(labels_file);
                message.append("\n");
                throw std::runtime_error(message);
            }
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            ::Rf_error("c++ exception (unknown reason)"); 
        }

		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
            int_read_buf_aux = (int) *read_buf_aux;
// 			if(char_to_int[*read_buf_aux] == ASCII_OTHER) continue;
			if(char_to_int[int_read_buf_aux] == ASCII_OTHER) continue;
			N++;
// 			if(char_to_int[*read_buf_aux]) n++;
			if(char_to_int[int_read_buf_aux]) n++;
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
    return EXIT_SUCCESS;
}

int read_labels_file(char *labels_file){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	long long i;// Iterator variable to be used in loops
	long long k; //Current table index
	char char_to_int[MAX_ASCII];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
	char *labels_aux = Y_tr;//Auxiliary pointer to array labels for increments
    int int_read_buf_aux;

	//Try to open file, giving an error message if it fails
    try {
        if(!(f_labels = fopen(labels_file,"r"))){
            string message = "Error in function read_labels_file when opening file ";
            message.append(labels_file);
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
    try {
        if(!read_buf){
            throw std::runtime_error("Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	//Initialize the char to int converter
	for(i=0;i<MAX_ASCII;i++) char_to_int[i] = ASCII_OTHER;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int[ASCII_0] = 0; char_to_int[ASCII_1] = 1;

	// Read the entire file
	i = 0; //Here i stands for the number of labels read so far
	k = 0; //Assume all observations have been ordered by table in the input file!!!
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
        try {
            if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
                string message = "Error in function read_labels_file while reading the file ";
                message.append(labels_file);
                message.append("\n");
                throw std::runtime_error(message);
                return EXIT_FAILURE;
            }
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            ::Rf_error("c++ exception (unknown reason)"); 
        }

		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
            int_read_buf_aux = (int) *read_buf_aux;
// 			if(char_to_int[*read_buf_aux] == ASCII_OTHER) continue;
			if(char_to_int[int_read_buf_aux] == ASCII_OTHER) continue;
			*labels_aux++ = char_to_int[int_read_buf_aux];
// 			*labels_aux++ = char_to_int[*read_buf_aux];
			nt[k] += char_to_int[int_read_buf_aux];
// 			nt[k] += char_to_int[*read_buf_aux];
			i++;
			if(i==cum_Nt[k+1]) k++;
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	// Sanity check to see if we successfully read the correct number of labels
	i = labels_aux-Y_tr;
    try {
        if(i != N){
                string message = "Error in function read_labels_file: incorrect number of labels read. Read ";
            message.append(AnotherToString(i));
            message.append(", correct number ");
            message.append(AnotherToString(N));
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
    return EXIT_SUCCESS;
}

int get_L(char *filename){
	FILE *f_dat = ((FILE*)0);
	int i, n_read;
// 	int j;
	char char_to_int[MAX_ASCII];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
    int int_read_buf_aux;

	//Try to open file, giving an error message if it fails
    try {
        if(!(f_dat = fopen(filename,"r"))){
            string message = "Error in function get_L when opening file ";
            message.append(filename);
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }
	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
    try {
        if(!read_buf){
            throw std::runtime_error("Error in function get_L: couldn't allocate memory for array read_buf\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Initialize the char to int converter
	for(i=0;i<MAX_ASCII;i++) char_to_int[i] = 0;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int[ASCII_NEWLINE] = 1;

	// Read the entire file, counting the number of lines
	L = 0;
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_dat);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
        try {
            if((n_read < READ_BUF_SIZ) && !feof(f_dat)){
                string message = "Error in function get_L while reading the file ";
                message.append(filename);
                message.append("\n");
                throw std::runtime_error(message);
            }
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            ::Rf_error("c++ exception (unknown reason)"); 
        }
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++) {
//             if(char_to_int[*read_buf_aux]) L++;
            int_read_buf_aux = (int) *read_buf_aux;
            if(char_to_int[int_read_buf_aux]) L++;
        }
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_dat)) break;
	}

	// Close file
	fclose(f_dat);
	// Free allocated memory
	free(read_buf);
    return EXIT_SUCCESS;
}



int read_dataset_file(char *filename, char *ptr){
	FILE *f_dat = ((FILE*)0);
	int i, n_read;
// 	int j;
	char char_to_int[MAX_ASCII];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
    int int_read_buf_aux;

	//Try to open file, giving an error message if it fails
    try {
        if(!(f_dat = fopen(filename,"r"))){
            string message = "Error in function get_L when opening file ";
            message.append(filename);
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }
	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
    try {
        if(!read_buf){
            throw std::runtime_error("Error in function get_L: couldn't allocate memory for array read_buf\n");
            return EXIT_FAILURE;
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Initialize the char to int converter
	for(i=0;i<MAX_ASCII;i++) char_to_int[i] = ASCII_OTHER;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int[ASCII_0] = 0; char_to_int[ASCII_1] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_dat);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
        try {
            if((n_read < READ_BUF_SIZ) && !feof(f_dat)){
                string message = "Error in function get_L while reading the file ";
                message.append(filename);
                message.append("\n");
                throw std::runtime_error(message);
            }
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            ::Rf_error("c++ exception (unknown reason)"); 
        }
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
            int_read_buf_aux = (int) *read_buf_aux;
// 			if(char_to_int[*read_buf_aux] == ASCII_OTHER) continue;
			if(char_to_int[int_read_buf_aux] == ASCII_OTHER) continue;
// 			*ptr++ = char_to_int[*read_buf_aux];
			*ptr++ = char_to_int[int_read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_dat)) break;
	}

	// Close file
	fclose(f_dat);
	// Free allocated memory
	free(read_buf);
    return EXIT_SUCCESS;
}

/* Do a first scan of the file containing the number of observations per table in order to compute the number of tables
 * */
int get_K(char *covariates_file){
	FILE *f_covariates;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[MAX_ASCII];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
    int int_read_buf_aux;

	// Initialise both counters to 0 (the variables are defined as global variables in wy.c)
	K = 0;

	//Try to open file, giving an error message if it fails
    try {
        if(!(f_covariates = fopen(covariates_file,"r"))){
            string message = "Error in function get_K when opening file ";
            message.append(covariates_file);
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
    try {
        if(!read_buf){
            throw std::runtime_error("Error in function get_K: couldn't allocate memory for array read_buf\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Initialize the char to int converter
	for(i=0;i<MAX_ASCII;i++) char_to_int[i] = ASCII_OTHER;
	// We only care about newlines
	char_to_int[ASCII_NEWLINE] = 0;

	// Read the entire file, counting the number of lines
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_covariates);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
        try {
            if((n_read < READ_BUF_SIZ) && !feof(f_covariates)){
                string message = "Error in function get_K while reading the file ";
                message.append(covariates_file);
                message.append("\n");
                throw std::runtime_error(message);
            }
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            ::Rf_error("c++ exception (unknown reason)"); 
        }
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
            int_read_buf_aux = (int) *read_buf_aux;
			//If the character is not a newline process the next character
// 			if(char_to_int[*read_buf_aux] == ASCII_OTHER) continue;
			if(char_to_int[int_read_buf_aux] == ASCII_OTHER) continue;
			K++;
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_covariates)) break;
	}

	//Close the file
	fclose(f_covariates);

	//Free allocated memory
	free(read_buf);
    return EXIT_SUCCESS;
}

int read_covariates_file(char *covariates_file){
	FILE *f_covariates;//Stream with file containing class labels
	int n_read;//Number of chars read
	long long i;// Iterator variable to be used in loops
	long long k;//Number of tables already processed
	//char c;// Iterator variable to be used in loops
	int c;// Iterator variable to be used in loops
	char char_to_int[MAX_ASCII];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
    int int_read_buf_aux;

	//Try to open file, giving an error message if it fails
    try {
        if(!(f_covariates = fopen(covariates_file,"r"))){
            string message = "Error in function read_covariates_file when opening file ";
            message.append(covariates_file);
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }


	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
    try {
        if(!read_buf){
            throw std::runtime_error("Error in function read_covariates_file: couldn't allocate memory for array read_buf\n");
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Initialize the char to int converter
	for(i=0;i<MAX_ASCII;i++) char_to_int[i] = ASCII_OTHER;
	// We only care about chars representing digits and newline
	for(c=ASCII_0; c<=ASCII_9; c++) char_to_int[c] = c - ASCII_0;
	char_to_int[ASCII_NEWLINE] = ASCII_OTHER_NEWLINE;

	// Read the entire file
	i = 0; k = 0;
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_covariates);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
        try {
            if((n_read < READ_BUF_SIZ) && !feof(f_covariates)){
                string message = "Error in function read_covariates_file while reading the file ";
                message.append(covariates_file);
                message.append("\n");
                throw std::runtime_error(message);
            }
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            ::Rf_error("c++ exception (unknown reason)"); 
        }
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is neither a digit nor a newline process the next char

            int_read_buf_aux = (int) *read_buf_aux;
// 			if(char_to_int[*read_buf_aux] == ASCII_OTHER) continue;
			if(char_to_int[int_read_buf_aux] == ASCII_OTHER) continue;
			// If the character is a newline, we have read a number already so we save it and go to next line
// 			if(char_to_int[*read_buf_aux] == ASCII_OTHER_NEWLINE){
			if(char_to_int[int_read_buf_aux] == ASCII_OTHER_NEWLINE){
				Nt[k++] = i;
				cum_Nt[k] = cum_Nt[k-1] + Nt[k-1];
				i = 0;
				continue;
			}
			// Otherwise the character is a digit, so we accumulate it into the current number
// 			i = 10*i + char_to_int[*read_buf_aux];
			i = 10*i + char_to_int[int_read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_covariates)) break;
	}

	// Sanity check to see if we successfully read the distribution of observations per table
	i = 0;
	for(k=0; k<K; k++) i += Nt[k];
    try {
        if(i != N){
            string message = "Error in function read_covariates_file: incorrect number of observations per table read. Total N ";
            message.append(AnotherToString(N));
            message.append(", Accumulated N in covariates file ");
            message.append(AnotherToString(i));
            message.append("\n");
            throw std::runtime_error(message);
        }
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

	//Close the file
	fclose(f_covariates);

	//Free allocated memory
	free(read_buf);

    return EXIT_SUCCESS;
}


inline int bucket_idx(double pval){
	int idx;
	idx = (int)floor(-log10(pval)/log10_p_step);
	if(idx<0) idx = 0;
	if(idx>NGRID) idx = NGRID;
	return idx;
}

int qsort_cmp_betas ( const void *x, const void *y ){
	if ( betas[*((long long *)x)] < betas[*((long long *)y)] ) return (-1);
	else return 1;
}













/* ----------------------------------------------------ENTRY POINT---------------------------------------------------------- */

// void sis_init(char *X_filename, char *Y_filename, char *C_filename, double target_fwer, long long l_max){

// 		printf("\tUSAGE: ./program_name X_file Y_file C_file alpha L_max base_filename [-postprocessing_folder path_to_pyfiltering.py] [-pval_file all_pvals_file]\n");


//the main computation of FastCMH is here
int computeFastCMH(char* xfilenameCpp, char* yfilenameCpp, char* cfilenameCpp, double alphaval, long long Lmaxlonglong){
	// Get time when program started to run
	t_init = measureClocks();

	// INITIALISATION
	tic = measureClocks();

    int initSuccess = EXIT_FAILURE;
    try {
        initSuccess = sis_init(xfilenameCpp, yfilenameCpp, cfilenameCpp, alphaval, Lmaxlonglong);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

    if (initSuccess==EXIT_SUCCESS){
        toc = measureClocks();
        time_initialisation = (toc-tic) / CLOCKS_PER_SEC;

        // Main functionality
        tic = measureClocks();
        compute_corrected_significance_threshold();
        toc = measureClocks();
        time_comp_threshold = (toc-tic) / CLOCKS_PER_SEC;
        tic = measureClocks();
        find_significant_intervals();
        toc = measureClocks();
        time_comp_significant_intervals = (toc-tic) / CLOCKS_PER_SEC;



        // Get time when program finished
        t_end = measureClocks();

        // Produce output and free memory
        sis_end();
    } else {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS; 
}



//run FDR on the FastCMH
//only alphaval is a parameter
void computeFDR_ForFastCMH(double alphaval){
        std::vector<long long> perm =  gilbertFDR(allTestablePval, allTestableTau, allTestableL, alphaval, useDependenceFDR);
        fdrPval = extractFdrPvalue(allTestablePval, perm);
        fdrTau = extractFdrTau(allTestableTau, perm);
        fdrL = extractFdrL(allTestableL, perm);
}



//create the timing list
//all variables are global, so no need for arguments
Rcpp::List createTimingList(){
    Rcpp::List timingList = Rcpp::List::create(Rcpp::Named("details") = timingString,
            Rcpp::Named("exec") = timeExecution,
            Rcpp::Named("init") = timeInitialisation,
            Rcpp::Named("fileIO") = timeFileIO,
            Rcpp::Named("compSigThresh") = timeComputeSigThreshold,
            Rcpp::Named("compSigInt") = timeComputeSigIntervals);
//             Rcpp::Named("peakMemUsage") = peakMemoryUsageInBytes);
    return timingList;
}



//Function that creates Rcpp list from vectors
//Want to move this bulky code out of the main
//only argument is filteredInt
//other values are all global
Rcpp::List createReturnListNoFDR(vector<Interval> filteredInt){
    //create list
    Rcpp::List returnList;

    //create timing list
    Rcpp::List timingList = createTimingList();

    //now wrap in Rcpp dataframe
    Rcpp::DataFrame filteredDf = extractDataFrameFromIntervalVector(filteredInt);

    //now wrap sigDf
    Rcpp::DataFrame sigDf = createDataFrameTauLPvalue(sigTau, sigL, sigPval);

    //now wrap all intervals
    Rcpp::DataFrame allTestableDf = createDataFrameTauLPvalue(allTestableTau, allTestableL, allTestablePval);



    //saveAllPvals is actually global
    if (saveAllPvals){
        returnList = Rcpp::List::create(Rcpp::Named("sig") = filteredDf,
                                            Rcpp::Named("unfiltered") = sigDf,
                                            Rcpp::Named("allTestable") = allTestableDf,
                                            Rcpp::Named("histObs") = histObs,
                                            Rcpp::Named("histFreq") = histFreq,
                                            Rcpp::Named("summary") = summaryString,
                                            Rcpp::Named("timing") = timingList);
    } else {
        returnList = Rcpp::List::create(Rcpp::Named("sig") = filteredDf,
                                            Rcpp::Named("unfiltered") = sigDf,
                                            Rcpp::Named("histObs") = histObs,
                                            Rcpp::Named("histFreq") = histFreq,
                                            Rcpp::Named("summary") = summaryString,
                                            Rcpp::Named("timing") = timingList);
    }

    return returnList;
}



//Function that creates Rcpp list from vectors WITH FDR
//Want to move this bulky code out of the main
//other values are all global
Rcpp::List createReturnListWithFDR(vector<Interval> filteredInt, vector<Interval> filteredIntFDR){
    //create list
    Rcpp::List returnList;

    //create timing list
    Rcpp::List timingList = createTimingList();

    //now wrap in Rcpp dataframe
    Rcpp::DataFrame filteredDf = extractDataFrameFromIntervalVector(filteredInt);

    //now wrap sigDf
    Rcpp::DataFrame sigDf = createDataFrameTauLPvalue(sigTau, sigL, sigPval);

    //now wrap all intervals
    Rcpp::DataFrame allTestableDf = createDataFrameTauLPvalue(allTestableTau, allTestableL, allTestablePval);


    //filter FDR values and extract dataframe
    Rcpp::DataFrame filteredDfFDR = extractDataFrameFromIntervalVector(filteredIntFDR);
    Rcpp::DataFrame fdrDf = createDataFrameTauLPvalue(fdrTau, fdrL, fdrPval);



    //now include FDR filtered
    //saveAllPvals is actually global
    if (saveAllPvals){
        returnList = Rcpp::List::create(Rcpp::Named("sig") = filteredDf,
                                            Rcpp::Named("unfilteredSig") = sigDf,
                                            Rcpp::Named("fdr") = filteredDfFDR,
                                            Rcpp::Named("unfilteredFdr") = fdrDf,
                                            Rcpp::Named("allTestable") = allTestableDf,
                                            Rcpp::Named("histObs") = histObs,
                                            Rcpp::Named("histFreq") = histFreq,
                                            Rcpp::Named("summary") = summaryString,
                                            Rcpp::Named("timing") = timingList);
    } else {
        returnList = Rcpp::List::create(Rcpp::Named("sig") = filteredDf,
                                            Rcpp::Named("unfilteredSig") = sigDf,
                                            Rcpp::Named("fdr") = filteredDfFDR,
                                            Rcpp::Named("unfilteredFdr") = fdrDf,
                                            Rcpp::Named("histObs") = histObs,
                                            Rcpp::Named("histFreq") = histFreq,
                                            Rcpp::Named("summary") = summaryString,
                                            Rcpp::Named("timing") = timingList);
    }

    return returnList;
}





//This is the main function called from R
//1. Unwraps variables from R (and recasts them)
//2. Calls computeFastCMH, which computes necessary values
//3. Filtering function is called
//4. Wraps variables for R
// [[Rcpp::export]]
Rcpp::List main_fastcmh2(Rcpp::String xfilenameR, 
                         Rcpp::String yfilenameR, 
                         Rcpp::String cfilenameR, 
                         Rcpp::NumericVector alphaR,
                         Rcpp::NumericVector lmaxR,
                         Rcpp::LogicalVector showProcessingR,
                         Rcpp::LogicalVector saveAllPvalsR,
                         Rcpp::LogicalVector doFDR_R, 
                         Rcpp::LogicalVector useDependenceFDR_R){

    //First thing to do: clear vectors
    sigTau.clear();
    sigL.clear();
    sigPval.clear();

    allTestableTau.clear();
    allTestableL.clear();
    allTestablePval.clear();

    histObs.clear();
    histFreq.clear();

    fdrTau.clear();
    fdrL.clear();
    fdrPval.clear();


    //----------------------------------------------//
    //
    const size_t MAXSTRINGSIZE = 2000 * sizeof(char);
    //1. Unwraps variables from R (and recasts them)
    //now conversions
    const char* xfilenameCppConst = xfilenameR.get_cstring();
    char xfilenameCpp[MAXSTRINGSIZE];
    strcpy(xfilenameCpp, xfilenameCppConst);

    const char* yfilenameCppConst = yfilenameR.get_cstring();
    char yfilenameCpp[MAXSTRINGSIZE];
    strcpy(yfilenameCpp, yfilenameCppConst);

    const char* cfilenameCppConst = cfilenameR.get_cstring();
    char cfilenameCpp[MAXSTRINGSIZE];
    strcpy(cfilenameCpp, cfilenameCppConst);

    int Lmaxint = Rcpp::as<int> (lmaxR);
    long long Lmaxlonglong = (long long) Lmaxint;

    double alphaval = Rcpp::as<double> (alphaR);

//     const char* sigintfilenameCppConst = sigintfilenameR.get_cstring();
//     char sigintfilenameCpp[MAXSTRINGSIZE]; 
//     strcpy(sigintfilenameCpp, sigintfilenameCppConst);

    //convert to boolean - global variables controlling printing to screen and
    //whether or not all pvals are saved
    showProcessing = Rcpp::as<bool>(showProcessingR);
    saveAllPvals = Rcpp::as<bool>(saveAllPvalsR);
    doFDR = Rcpp::as<bool>(doFDR_R);
    useDependenceFDR = Rcpp::as<bool>(useDependenceFDR_R);


    //----------------------------------------------//
    //2. Calls computeFastCMH, which computes necessary values
    int fastcmh_SUCCESS = EXIT_FAILURE;
    try {
        //here will change to EXIT_SUCCESS, if there is success
        fastcmh_SUCCESS = computeFastCMH(xfilenameCpp, yfilenameCpp, cfilenameCpp, alphaval, Lmaxlonglong);
    } catch(std::exception &ex) {	
        forward_exception_to_r(ex);
//         return EXIT_FAILURE;
//         Do nothing
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }

    //create an empty returnList, which will be modified below
    Rcpp::List returnList;

    if (fastcmh_SUCCESS==EXIT_SUCCESS){
        //----------------------------------------------//
        //3. Filtering function is called
        //TODO: fix size of string? 512?
        //now filter
        if (showProcessing){
            Rcpp::Rcout << "Filtering overlapping intervals...\n\n " << std::endl;
//             printf("Filtering overlapping intervals...\n\n");
        }
        vector<Interval> filteredInt = cpp_filterIntervalsFromMemory(sigTau, sigL, sigPval);

        //----------------------------------------------//
        //4. Do FDR (maybe)
        vector<Interval> filteredIntFDR;
        if (doFDR){
            computeFDR_ForFastCMH(alphaval);
            filteredIntFDR = cpp_filterIntervalsFromMemory(fdrTau, fdrL, fdrPval);
        }

        //----------------------------------------------//
        //5. Wraps variables for R

        if (doFDR){
            returnList = createReturnListWithFDR(filteredInt, filteredIntFDR);
        } else {
            returnList = createReturnListNoFDR(filteredInt);
        } 
    } else {
        //return a list with an error message
        returnList = createErrorReturnList();
    }


    //----------------------------------------------//
    return returnList;
//     Rcpp::Rcout << "exit success from computeFastCMH" << std::endl;
}
