/*
 * Author: Luke Miratrix (lmiratrix@stat.harvard.edu)
 * This is a rewrite of code from the ngram package of Georgiana Ifrim (georgiana.ifrim@gmail.com)
 *
 */

#include <Rcpp.h>

#include <cfloat>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <functional>
#include <map>
#include <set>
#include <iterator>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include "sys/time.h"
#include <list>

// #if defined(__i386__) && _WINDOWS_
#ifdef _WIN32
 	#define LONG_DOUBLE double
#else
 	#define LONG_DOUBLE long double
#endif

#define SINGLE_FILE 1
#define MULTI_FILE 2
#define FILE_INPUT_MODE MULTI_FILE // true means single-file, false means read in labeling directory

#define BFS 0
#define DFS 1

#define CHAR_TOKEN_T 1
#define WORD_TOKEN 0

// this makes sure that a gradient has to be at least this much better than prior one
// to get counted.
#define TINY_WEIGHT 0

//#define TINY_WEIGHT 0.00000001

#define LP_INFINITY 10

#define INTERCEPT_STRING "*intercept*"

using namespace std;

using namespace Rcpp;


// #ifdef NDEBUG
// # define // [Fix] my_assert(EX)
// #else
// # define // [Fix] my_assert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))
// #endif
//
// void __assert (const char *msg, const char *file, int line) {
//     char buffer [100];
//     snprintf( buffer, 100, "Assert Failure: %s at %s line #%d", msg, file, line );
//     ::Rf_error( buffer );
// }


//void // [Fix] my_assert( string name, bool what ) {
//	if (!what) {
//		::Rf_error( ("assert failure: " + name).c_str() );
//	}
//}







// The search space rooted at this ngram.
class space_t {

public:

	// docid of last document added via add()
	int last_docid;

	// Pointer to previous ngram.
	space_t *prev;

	// Label of node where extension is done.
	string ne;

	// Vector of ngrams which are extensions of current ngram.
	std::vector <space_t *> children;

	// Number of tokens of ngram
	unsigned int size;

	// Gradient value for this ngram.
	LONG_DOUBLE gradient;
	LONG_DOUBLE upos;
	LONG_DOUBLE uneg;

	// Full Ngram name
	std::string  ngram;

	// Ngram support, e.g. docids where it occurs in the collection.
	std::vector <unsigned int> doc_support;

	// total number of times ngram appears across corpus.
	unsigned int total_support;

	// total number of documents ngram appears in
	unsigned int total_docs;

	// Ngram support weights (i.e., number of times ngram appears in each support document)
	std::vector <int> weight;

	// normalization.  All weights are divided by this.
	LONG_DOUBLE Z;

	// penalize rule?  (generally don't for intercept, do for all others)
	bool penalize;

private:
	// Total list of occurrences in the entire collection for this ngram. A sort of expanded inverted index.
	// negative numbers are document ids (flip and subtract 1).  Pos numbers are character locations
	// for _endpoint_ of the ngram.
	std::vector <int>       loc;

	// true if calc_support_weights() has been called on this.
    bool converted;

	// has shrink been called (which destroys locations of features)
    bool shrunk;

public:

	friend bool operator < (const space_t &r1, const space_t &r2) {
		return r1.ngram < r2.ngram;
	}


    bool isConverted() {
        return converted;
    }

    bool isShrunk() {
        return shrunk;
    }


    std::vector <int>  &get_loc() {
        // [Fix] my_assert( !shrunk );
//		// [Fix] my_assert( "hasn't shrunk in get_loc",  !shrunk );

    	return loc;
    }


	// Add a doc_id and position of occurrence to the total list of occurrences,
	// for this ngram.
	// Encode the doc_id in the vector of locations.
	// Negative entry means new doc_id.
	void add (unsigned int docid, int pos) {
		//Rcout << "adding " << docid << "/" << pos << " to " << ne << endl;
		if (last_docid != (int)docid) {
			loc.push_back (-(int)(docid+1));
			total_docs++;
		}
		loc.push_back (pos);
		last_docid = (int)docid;
		total_support++;
		//Rcout << "\ndocid:" << docid << " pos:" << pos;
	}

	// figure out which doc ids are the support of this space, and also how
	// many times the ngram is found in each doc in the support.
	// Also compute the normalizing constant Z for this feature.
	void calc_support_weights( LONG_DOUBLE Lp, bool binary_features, bool no_regularization ) {
        // [Fix] my_assert( !shrunk );
        // [Fix] my_assert( !converted );
		int rcntr = 0;
		unsigned int my_total_support = 0;
		doc_support.clear();
		weight.clear();
		for (unsigned int i = 0; i < loc.size(); ++i) {
			//Rcout << loc[i] << " rcntr " << rcntr << endl;
			if (loc[i] < 0) {
				doc_support.push_back( (unsigned)(-loc[i] - 1 ));
				weight.push_back( 0 );
				rcntr++;
			} else {
				if ( binary_features ) {
					if ( weight[rcntr-1] == 0 ) {
						weight[rcntr-1] = 1;
						my_total_support++;
					}
				} else {
					weight[rcntr-1]++;
					my_total_support++;
				}
			}
		}
		//Rcout << total_support << " " << my_total_support << endl;
		if ( binary_features ) {
			// [Fix] my_assert( total_docs == my_total_support );

		} else {
			// [Fix] my_assert( total_support == my_total_support );
		}
		// now re-weight based on total support (this is L2 normalization)
		Z = 0.0;
        if ( no_regularization ) {
            Z = 1.0;
        } else if ( Lp >= LP_INFINITY ) {
			for ( unsigned int i = 0; i < weight.size(); i++ ) {
				if ( weight[i] > Z ) {
					Z = weight[i];
				}
			}
		} else {
			for ( unsigned int i = 0; i < weight.size(); i++ ) {
				Z += pow( weight[i], Lp );
			}
			Z = pow( (LONG_DOUBLE)Z, 1/Lp );
		}

		//Rcout << "Condensing vectors" << endl;
		std::vector<int>(weight).swap (weight); // shrink
		std::vector<unsigned int>(doc_support).swap (doc_support); // shrink
		converted = true;
	}

	// Shrink the list of total occurrences to contain just support doc_ids.
	void shrink () {
		// [Fix] my_assert( !shrunk );
		std::vector<int> tmp;
		for (unsigned int i = 0; i < loc.size(); ++i) {
			if (loc[i] < 0) tmp.push_back (loc[i]);
		}
		loc = tmp;
		std::vector<int>(loc).swap (loc); // shrink
		last_docid = -1;
		shrunk = true;
	}

	// Return the support of current ngram.
	// Simply count the negative loc as doc_ids.
	unsigned int support () {
		// [Fix] my_assert( converted );
		return total_support;
		/*unsigned int result = 0;
         for (unsigned int i = 0; i < loc.size(); ++i) {
         if (loc[i] < 0) ++result;
         }
         return result;*/
	}

	void set_ngram_string( ) {
		ngram = ne;
		for (space_t *t = prev; t; t = t->prev) {
			//	ngram = t->ne + ":" + ngram;
			ngram = t->ne + " " + ngram;
		}
	}

	void set_ne( const string &new_ne ) {
		ne = string(new_ne);
		set_ngram_string();
	}

	void debug_print_support( ) {
		Rcout << ngram << ": ";
		if ( total_support < 50 ) {
			for (unsigned int i = 0; i < loc.size(); ++i) {
				Rcout << loc[i] << " ";
			}
		} else {
			Rcout << total_support << " documents";
		}
		Rcout << "\n";
	}

	void print_rule( ) {
		Rcout << "RULE: '" << ngram << "' gr: " << gradient << "\n\tSupport:";
		if ( ngram != INTERCEPT_STRING ) {
			unsigned int cut = doc_support.size();
			if ( cut > 10 ) { cut = 10; }
			for ( unsigned int i = 0; i < cut; i++ ) {
				Rcout << " " << doc_support[i] << "(" << weight[i] << ")";
			}
			if ( doc_support.size() > cut ) {
				Rcout << " ...";
			}
		}
		Rcout << "\n\tZ: " << Z;
		Rcout << endl;
	}

	void print_space( ) {
		Rcout << "SPACE: '";
		Rcout << ngram << "' last docID: " << last_docid << " prev: " << prev << "\n\tFull Support:";
		for ( unsigned int i = 0; i < loc.size(); i++ ) {
			Rcout << " " << loc[i] ;
		}
		Rcout << endl;
	}



	// Constructor for space_t.   Copy clone to us, and set our parent to parent.
	space_t( space_t &clone, space_t *parent, const string &tok );
	//space_t( const string &tok );
	space_t();


	Rcpp::List to_Rcpp_list() {
		Rcpp::NumericVector suppV(doc_support.size());
		for ( unsigned int i = 0; i < doc_support.size(); i++ ) {
			suppV[i] = doc_support[i];
		}
		Rcpp::NumericVector weightV(weight.size());
				for ( unsigned int i = 0; i < weight.size(); i++ ) {
					weightV[i] = weight[i];
		}

		Rcpp::List lst = Rcpp::List::create( Rcpp::Named("ngram") = CharacterVector(ngram),
	                                        Rcpp::Named("support") = suppV,
	                                        Rcpp::Named("weight") = weightV );

		return lst;

	}


	/**
	 * How many wildcards are at the tail of this feature?
	 */
	unsigned int num_consec_gaps( ) {
		unsigned int num_consec_gaps = 0;
		for (space_t *t = this; t; t = t->prev) {
			if (t->ne.compare("*") == 0)
				num_consec_gaps ++;
			else break;
		}

		return num_consec_gaps;
	}


};

space_t::space_t( space_t &clone, space_t *parent, const string &tok ) {
	last_docid = clone.last_docid;
	prev = parent;
	ne = tok;
	children.clear();
	gradient = clone.gradient;
	//size = 0;
	ngram=clone.ngram;
	doc_support.clear();
	weight.clear();
	Z = 0;
	penalize = true;
	converted = false;
	loc = clone.loc;
	total_support = clone.total_support;
	total_docs = clone.total_docs;
	set_ngram_string();
	std::vector<int>(loc).swap(loc);
	size = parent->size + 1;
	shrunk = false;
}

space_t::space_t() {
	uneg = upos = 0;
	last_docid = -1;
	prev = 0;
	total_support = 0;
	total_docs = 0;
	ne = "";
	loc.clear();
	children.clear();
	gradient = 0;
	//size = 0;
	ngram="";
	doc_support.clear();
	weight.clear();
	Z = 0;
	penalize = true;
	converted = false;
	size = 1;
	shrunk = false;
}



ostream & operator << ( ostream & os, const space_t & s ) {

	return os << "[" << s.ngram << "]";

}
ostream & operator << ( ostream & os, const space_t * s ) {

	if ( s != NULL ) {
		return os << "[" << s->ngram << "]";
	} else {
		return os << "NONE";
	}

}




const char *DELIMETERS = "\t ";



class SeqLearner {


	//private:
public:


	// Entire collection of documents, each doc represented as a string.
	// The collection is a vector of strings.
	std::vector < string > corpus;


	// True classes.
	std::vector < int >                 y;

	// Document weights (how much each document 'counts')   1 for everything is the norm.
	//std::vector < LONG_DOUBLE >         doc_weight;

	//  number of positive and negative examples
	unsigned int num_pos;
	unsigned int num_neg;

	// The fraction: 1 / 1 + exp^(yi*beta^t*xi) in the gradient computation.
	std::vector < LONG_DOUBLE >         exp_fraction;
	// Per document, sum of best beta weights, beta^t * xi = sum_{j best beta coord} gradient_j
	std::vector < LONG_DOUBLE >              sum_best_xbeta;


	// Regularized loss function: loss + lambda * elnreg
	// SLR loss function: log(1+exp(-yi*beta^t*xi))
	// Squared Hinge SVM loss function: sum_{i|1-yi*beta^t*xi > 0} (1 - yi*beta^t*xi)^2
	//LONG_DOUBLE loss;

	// cache of all rules used in model.
	std::map <string,space_t *> used_rule_cache;

	std::map <string, LONG_DOUBLE> final_model_cache;

	std::map <string, LONG_DOUBLE> features_cache;
	//map<string, LONG_DOUBLE>::iterator features_it;

	// Objective function. For now choice between logistic regression and l2 (Squared Hinge Loss) or l1 (Hinge loss) SVM loss.
	unsigned int objective;
	// Regularizer value.
	LONG_DOUBLE C;
	// Weight on L1 vs L2 regularizer.
	LONG_DOUBLE alpha;

	// p for the Lp norm on the regularization (and complemntary q)
	LONG_DOUBLE Lp;
	LONG_DOUBLE Lq;

    // true means don't use counts, just use existance of feature in document.
    bool binary_features;

    // do not regularize the features
    bool no_regularization;

	// The sum of squared values of all non-zero beta_j.
	LONG_DOUBLE sum_squared_betas;

	// The sum of abs values of all non-zero beta_j.
	LONG_DOUBLE sum_abs_betas;

    // set of unigrams with sufficient support (to check possible extensions of ngram phrase)
	std::set <string> single_node_minsup_cache;
    std::map <string, space_t> unigrams;

	// tokens that are banned from consideration.  No phrase
	// can cross them, and they cannot be components of phrases.
	std::set <string> banned_words;


	// The intercept rule
	space_t		 int_rule;

	// Pointer to current rule.
	space_t      *rule;
	// Current suboptimal gradient.
	LONG_DOUBLE       tau;



	// Whether to force betas to be nonnegative (1 means yes)
	bool pos_only;


	// How much to inflate the weight of positive features.
	LONG_DOUBLE positive_weight;

	// Max length for an ngram.
	unsigned int maxpat;
	// Min length for an ngram.
	unsigned int minpat;
	// Min supoort for an ngram.
	unsigned int minsup;

	// Allow features with gaps. Treat the gap as an additional unigram.
	// Max # of consecutive gaps allowed in a feature.
	unsigned int maxgap;

	// Total number of nodes in stored ngram tree
	unsigned int total_nodes;

	// Total number of times the pruning condition checked in current step.
	unsigned int total;
	// Total number of times the pruning condition is satisfied in current step.
	unsigned int pruned;
	// Total number of times the best rule is rewritten in current step.
	unsigned int rewritten;

	// Convergence threshold on aggregated change in score predictions.
	// Can be used as a regularization parameter.
	LONG_DOUBLE convergence_threshold;

	// Verbosity level: 0 - print no information,
	// 					1 - print profiling information,
	//					2 - print statistics on model and loss ftc per iteration
	//					> 2 - print details about search for best n-gram and pruning process
	int verbosity;
	int step_verbosity;

	// Type of token
	bool token_type;

	// Traversal strategy: BFS or DFS.
	bool traversal_strategy;

	// Profiling variables.
	struct timeval t;
	struct timeval t_origin;


	// number of iterations run for gradient descent
	unsigned int itr;

	//LONG_DOUBLE LDBL_MAX = numeric_limits<LONG_DOUBLE>::max();

	// for keeping history of model
	std::vector < LONG_DOUBLE >    history_steps;
	std::vector < string >         history_ngrams;



	SeqLearner() {
		// some defaults
        verbosity = 1;
        objective = 2;
		pos_only = 0;
		positive_weight = 1.0;
		maxpat = 0xffffffff;
		minpat = 1;
		minsup = 1;
		maxgap = 0;
        token_type = WORD_TOKEN;
		traversal_strategy = BFS;
		// Set the convergence threshold as in paper by Madigan et al on BBR.
		convergence_threshold = 0.005;
        C = 0.0;
        alpha = 1;
        Lp = Lq = 2;
        step_verbosity = 1;

        banned_words.clear();
        y.clear();
        corpus.clear();

        history_steps.clear();
        history_ngrams.clear();

        reset_tree();

        gettimeofday(&t_origin, NULL);

	}

    // have the unigrams been culled?
    bool has_unigram_cache() {
        return !unigrams.empty();
    }


    // Reset tree for a new textreg run.
    void reset_tree() {
        unigrams.clear();
        single_node_minsup_cache.clear();

        gettimeofday(&t_origin, NULL);
    }


    void set_Lp( double new_Lp ) {
		// [Fix] my_assert( new_Lp >= 1 );

    	Lp = new_Lp;
    	if ( Lp >= LP_INFINITY ) {
            Lp = LP_INFINITY;
    		Lq = 1.0;
    	} else {
    		Lq = 1.0 / (1.0 - 1.0/Lp);
    	}

        if ( verbosity >= 2 ) {
            Rcout << "\nLp set to " << Lp << " with Lq = " << Lq << endl;
        }
    }

	void print_settings( unsigned int maxitr ) {

		Rcout << "\nParameters used: " << "\n\tobjective fct: " << objective << " (0 = L1LR, 1 = SVM, 2 = Hinge w/ Lasso)\n\tT: " << maxitr
        << "\n\tminpat: " << minpat << " maxpat: " << maxpat << " minsup: " << minsup
        << "\n\tpos weight: " << positive_weight << " (inflate positive features)"
        << "\n\tmaxgap: " << maxgap << "\n\ttoken_type: " << (token_type == WORD_TOKEN ? "word" : "character")
        << "\n\ttraversal_strategy: " << (traversal_strategy == BFS ? "BFS" : "DFS")
        << "\n\tconvergence_threshold: " 	<< convergence_threshold
        << "\n\tC (regularizer value): " << C
        << "\n\tLp / Lq: " << Lq << " / " << Lp // switched on purpose.  Legacy code.
        << "\n\tb (binary only): " << binary_features
        << "\n\tn (no regularization): " << no_regularization
        << "\n\tp (positive only): " << pos_only
        << "\n\talpha (weight on l1_vs_l2_regularizer): " << alpha  << "\n\tverbosity: " << verbosity
        << "\n\tFile Mode: " << FILE_INPUT_MODE << " (single file = " << SINGLE_FILE << ")"<< endl;
        Rcout.flush();
	}


	/*
	 * Print out the loss for each document.  Also print out total loss (with penalty).
	 */
	void print_estimates( vector<LONG_DOUBLE>& sum_betas, bool print_cache ) {
		LONG_DOUBLE loss = 0;
		LONG_DOUBLE tloss;
		Rcout << "n\ty\txbeta\t\tloss" << endl;

		for ( unsigned int i = 0; i < corpus.size(); i++ ) {
			// Compute loss.
			tloss = 0;
			if (objective == 0) { //SLR
				exp_fraction[i] = calc_exp_fraction( y[i], sum_betas[i] );
				if (-y[i] * sum_betas[i] > 8000) {
					tloss = log(LDBL_MAX);
				} else {
					tloss = log(1 + exp(-y[i] * sum_betas[i]));
				}
			} //end SLR
			if (objective == 1) { //L1-SVM (Hinge loss)
				if (1 - y[i] * sum_betas[i] > 0) {
					tloss = 1 - y[i] * sum_betas[i];
				}
			} //end L1-SVM
			if (objective == 2) { //L2-SVM
				if (1 - y[i] * sum_betas[i] > 0) {
					tloss = pow(1 - y[i] * sum_betas[i], 2);
				}
			} //end L2-SVM
			loss += tloss;
			Rcout  << i << "\t" << y[i] << "\t" << sum_betas[i] << "\t" << tloss << endl;
		}
		Rcout << "Total loss: " << loss << endl;

		if ( print_cache ) {
			Rcout << "Feature cache betas:" << endl;
			if ( features_cache.empty() ) {
				Rcout << "\tempty" << endl;
			}
			LONG_DOUBLE sum_squared_betas = 0.0;
			LONG_DOUBLE sum_abs_betas = 0.0;
			for(std::map<string,LONG_DOUBLE>::iterator iter = features_cache.begin(); iter != features_cache.end(); ++iter) {
				sum_squared_betas += pow( iter->second, 2 );
				sum_abs_betas += abs( iter->second );
				Rcout << "\t" << iter->first << "\t" << iter->second << endl;
			}
		}
		if ( C != 0 ) {
			LONG_DOUBLE pen = calc_penalty( sum_abs_betas, sum_squared_betas, 0, 0 );
			Rcout << "Loss + C*Penalty: " <<  setprecision(3) << loss << " + " << setprecision(3) << C << " * " << setprecision(3) << (pen/C) << " = " << (pen+loss) << endl;
		}
	}



	/**
	 * Calculate some summary statistics for the given feature 'rule'
	 */
	void print_rule_stats( space_t *rule, LONG_DOUBLE beta, std::ostream &os ) {
		int pos_count = 0;
		//int pos_word_count = 0;
		int neg_count = 0;
		//int neg_word_count = 0;
		for ( unsigned int i = 0; i < rule->doc_support.size(); i++ ) {
			unsigned int dc = rule->doc_support[i];
			if ( y[dc] == 1 ) {
				pos_count++;
				//pos_word_count += rule->weight[i];
			} else {
				neg_count++;
				//neg_word_count += rule->weight[i];
			}
		}
		float per = (float)pos_count / (rule->total_docs);
		float perpos = (float)pos_count / num_pos;
		float perneg = (float)neg_count / num_neg;
        //	float odds = (float)pos_count * ((float)(num_neg - neg_count)) / ( neg_count * ((float)(num_pos - pos_count ) ) );
	//	float word_per = (float)pos_word_count / (rule->total_support);
		//float word_odds = (float)pos_word_count * ((float)(rule->total_support - neg_word_count)) / ( neg_word_count * ((float)(rule->total_support - pos_word_count ) ) );

		os << setprecision(3) << beta << '\t' << rule->Z << '\t'  << setprecision(3) << (beta / rule->Z) << '\t' << rule->total_support << "\t" << rule->total_docs << "\t";
		os << pos_count << '\t' << neg_count << '\t' << setprecision(3) << per << '\t' << setprecision(3) << perpos << '\t' << setprecision(3) << perneg << '\t';
	//	os << pos_word_count << '\t' << neg_word_count << '\t' << setprecision(3) << word_per << '\t';
		os << rule->ngram << endl;
	}

	/**
	 * Print out all the rules in the model with details about those rules
	 * in a tab-separated form.
	 */
	void print_full_model( std::map<string,LONG_DOUBLE> &model, unsigned int itr, unsigned int maxitr, std::ostream &os  ) {
		print_ban_list( os );
		os << "END BAN\n";
		os << "# Positive: " << num_pos << endl;
		os << "# Negative: " << num_neg << endl;
		os << "Total documents: " << corpus.size() << endl;
		os << "Iterations: " << itr << "/" << maxitr << endl;
		os << "START MODEL" << endl;
		os << "beta\tZ\tweight\tsupport\tdoc_sup\tpos\tneg\tper\tperpos\tperneg\tpos w\tneg w\tper w\tngram\n";
		// negative ones only
		for(std::map<string,LONG_DOUBLE>::iterator iter = model.begin(); iter != model.end(); ++iter) {
			// [Fix] my_assert( used_rule_cache.count( iter->first )  > 0 );
			space_t *rule = used_rule_cache[iter->first];
			//rule->debug_print_support();
			if ( iter->second <= 0 ) {
				print_rule_stats( rule, iter->second, os );
			}
		}
		// now positive
		for(std::map<string,LONG_DOUBLE>::iterator iter = model.begin(); iter != model.end(); ++iter) {
			space_t *rule = used_rule_cache[iter->first];
			if ( iter->second > 0 ) {
				print_rule_stats( rule, iter->second, os );
			}
		}

	}


	/**
	 * Make Rcpp List of all the updates (time series of path taken)
	 */
	Rcpp::List make_search_path( ) {
		if ( verbosity >= 10 ) {
			Rcout << "making search path" << endl;
		}
		Rcpp::NumericVector stepsV( history_steps.size() );

		Rcpp::StringVector ngramsV( history_steps.size() );

		for( unsigned int i = 0; i < history_steps.size(); i++ ) {
			stepsV[i] = history_steps[i];
			ngramsV[i] = history_ngrams[i];
		}

		Rcpp::List rulesV( 2 );
		rulesV[0] = stepsV;
		rulesV[1] = ngramsV;
		if ( verbosity >= 10 ) {
					Rcout << "finished making search path" << endl;
		}

		return rulesV;
	}



	/**
	 * Make Rcpp List of all the rules
	 */
	Rcpp::List make_rule_set(std::map<string,LONG_DOUBLE> &model ) {
		Rcpp::List rulesV( model.size() );

		int cntr = 0;
		for(std::map<string,LONG_DOUBLE>::iterator iter = model.begin(); iter != model.end(); ++iter) {
			// [Fix] my_assert(  used_rule_cache.count( iter->first )  > 0 );
			space_t *rule = used_rule_cache[iter->first];

			rulesV[ cntr ] = rule->to_Rcpp_list();
			cntr++;
		}
		//		mylist.names() = CharacterVector::create("a","b");

		return rulesV;
	}


	/**
	 * Make printout of the selected rules with their correponding Z
	 * frequency counts, etc.
	 *
	 * (This could probably be pushed into R instead of being in C++.)
	 */
	DataFrame make_full_model_dataframe( std::map<string,LONG_DOUBLE> &model ) {

		Rcpp::NumericVector betaV;
		Rcpp::NumericVector ZV;
		Rcpp::StringVector ngramV;
		Rcpp::IntegerVector posCountV, negCountV, numPosV, numNegV, supportV, totalDocV;
		//, posWordCountV, negWordCountV;

		for(std::map<string,LONG_DOUBLE>::iterator iter = model.begin(); iter != model.end(); ++iter) {
					// [Fix] my_assert( used_rule_cache.count( iter->first )  > 0 ); // "make_full_model_dataframe number of rules is more than 0",
					space_t *rule = used_rule_cache[iter->first];
					LONG_DOUBLE beta = iter->second;

					int pos_count = 0;
				//	int pos_word_count = 0;
					int neg_count = 0;
				//	int neg_word_count = 0;
					for ( unsigned int i = 0; i < rule->doc_support.size(); i++ ) {
							unsigned int dc = rule->doc_support[i];
							if ( y[dc] == 1 ) {
									pos_count++;
//									pos_word_count += rule->weight[i];
								} else {
									neg_count++;
//									neg_word_count += rule->weight[i];
								}
					}

//					float per = (float)pos_count / (rule->total_docs);
//					float perpos = (float)pos_count / num_pos;
//					float perneg = (float)neg_count / num_neg;
				//	float word_per = (float)pos_word_count / (rule->total_support);


					ngramV.push_back( rule->ngram );
					betaV.push_back( beta );
					ZV.push_back( rule->Z );
					supportV.push_back( rule->total_support );
					totalDocV.push_back( rule->total_docs );
					posCountV.push_back( pos_count );
					negCountV.push_back( neg_count );
//					numPosV.push_back( num_pos );
//					numNegV.push_back( num_neg );
				//	posWordCountV.push_back( pos_word_count );
				//	negWordCountV.push_back( neg_word_count );
		}

		return DataFrame::create(_["ngram"]=ngramV, _["beta"]= betaV, _["Z"]= ZV,
				_["support"] = supportV, _["totalDocs"] = totalDocV,
				//_["numPos"]=numPosV, _["numNeg"]=numNegV,
				_["posCount"]=posCountV, _["negCount"]=negCountV );
				//_["posWordCount"] = posWordCountV, _["negWordCount"] = negWordCountV );
	}



	void print_out_model( std::map<string,LONG_DOUBLE> &model, std::ostream &os ) {
		for(std::map<string,LONG_DOUBLE>::iterator iter = model.begin(); iter != model.end(); ++iter) {
			os << iter->second << '\t' << iter->first << endl;
		}
	}





	void print_ban_list(std::ostream &os) {
		for(std::set<string>::iterator iter = banned_words.begin();
            iter != banned_words.end(); ++iter) {
			os << *iter << "\n";
		}
	}






	/*
	 * Add a document to the corpus if label is non-zero
	 */
    void add_document ( std::string doc, int label ) {
        if ( verbosity > 9 ) {
        	Rcout << "Adding document w/ label = " << label << endl;
        }
        if ( label != 0 ) {
            corpus.push_back( doc );
            y.push_back( label );
        }
	}

    void add_banned_word ( std::string word ) {
        banned_words.insert( word );
	}


	/*
	 * File format: Each line is just the corpus text.
	 * with corresponding line in label vector of 1 or -1 or 0
     *
     */
	bool read_in_data (const char *filename, Rcpp::NumericVector &labelV ) {
		// Set the max line/document size to (10Mb).
		const int kMaxLineSize = 10000000;
		char *line = new char[kMaxLineSize];
		string doc;

		num_pos = 0;
		num_neg = 0;

		if ( verbosity > 0 ) {
			Rcout << "\nLoading the data file from '" << filename << "'\n";
		}

		gettimeofday(&t_origin, NULL);
		std::ifstream ifs (filename);
		if (! ifs) return false;

		int line_cntr = 0;
		while (ifs.getline (line, kMaxLineSize)) {
			// check for comments (but keep empty lines)
			//	 if (line[0] == '\0' || line[0] == ';') continue;
			int len = strlen(line);
			if ( len > 0 && line[0] == ';') continue;
			if ( len > 0 && line[strlen(line) - 1 ] == '\r') {
				line[strlen(line) - 1 ] = '\0';
				len--;
			}
			if ( len == 0 ) {
				Rcout << "WARNING: empty line on line " << (1+line_cntr) << endl;
			}

			// Prepare doc. Assumes 'line' is the original text, e.g. no bracketing of original doc.
			doc.assign(line);
			add_document( doc, labelV[line_cntr] );

			//Rcout << "\nscanning docid: " << transaction.size() - 1 << " class y:" << _y << " doc:" << doc <<"*";
			Rcout.flush();

			line_cntr++;
		}

		delete [] line;

		ifs.close();

		if ( verbosity > 0 ) {
			Rcout << "Read Number of documents = " << corpus.size() << endl;
			gettimeofday(&t, NULL);
			Rcout << "( " << (t.tv_sec - t_origin.tv_sec) << " seconds; " << (t.tv_sec - t_origin.tv_sec) / 60.0 << " minutes )\n";
			Rcout.flush();
		}

		// [Fix] my_assert( labelV.size() == line_cntr ); //"right number of labels",
		return true;
	}




    void finish_initializing() {
		num_pos = 0;
        num_neg = 0;
        // [Fix] my_assert( y.size() == corpus.size() ); // "label vec match corp size",

        for ( unsigned int i = 0; i < y.size(); i++ ) {
            if (y[i] == 1) num_pos++ ;
            if (y[i] == -1) num_neg++ ;
        }

        if ( verbosity > 0 ) {
        	Rcout << "There are " << banned_words.size() << " banned words\n";

        	Rcout << "Number of documents = " << corpus.size() << endl;

        	Rcout << "\n# positive samples: " << num_pos;
        	Rcout << "\n# negative samples: " << num_neg << "\n";

        	gettimeofday(&t, NULL);
        	Rcout << "( " << (t.tv_sec - t_origin.tv_sec) << " seconds; " << (t.tv_sec - t_origin.tv_sec) / 60.0 << " minutes )\n";
        	Rcout.flush();
        }

		// [Fix] my_assert(  objective == 0 || objective == 1 || objective == 2 ); //"has valid objective",


    }




	LONG_DOUBLE calc_exp_fraction( LONG_DOUBLE y, LONG_DOUBLE xbeta ) {
		if (y * xbeta > 8000) {
			return 0;
		} else {
			return  1.0 / (1 + exp(y * xbeta ));
		}

	}


	// compute gradient of given featue 'space' at xbeta_0 without including the regularization terms
	// loc is list of documents that have a specific feature
	space_t *calc_gradient( space_t *space, std::vector<LONG_DOUBLE> &xbeta_0  ) {
		LONG_DOUBLE dbeta = 0.0;

		if ( verbosity > 4 ) {
			Rcout << "Calc gradient of " << space->ngram << endl;
		}
		// [Fix] my_assert( space->isConverted() ); // "calc gradient using converted space",

		vector<int> &weights = space->weight;
		std::vector <unsigned int>& support = space->doc_support;

		space->gradient = 0;
		space->upos = 0;
		space->uneg= 0;

		// Compute the gradient and the upper bound on gradient of extensions.
		for (unsigned int i = 0; i < support.size(); ++i) {
			unsigned int j = support[i];

			// Choose objective function. 0: SLR, 1: SVM.
			if (objective == 0) { //SLR
				// From differentiation we get a - in front of the sum_i_to_N
				dbeta = y[j] * exp_fraction[j] * weights[i];
				//				if ( y[j] > 0 ) {
				//					grad.upos -= y[j] * exp_fraction[j] * weights[i];
				//				} else {
				//					grad.uneg -= y[j] * exp_fraction[j] * weights[i];
				//				}
			} else if (objective == 1) { //L1-SVM (Hinge loss)
				if (1 - y[j] * xbeta_0[j] > 0) {
					dbeta =  y[j] * weights[i];
					//					if (y[j] > 0) {
					//						grad.upos -= y[j] * weights[i];
					//					} else {
					//						grad.uneg -= y[j] * weights[i];
					//					}
				} else {
					dbeta = 0;
				}
			} else if (objective == 2) { //L2-SVM
				if (1 - y[j] * xbeta_0[j] > 0) {
					dbeta = 2 * y[j] * (1 - y[j] * xbeta_0[j]) * weights[i];
					//					if (y[j] > 0) {
					//						grad.upos -= 2 * y[j] * (1 - y[j] * xbeta_0[j]) * weights[i];
					//					} else {
					//						grad.uneg -= 2 * y[j] * (1 - y[j] * xbeta_0[j]) * weights[i];
					//					}
				} else {
					dbeta = 0;
				}
			}

			if ( dbeta != 0 ) {
				// add to gradient
				if ( y[j] == 1 && positive_weight != 1 ) {
					dbeta *= positive_weight;
				}
				space->gradient -= dbeta;

				if ( Lp== 1 ) {
					if ( y[j] > 0 ) {
						space->upos = max( space->upos, dbeta );
					} else {
						space->uneg = max( space->uneg, -dbeta );
					}
				} else if ( Lp >= LP_INFINITY ) {
					if ( y[j] > 0 ) {
						space->upos += dbeta;
					} else {
						space->uneg += dbeta;
					}
				} else{
					if ( y[j] > 0 ) {
						space->upos += pow( dbeta, Lq );
					} else {
						space->uneg += pow( -dbeta, Lq );
					}
				}
			} // end upos/uneg update
		}

		if ( verbosity > 3 ) {
			Rcout << "found mass: " << space->gradient << endl;
		}

		// normalize
		space->gradient /= space->Z;
		if ( Lp > 1 && Lp < LP_INFINITY ) {
			space->upos = pow( space->upos, 1/Lq );
			space->uneg = pow( space->uneg, 1/Lq );
		}
		if ( verbosity > 3 ) {
			space->print_rule( );
			Rcout << "\tGradient (w/o reg) = " << space->gradient << " : " << space->uneg << " (neg) / " << space->upos << " (pos)" << endl;
		}
		return space;
	}



	// compute gradient for intercept
	LONG_DOUBLE calc_int_gradient( std::vector<LONG_DOUBLE> &xbeta_vector ) {
		LONG_DOUBLE grad = 0.0;
		LONG_DOUBLE sc = 0.0;
		//Rcout << "int grad calc" << endl;
		for (unsigned int j = 0; j < corpus.size(); ++j) {
			sc =  ( y[j] == 1 ) ? positive_weight : -1;
			if (objective == 0) { //SLR
				// From differentiation we get a - in front of the sum_i_to_N
				grad -= sc * calc_exp_fraction( y[j], xbeta_vector[j] );
			} //end SLR

			if (objective == 1) { //L1-SVM (Hinge loss)
				if (1 - y[j] * xbeta_vector[j] > 0) {
					grad -= sc;
				}
			} //end SVM

			if (objective == 2) { //L2-SVM
				if (1 - y[j] * xbeta_vector[j] > 0) {
					grad -= 2 * sc * (1 - y[j] * xbeta_vector[j]);
					//Rcout << j << "\t" << y[j] << "\t" << xbeta_vector[j] << "\t" << (-2 * y[j] * (1 - y[j] * xbeta_vector[j])) << "\t" << grad << endl;
				}
			} //end SVM

		}
		return grad;
	}



	void print_diagnostics_for_prune( space_t *space ) {
		Rcout << "\ncurrent ngram rule: '" << space->ngram << "'";
		Rcout << "\n\tlocation size: " << space->total_support;
		//for (unsigned int i = 0; i < space->loc.size(); ++i)
		//	Rcout <<  space->loc[i] << " ";
		Rcout << "\n\tgradient (before regularizer): " << space->gradient;
		Rcout << "\n\tupos (before regularizer): " << space->upos;
		Rcout << "\n\tuneg (before regularizer): " << space->uneg;
		Rcout << "\n\ttau: " << tau << endl;
	}


	void update_best_rule( space_t *space ) {
		if (verbosity >= 3) {
			Rcout << "updating current best ngram rule to " << space->ngram << " \tgradient: " << space->gradient << "\n";
		}

		rewritten++;

		tau = std::abs(space->gradient) + TINY_WEIGHT;
		rule = space;
	}

	void update_best_rule_to_intercept( LONG_DOUBLE gradient ) {
		tau = std::abs(gradient) + TINY_WEIGHT;
		int_rule.gradient = gradient;
		rule = &int_rule;
	}



	int calc_sign( LONG_DOUBLE val ) {
		return val > 0 ? +1 : (val < 0 ? -1 : 0);
		//		abs(features_it->second)/features_it->second
	}

	/*
	 * Shrink val amount to 0, but stick at 0 if needed.
	 */
	LONG_DOUBLE shrink_to_zero( LONG_DOUBLE val, LONG_DOUBLE amount ) {
		if ( val > amount ) {
			return val - amount;
		} else if ( val < -amount ) {
			return val + amount;
		} else {
			return 0;
		}
	}


	LONG_DOUBLE penalize( LONG_DOUBLE grad, LONG_DOUBLE beta, LONG_DOUBLE alpha, LONG_DOUBLE C ) {
		if ( beta == 0 ) {
			return shrink_to_zero( grad, alpha*C );
		} else if ( beta < 0 ) {
			return grad - C*alpha + (1-alpha)*beta;
		} else {
			return grad + C*alpha + (1-alpha)*beta;
		}
	}




	/*
	 * For current ngram, compute the gradient value and update the current optimal ngram
	 * if current one ('space') is better.
	 *
	 * cur_beta = features_cache.end() means space is still at 0 for its coefficient.
	 */
	void calc_gradient_and_update( space_t *space, map<string, LONG_DOUBLE>::iterator cur_beta ) {
		total++;

		// [Fix] my_assert( space->isConverted() );// "calc-and-update using converted space",
		// [Fix] my_assert( space->total_support >= minsup );

		if ( verbosity > 4 ) {
			space->print_space();
		}

		// calc gradient for current ngram, _not_ including regularization
		calc_gradient(space, sum_best_xbeta );

		// now include regularization, if needed
		if ( C != 0 ) {
			if ( cur_beta != features_cache.end() && cur_beta->second != 0 ) {
				// update the gradient: gradient += C * [alpha*sign(beta_j) + (1-alpha)*beta_j];
				//Rcout << "penalizing gradient: " << space->gradient << "\t beta=" << cur_beta->second << "\t alpha=" << alpha << "\t C=" << C << endl;
				space->gradient = penalize( space->gradient, cur_beta->second, alpha, C );
			} else {
				// we are at 0, so shrink gradient to 0 by the L1.
				//Rcout << "shrinking gradient: " << space->gradient << "\t" << alpha << "\t" << C << endl;
				space->gradient = shrink_to_zero( space->gradient, alpha*C );
			}

			if (verbosity > 3) {
				Rcout << "gradient after regularizer: " << space->gradient << endl;
			}
			if ( verbosity > 5) {
				print_diagnostics_for_prune( space );
			}
		}

		if ( (space->size >= minpat) && (std::abs (space->gradient) > tau)) {
			if ( !pos_only || space->gradient < 0 || (cur_beta != features_cache.end() ) ) {
				update_best_rule( space );
			}
		}
	}

	void calc_gradient_and_update( space_t *space ) {
		calc_gradient_and_update( space, features_cache.end() );
	}



	/**
	 * Go through current model and calculate gradient for each non-zero feature
	 * and update current best gradient.
	 *
	 * NOTE: After this, the BFS or DFS can assume everything is at 0, since we got all of our
	 *   nonzero spots.
	 */
	void check_model_features() {
		map<string, LONG_DOUBLE>::iterator feat_it = features_cache.begin();
		for ( feat_it = features_cache.begin(); feat_it != features_cache.end(); feat_it++ ) {
			// [Fix] my_assert( used_rule_cache.count( feat_it->first )  > 0 );
			space_t *space = used_rule_cache[ feat_it->first ];

			calc_gradient_and_update( space, feat_it );
		}
	}


	// For current ngram, check pruning conditions for its AT-ZERO children
	// (any children already in model will be automatically checked at start of cycle)
	bool can_prune (space_t *space ) {

		// Assume we can prune all the children until bound says otherwise.
		bool flag_can_prune = 1;

		if ( std::max (space->upos, space->uneg) - C*alpha > tau ) {
			flag_can_prune = 0;
		}

		if ( flag_can_prune ) {
			++pruned;
			if (verbosity > 3) {
				Rcout << "Pruned all at-zero children of " << space->ngram << "!\n";
			}
			return true;
		} else {
			return false;
		}
	}





	/**
	 * This method expands a ngram and discovers all its children.
	 *
	 * Note: other methods dumps the info to calculate this in order to save space.
	 * When this is done, this is called a shrunk space.
	 *
	 * Given a NON-SHRUNK space 'space' generate a map of all possible
	 * following unigrams (and, possibly '*') for the n-gram.
	 *
	 * Do this by iterating across all locations of where the ngram 'space'
	 * is in the doc database and finding the subsequent unigrams.
	 */
	std::map<string, space_t> find_candidates( space_t *space ) {
		//Rcout << "\nnext.empty";
		// Prepare possible extensions.
		unsigned int docid = 0;
		// Candidates obtained by extension.
		std::map<string, space_t> candidates;
		std::vector <int>& loc = space->get_loc();

		if ( verbosity > 4 ) {
			Rcout << "\nFinding candidates for " << space->ngram << endl;
		}

		// Iterate through the inverted index of the current feature.
		for (unsigned int i = 0; i < loc.size(); ++i) {

			// check if starting new document.  Stash document id if so.
			if (loc[i] < 0) {
				docid = (unsigned int)(-loc[i]) - 1;
				//Rcout << "\ndocid: " << docid << " with size " << corpus[docid].size();
				continue;
			}

			// the start position of this unigram in the document
			unsigned int pos = loc[i];

			string next_unigram;
			// If not re-initialized, it should fail.
			unsigned int pos_next_unigram = 0;

			//Rcout << "\npos: " << pos;
			if (pos < corpus[docid].size() - 1) {
				//Rcout << "\npotential extension...";
				if (token_type==WORD_TOKEN) {
					// Word level token.
					size_t pos_next_space = corpus[docid].find(' ', pos + 1);
					//Rcout << "next space: " << pos_next_space << endl;
					// big numbers for find means didn't find, I believe??
					if (  pos_next_space <= corpus[docid].size() - 2) {
						//unigram.assign(transaction[docid].substr(pos, pos_next_space - pos));
						size_t pos_next_next_space = corpus[docid].find(' ', pos_next_space + 1);
						//Rcout << "next next space: " << pos_next_next_space << endl;
						if (pos_next_next_space < string::npos ) {
							next_unigram.assign(corpus[docid].substr(pos_next_space + 1, pos_next_next_space - pos_next_space - 1));
							pos_next_unigram = pos_next_space + 1;
							//Rcout << " next unigram" << next_unigram << " at " << pos_next_unigram << endl;
						} else {
							next_unigram.assign( corpus[docid].substr(pos_next_space + 1, corpus[docid].size() - 1) );
							pos_next_unigram = pos_next_space + 1;
							//Rcout << " next unigram (end string) " << next_unigram << " at " << pos_next_unigram << endl;
						}
					} else continue;
				} else { // CHARACTER TOKENS
					//unigram = transaction[docid][pos];
					// Char (i.e. byte) level token. Skip spaces.
					if (!isspace(corpus[docid][pos + 1])) {
						//Rcout << "\nnext char is not space";
						next_unigram = corpus[docid][pos + 1];
						pos_next_unigram = pos + 1;
					} else {
						if (pos + 2 <= corpus[docid].size() - 1) {
							if (isspace(corpus[docid][pos + 2])) {
								Rcerr <<"\nFATAL...consecutive spaces...exit";
								::Rf_error( "Consecutive spaces in passed corpus.  Remove them and rerun." );
							}
							//Rcout << "\nnext next char is not space";
							next_unigram = corpus[docid][pos + 2];
							pos_next_unigram = pos + 2;
						} else continue;
					}
				}

				// if token is not banned
				if ( banned_words.count( next_unigram ) == 0 ) {

					// make sure the token itself is not too rare, since if it were then the multi-token expression would be even more rare.
					if (minsup == 1 || single_node_minsup_cache.find (next_unigram) != single_node_minsup_cache.end()) {
						// next line can potentially add a new candidate, or append an existing one.
						//Rcout << "candidate added: +'" << next_unigram << "' doc#" << docid << " - pos#" << pos_next_unigram << endl;
						candidates[next_unigram].add (docid, pos_next_unigram);
					}

					if (maxgap && space->num_consec_gaps() < maxgap ) {
						// If we allow gaps, we treat a gap as an additional unigram "*".
						// Its inverted index will simply be a copy pf the inverted index of the original features.
						// Example, for original feature "a", we extend to "a *", where inverted index of "a *" is simply
						// a copy of the inverted index of "a", except for positions where "a" is the last char in the doc.
						candidates["*"].add (docid, pos_next_unigram);
					}
				}
			}
		}
		return candidates;
	} // end find_candidates



	/*
	 * Drop the locations of the ngram 'space' and add all of its ngram children
	 * as leaves to its tree, if those children have enough support.
	 *
	 * This keeps all the extensions of the current feature for future iterations.
	 * If we need to save memory we could sacrifice this storage.
	 */
	void extend_space_tree( space_t *space, std::map<string, space_t> &candidates ) {
		// Keep only doc_ids for occurrences of current ngram.
		space->shrink ();
		std::vector <space_t *>& children = space->children;

		for (std::map <string, space_t >::iterator it = candidates.begin();
             it != candidates.end(); ++it) {

			if ( it->second.total_support >= minsup ) {
				space_t* c = new space_t( it->second, space, it->first );
				total_nodes++;

				children.push_back (c);
			}
		}

		if (children.empty()) {
			children.push_back (0);
		}
		std::vector<space_t *>(children).swap (children);
	}



	void check_child( space_t *child, std::vector<space_t *>& new_space ) {

		if ( !child->isConverted() ) {
			if ( verbosity > 3 ) {
				Rcout << "Converting space and calculating support and weights for " << child << "\n";
			}
			child->calc_support_weights( Lp, binary_features, no_regularization );
		}

		if ( child->total_support < minsup ) {
			return;
		}

		// check if we already have beta.
		map<string, LONG_DOUBLE>::iterator features_it = features_cache.find(child->ngram);
		if (features_it != features_cache.end()) {
			// we have calculated at beginning.  Gradient info already done.
			//Rcout << "skipping calulating prev. term " << child->ngram << endl;
		} else if ( child->prev != NULL && child->total_support == child->prev->total_support ) {
			// we have same support as our parent so we must have same gradient information too.
			// no need to recalculate
			//Rcout << "copy over gradient for " << child->ngram << endl;
			child->gradient = child->prev->gradient;
			child->upos = child->prev->upos;
			child->uneg = child->prev->uneg;
		} else {
			// first check to see if we are the new best bet.
			calc_gradient_and_update( child, features_it );
		}

		// Now if our children have potential, enter ourself in to next round of expansions.
		// We will do the gradient check for pruning later.
		if ( child->size < maxpat ) {
			new_space.push_back (child);
		}
	}



	// 'space' is the parent ngram which has already been checked for being best, and so forth.
    // We are supposed to be examining all of 'space's children as one part of the bfs level traversal.
	// First step is to see if we can cut all of 'space's children and skip a bunch of work.
	// If not, for each child, calc gradient and update best rule.  Also add them to the next level
	// (We are being lazy and calc'ing if we can prune as late as possible so the current best gradient can be the biggest)
	// The growth is done breadth-first rather than depth-first, e.g. grow all bi-grams, than all tri-grams, etc.
	void span_bfs (space_t *space, std::vector<space_t *>& new_space, unsigned int size) {

		// can we prune all our children right away?
		if ( can_prune( space ) ) {
			return;
		}

		std::vector <space_t *>& children = space->children;

		// MEMORY VERSION
		if ( children.empty() ) {
			// compute all possible candidate extensions.  These are our children.
			std::map<string, space_t> candidates = find_candidates( space );
			extend_space_tree( space, candidates );
		}

		if (children.size() == 1 && children[0] == 0) {
			// we have extended before, and found no viable children.  Bail.
			return;
		}

		// Iterate through all our children, see if any are good.  Enter them into the next level
		// of the bfs if they have potential.
		for (std::vector<space_t*>::iterator it = children.begin(); it != children.end(); ++it) {
			check_child( (*it), new_space );
		}

		// Recalc every time version
		// compute all possible candidate extensions
		/*		std::map<string, space_t> candidates = find_candidates( space );
         for (std::map <string, space_t >::iterator it = candidates.begin();	it != candidates.end(); ++it) {

         if ( check_child( it->second, new_space, size ) ) {

         space_t* c = new space_t( space, it->first );
         total_nodes++;
         c->loc = it->second.loc;
         std::vector<int>(c->loc).swap(c->loc);

         children.push_back (c);

         //Rcout << "\nnode extended: " << space->ne;
         //Rcout << "\nnode for extention: " << it->first;
         }

         if (children.empty()) {
         children.push_back (0);
         }
         std::vector<space_t *>(children).swap (children);
         }
         for (std::vector<space_t*>::iterator it = children.begin(); it != children.end(); ++it) {

         }
         }
         */
	} //end span_bfs







	// Try to grow the ngram to next level, and prune the appropriate extensions.
	// The growth is done deapth-first rather than breadth-first, e.g. grow each candidate to its longest unpruned sequence
	void span_dfs (space_t *space, unsigned int size) {

		Rcerr << "WARNING: span_dfs is NOT UPDATED FUNCTION" << endl;

	}


	/**
	 * Calc penalty regularazation term.
	 * If new_beta_coef and old_beta_coef are not 0, then it will add new and subtract old from the sums of betas
	 * to get appropriate penalty.
	 */
	LONG_DOUBLE calc_penalty( LONG_DOUBLE sum_abs_betas, LONG_DOUBLE sum_squared_betas, LONG_DOUBLE new_beta_coef, LONG_DOUBLE old_beta_coef ) {
		if (old_beta_coef == 0 ) {
			return C * (alpha * (sum_abs_betas + abs(new_beta_coef)) + (1 - alpha) * 0.5 * (sum_squared_betas + pow(new_beta_coef, 2)));
		} else {
			return C * (alpha * (sum_abs_betas - abs(old_beta_coef) + abs(new_beta_coef))
                        + (1 - alpha) * 0.5 * (sum_squared_betas - pow(old_beta_coef, 2) + pow(new_beta_coef, 2)));

		}
	}

	/* calculates loss given (1) the prediction vector of X'b and
     *	(2) the currently-being-updated feature as indicated by features_it and delta_beta.
     *
     *	Prediction is X'b
     * note: if features_it is at end, and delta_beta != 0, then beta_coef is for a new feature.
     *    We don't need to know which feature corresponds to delta_beta in this case to calc loss given the X'b.
     *    delta_beta captures change in beta for feature (included in X'b calc already).   old value can be retrieved by features_it
     */
	LONG_DOUBLE calc_loss( 	std::vector<LONG_DOUBLE> &sum_betas,
                          map<string, LONG_DOUBLE>::iterator &features_it,
                          LONG_DOUBLE delta_beta,
                          bool penalize ) {
		LONG_DOUBLE loss = 0;
		LONG_DOUBLE delt = 0.0;
		for (unsigned int i = 0; i < corpus.size();  ++i) {
			// Compute loss.
			delt = 0;
			if (objective == 0) { //SLR
				exp_fraction[i] = calc_exp_fraction( y[i], sum_betas[i] );

				if (-y[i] * sum_betas[i] > 8000) {
					delt = log(LDBL_MAX);
				} else {
					delt = log(1 + exp(-y[i] * sum_betas[i]));
				}
			} //end SLR
			if (objective == 1) { //L1-SVM (Hinge loss)
				if (1 - y[i] * sum_betas[i] > 0)
					delt =  (1 - y[i] * sum_betas[i]);
			} //end L1-SVM
			if (objective == 2) { //L2-SVM
				if ( 1 > y[i] * sum_betas[i] )
					delt =  pow(1 - y[i] * sum_betas[i], 2);
			} //end L2-SVM
			//Rcout << "\n y[i]: " << y[i];
			//Rcout << "\nsum_betas[i]: " << sum_betas[i];
			//Rcout << "\nloss: " << loss;
	        if ( y[i] == 1 && positive_weight != 1 ) {
	        	loss += positive_weight * delt;
	        } else {
	        	loss += delt;
	        }
		}

		// Update the log-likelihood with the Elastic Net regularization term.
		if ( penalize && C != 0 ) {
			if (sum_squared_betas == 0 || features_it == features_cache.end() ) {
				// This is the first ngram selected or is a new ngram
				loss = loss + calc_penalty( sum_abs_betas, sum_squared_betas, delta_beta, 0 );
			} else {
				// This feature was selected before.  Adjust to take old value into account
				LONG_DOUBLE new_beta_coef =  features_it->second + delta_beta;
				loss = loss + calc_penalty( sum_abs_betas, sum_squared_betas, new_beta_coef, features_it->second );
			}
		}// end handling possible regularization.

		return loss;
	}


	/**
	 * No longer used.
	 * Don't understand dependence on n --- summing over all docs.  Also, why not just look at beta itself?
	 */
	LONG_DOUBLE calc_current_range( vector<LONG_DOUBLE>& sum_best_beta_n0, vector<LONG_DOUBLE>& sum_best_beta_n2 ) {
		LONG_DOUBLE current_range_size = 0.0;
		for (unsigned int i = 0; i < corpus.size();  ++i) {
			if (step_verbosity > 5 && i == 0) {
				Rcout << "\nsum_best_beta_n0[i]: " << sum_best_beta_n0[i];
				Rcout << "\nsum_best_beta_n2[i]: " << sum_best_beta_n2[i];
			}
			current_range_size += abs(sum_best_beta_n2[i] - sum_best_beta_n0[i]);
		}
		if (step_verbosity > 5) {
			Rcout << "\ncurrent range size |beta_n2 - beta_n0|: " << current_range_size << endl;
		}
		return current_range_size;
	}


	// Calculate the CHANGE in beta for a specific feature r given the xbeta vector
	//
	// Consider two beta vectors beta^a and beta^b that differ only on rule r.  I.e.
	// beta^a_k = beta^b_k for k \neq r, and beta^a_r \neq beta^b_r
	// Then since, for doc i,
	// y_i^a = xbeta^a[i] = mu + sum_{j != r} x_{ij}*beta^a_j + x_{ir}*beta^a_r
	// beta^a_r - beta^b_r = (y_i^b - y_i^a)/x_{ir}  for any i s.t. x_{ir} != 0
	//
	// xbeta_00 is the starting point xbeta for the first doc with nonzero weight on
	// the feature.  I.e. xbeta_00 = xbeta^b[rule.doc_support[0]] for some xbeta^b
	//
	// Note: x_{ir} = weight_{ir} / Z_r, with Z the normalizing constant.
	LONG_DOUBLE extract_delta_beta( space_t &rule, vector<LONG_DOUBLE> &xbeta, LONG_DOUBLE xbeta_00 ) {
		return rule.Z * (xbeta[rule.doc_support[0]] - xbeta_00) / rule.weight[0];
		// incorp Zs
	}

	// Search for the beta in the range beta_n-1, beta_mid_n-1_n, beta_n, beta_mid_n_n+1, beta_n+1
	// that minimizes the objective function. It suffices to compare the 3 points beta_mid_n-1_n, beta_n, beta_mid_n_n+1,
	// as the min cannot be achieved at the extreme points of the range.
	// Take the 3 point range containing the point that achieves minimum loss.
	// Repeat until the 3 point range is too small, or a fixed number of iterations is achieved.
	//
	// This finds the minimal loss point for a given rule
	// It is a line search method. Search for step size that minimizes loss.
	// Compute loss for: (1) beta_n1, the middle point of the range, and
	// also (2), midpoint of (beta_n0, beta_n1) and (3), midpoint for (beta_n1, beta_n2)
	// Compare the loss for the 3 points, and choose resulting range of 3 points
	// which contains the minimum. Repeat until the range spanned by the 3 points is small enough,
	// e.g. the range approximates well the vector where the loss function is minimized.
	// Stash this final middle point of the final best range in 'sum_best_beta_opt'
	//
    // Side effect: sum_best_beta_opt get's set to the new best place after model update.
    //
	// param:
	//  xbeta_00:  the original starting point value of x*beta for the first document in the support of the rule
	//             both above used to back-calculate the step-size and changes to the loss function if regularized.
	void find_best_range(	vector<LONG_DOUBLE>& sum_best_beta_n0,
                         vector<LONG_DOUBLE>& sum_best_beta_n1,
                         vector<LONG_DOUBLE>& sum_best_beta_n2,
                         space_t& rule,
                         LONG_DOUBLE &xbeta_00,
                         vector<LONG_DOUBLE>* sum_best_beta_opt,
                         bool penalize ) {

		vector<LONG_DOUBLE> sum_best_beta_mid_n0_n1(corpus.size());
		vector<LONG_DOUBLE> sum_best_beta_mid_n1_n2(corpus.size());

		LONG_DOUBLE min_range_size = 1e-4;
		//LONG_DOUBLE current_range_size = 0;
		int current_interpolation_iter = 0;

		LONG_DOUBLE loss_mid_n0_n1 = 0;
		LONG_DOUBLE loss_mid_n1_n2 = 0;
		LONG_DOUBLE loss_n1 = 0;

		//current_range_size = calc_current_range( sum_best_beta_n0, sum_best_beta_n2 );

		LONG_DOUBLE beta_coef_n0 = extract_delta_beta( rule, sum_best_beta_n0, xbeta_00 );
		LONG_DOUBLE beta_coef_n1 = 0;
		LONG_DOUBLE beta_coef_mid_n0_n1 = 0;
		LONG_DOUBLE beta_coef_mid_n1_n2 = 0;
		LONG_DOUBLE beta_coef_n2 = extract_delta_beta( rule, sum_best_beta_n2, xbeta_00 );

		map<string, LONG_DOUBLE>::iterator features_it = features_cache.find(rule.ngram);

		if ( step_verbosity > 3 ) {
			Rcout << "beta_coef n0 = " << beta_coef_n0 << "\nbeta_coef_n2 = " << beta_coef_n2 << endl;
		}

		// Start interpolation loop.
		while (abs( beta_coef_n2 - beta_coef_n0 ) > min_range_size ) {

			if (step_verbosity > 3) {
				Rcout << "\nCurrent interpolation iteration: " << current_interpolation_iter;
			}

			for (unsigned int i = 0; i < corpus.size();  ++i) { //loop through training samples
				sum_best_beta_mid_n0_n1[i] = (sum_best_beta_n0[i] + sum_best_beta_n1[i]) / 2;
				sum_best_beta_mid_n1_n2[i] = (sum_best_beta_n1[i] + sum_best_beta_n2[i]) / 2;
			}
			if (step_verbosity > 5) {
				Rcout << "\nsum_best_beta_mid_n0_n1[0]: " << sum_best_beta_mid_n0_n1[0];
				Rcout << "\nsum_best_beta_mid_n1_n2[0]: " << sum_best_beta_mid_n1_n2[0];
			}

			if ( C != 0) {
				beta_coef_n1 = extract_delta_beta( rule, sum_best_beta_n1, xbeta_00 );
				beta_coef_mid_n0_n1 = extract_delta_beta( rule, sum_best_beta_mid_n0_n1, xbeta_00 );
				beta_coef_mid_n1_n2 = extract_delta_beta( rule, sum_best_beta_mid_n1_n2, xbeta_00 );
			}

			loss_n1 = calc_loss( sum_best_beta_n1, features_it, beta_coef_n1, penalize );
			loss_mid_n0_n1 = calc_loss( sum_best_beta_mid_n0_n1, features_it, beta_coef_mid_n0_n1, penalize );
			loss_mid_n1_n2 = calc_loss( sum_best_beta_mid_n1_n2, features_it, beta_coef_mid_n1_n2, penalize );


			// Focus on the range that contains the minimum of the loss function.
			// Compare the 3 points beta_n, and mid_beta_n-1_n and mid_beta_n_n+1.
			if (loss_n1 <= loss_mid_n0_n1 && loss_n1 <= loss_mid_n1_n2) {
				// Min is beta_n1.
				if (step_verbosity > 5) {
					Rcout << "\nmin is sum_best_beta_n1";
					Rcout << "\n\tloss_mid_n0_n1: " << loss_mid_n0_n1 << "\tloss_n1: " << loss_n1 << "\tloss_mid_n1_n2: " << loss_mid_n1_n2;
				}
				// Make the beta_n0 be the beta_mid_n0_n1 and the beta_n2 be the beta_mid_n1_n2
				sum_best_beta_n0.assign(sum_best_beta_mid_n0_n1.begin(), sum_best_beta_mid_n0_n1.end());
				sum_best_beta_n2.assign(sum_best_beta_mid_n1_n2.begin(), sum_best_beta_mid_n1_n2.end());
			} else if (loss_mid_n0_n1 <= loss_n1 && loss_mid_n0_n1 <= loss_mid_n1_n2) {
				// Min is beta_mid_n0_n1.
				if (step_verbosity > 5) {
					Rcout << "\nmin is sum_best_beta_mid_n0_n1";
					Rcout << "\n\tloss_mid_n0_n1: " << loss_mid_n0_n1 << "\tloss_n1: " << loss_n1 << "\tloss_mid_n1_n2: " << loss_mid_n1_n2;
				}
				// Make the beta_n2 be the beta_n1 and the beta_n1 be the beta_mid_n0_n1.
				sum_best_beta_n2.assign(sum_best_beta_n1.begin(), sum_best_beta_n1.end());
				sum_best_beta_n1.assign(sum_best_beta_mid_n0_n1.begin(), sum_best_beta_mid_n0_n1.end());
			} else if (loss_mid_n1_n2 <= loss_n1 && loss_mid_n1_n2 <= loss_mid_n0_n1) {
				// Min is beta_mid_n1_n2.
				if (step_verbosity > 5) {
					Rcout << "\nmin is sum_best_beta_mid_n1_n2";
					Rcout << "\n\tloss_mid_n0_n1: " << loss_mid_n0_n1 << "\tloss_n1: " << loss_n1 << "\tloss_mid_n1_n2: " << loss_mid_n1_n2;
				}
				// Make the beta_n0 be the beta_n1 and the beta_n1 be the beta_mid_n1_n2.
				sum_best_beta_n0.assign(sum_best_beta_n1.begin(), sum_best_beta_n1.end());
				sum_best_beta_n1.assign(sum_best_beta_mid_n1_n2.begin(), sum_best_beta_mid_n1_n2.end());
			}

			++current_interpolation_iter;

			// comment: changed so we converge when the beta coef endpoints gets close.
			// the calc_current_range depends on n, and seems weird to me.  Am I missing something?
			beta_coef_n0 = extract_delta_beta( rule, sum_best_beta_n0, xbeta_00 );
			beta_coef_n2 = extract_delta_beta( rule, sum_best_beta_n2, xbeta_00 );
			if ( step_verbosity > 3 ) {
				Rcout << "\nrange = " << beta_coef_n0 << " - " << beta_coef_n2 << " = " << (beta_coef_n0-beta_coef_n2) << endl;
			}
			//current_range_size = calc_current_range( sum_best_beta_n0, sum_best_beta_n2 );
		} // end while loop.


		if ( step_verbosity > 3 ) {
			Rcout << "Finishing find_best_range and assigning midpoint to sum_best_beta_opt" << endl;
			Rcout << "The midpoint of the range found\n";
			print_estimates( sum_best_beta_n1, 1 );
		}

		// Keep the middle point of the best range.
		sum_best_beta_opt->assign(sum_best_beta_n1.begin(), sum_best_beta_n1.end() );

		//sum_best_beta_opt->clear();
		//for (unsigned int i = 0; i < transaction.size();  ++i) {
		//	sum_best_beta_opt->push_back(sum_best_beta_n1[i]);
		// Trust this step only by a fraction.
		//sum_best_beta_opt->push_back(0.5 * sum_best_beta[i] + 0.5 * sum_best_beta_n1[i]);
		//}
		//Rcout << "\n end find_best_range()!";
	} // end find_best_range().




	// looking at moving along feature as defined by rule
	// find three points, n0, n1, n2 with the minimum loss somewhere in that range
	void find_bracketing_range( space_t &rule,
                               vector<LONG_DOUBLE>& sum_best_beta_n0, vector<LONG_DOUBLE>& sum_best_beta_n1, vector<LONG_DOUBLE>& sum_best_beta_n2,
                               bool penalize, bool pos_beta_only ) {
		// Starting value for parameter in step size search.
		// Set the initial epsilon value small enough to guarantee
		// log-likelihood increases in the first steps.
		LONG_DOUBLE exponent = ceil(log10(abs(rule.gradient)));
		LONG_DOUBLE epsilon = min((float)1e-3, pow((float)10, (float)-exponent));
		LONG_DOUBLE step = 0.0;

		if (step_verbosity > 3) {
			Rcout << "\nrule.ngram: " << rule.ngram << "\nrule.gradient: " << rule.gradient << "\nexponent of epsilon: " << -exponent << "\nepsilon: " << epsilon << endl;
		}

		// To Keep track of loss at the three points, n0, n1, n2.
		LONG_DOUBLE loss_n0 = 0;
		LONG_DOUBLE loss_n1 = 0;
		LONG_DOUBLE loss_n2 = 0;

		// Binary search for epsilon. Similar to bracketing phase in which
		// we search for some range with promising epsilon.
		// The second stage finds the epsilon or corresponding weight vector with smallest l2-loss value.

		// As long as the l2-loss decreases, double the epsilon.
		// Keep track of the last three values of beta, or correspondingly
		// the last 3 values for the scalar product of beta and xi.
		int n = 0;

		// find prev incarnation of rule
		map<string, LONG_DOUBLE>::iterator  features_it = features_cache.find(rule.ngram);
		loss_n1 = loss_n2 = calc_loss( sum_best_beta_n1, features_it, 0, true );

		// find current beta for rule to check crossing 0 boundary
		LONG_DOUBLE old_beta_Z = 0.0;
		if ( features_it != features_cache.end() ) {
			old_beta_Z = features_it->second / rule.Z;
			if ( step_verbosity > 1 ) {
				Rcout << "Old value of beta = " << features_it->second << endl;
			}
		}
		bool hit_wall = 0;
		LONG_DOUBLE beta_coeficient_update_Z = 0;   // cumulative update to beta divided by norm constant Z
		do {
			// For each location (e.g. docid), update the score of the documents containing best rule.
			// E.g. update beta^t * xi.
			step = pow((float)2, (float)(n * 1.0)) * epsilon * rule.gradient / rule.Z;
			if ( pos_beta_only && ((old_beta_Z - step - beta_coeficient_update_Z) <= 0) ) {
				Rcout << "Stopped epsilon doubling due to hitting 0." << endl;
				step = old_beta_Z - beta_coeficient_update_Z;
				hit_wall = 1;
			}
			beta_coeficient_update_Z -= step;
			if (step_verbosity > 5) {
				Rcout << "\nepsilon doubling itr # n = " << n << "\t step=" << step << endl;
			}
			for (unsigned int i = 0; i < rule.doc_support.size(); ++i) {
				sum_best_beta_n0[rule.doc_support[i]] = sum_best_beta_n1[rule.doc_support[i]];
				sum_best_beta_n1[rule.doc_support[i]] = sum_best_beta_n2[rule.doc_support[i]];
				sum_best_beta_n2[rule.doc_support[i]] = sum_best_beta_n1[rule.doc_support[i]] - step*rule.weight[i];
			}
			if (step_verbosity > 5 ) {
				Rcout << "sum_best_beta_n0[rule.doc_support[0]: " << sum_best_beta_n0[rule.doc_support[0]] << "\nsum_best_beta_n1[rule.doc_support[0]: " << sum_best_beta_n1[rule.doc_support[0]] << "\nsum_best_beta_n2[rule.doc_support[0]: " << sum_best_beta_n2[rule.doc_support[0]] << endl;
			}

			loss_n0 = loss_n1;
			loss_n1 = loss_n2;
			loss_n2 = calc_loss( sum_best_beta_n2, features_it, beta_coeficient_update_Z * rule.Z, penalize );

			if (step_verbosity > 5 ) {
				Rcout << "loss_n0: " << loss_n0 << "\nloss_n1: " << loss_n1 << "\nloss_n2: " << loss_n2 << endl;
			}
			++n;
		} while (!hit_wall && loss_n2 < loss_n1);


		if (step_verbosity > 3) {
			Rcout << "Finished doubling epsilon! The monotonicity loss_n+1 < loss_n is broken! Reached iter #" << n << "\n";
			Rcout << "\tfinal losses: loss_n0: " << loss_n0 << "\tloss_n1: " << loss_n1 << "\tloss_n2: " << loss_n2 << endl;
			Rcout << "\txbeta_00s: " << sum_best_beta_n0[rule.doc_support[0]] << "\t" << sum_best_beta_n1[rule.doc_support[0]] << "\t" << sum_best_beta_n2[rule.doc_support[0]] << endl;
		}

		// if we immediately passed low point, set n1 to midpoint of n0 and n2 so we have different points
		// (currently n1 = n0)
		if ( n == 1 ) {
			for (unsigned int i = 0; i < rule.doc_support.size(); ++i) {
				sum_best_beta_n1[rule.doc_support[i]] = (sum_best_beta_n0[rule.doc_support[i]] + sum_best_beta_n2[rule.doc_support[i]])/2;
			}
			loss_n1 = calc_loss( sum_best_beta_n2, features_it, beta_coeficient_update_Z * rule.Z, penalize );
			if ( step_verbosity > 3 ) {
				Rcout << "Adjusted n1 to be midpoint of n0 and n2 due to only having 1 iteration" << endl;
				Rcout << "\trevised losses: loss_n0: " << loss_n0 << "\tloss_n1: " << loss_n1 << "\tloss_n2: " << loss_n2 << endl;
				Rcout << "\txbeta_00s: " << sum_best_beta_n0[rule.doc_support[0]] << "\t" << sum_best_beta_n1[rule.doc_support[0]] << "\t" << sum_best_beta_n2[rule.doc_support[0]] << endl;

			}

		}


	}



	// Line search method. Binary search for optimal step size. Calls find_best_range(...).
	// sum_best_beta will get set to the scalar product beta_best^t*xi for each doc xi.
	// Instead of working with the new weight vector beta_n+1 obtained as beta_n - epsilon * gradient(beta_n)
	// we work directly with the scalar product.
	// We output the sum_best_beta_opt which contains the scalar product of the optimal beta found,
	// by searching for the optimal
	// epsilon, e.g. beta_n+1 = beta_n - epsilon_opt * gradient(beta_n)
	// epsilon is the starting value
	// rule contains info about the gradient at the current iteration
	// xbeta_0 is the Xbeta vector (one value for each document) of the starting point.
	void binary_line_search(space_t& rule, vector<LONG_DOUBLE>& xbeta_0, vector<LONG_DOUBLE>* sum_best_beta_opt) {

		// Keep track of scalar product at points beta_n-1, beta_n and beta_n+1.
		// They are denoted with beta_n0, beta_n1, beta_n2.
		vector<LONG_DOUBLE> sum_best_beta_n0(xbeta_0);
		vector<LONG_DOUBLE> sum_best_beta_n1(xbeta_0);
		vector<LONG_DOUBLE> sum_best_beta_n2(xbeta_0);

		find_bracketing_range( rule, sum_best_beta_n0, sum_best_beta_n1, sum_best_beta_n2, rule.penalize, pos_only );

		if ( step_verbosity > 3 ) {
			// print out the three bracketing vectors
			for ( unsigned int i = 0; i < corpus.size(); i++ ) {
				Rcout << sum_best_beta_n0[i] << "\t" << sum_best_beta_n1[i] << "\t" << sum_best_beta_n2[i] << endl;
			}
		}


		find_best_range(sum_best_beta_n0, sum_best_beta_n1, sum_best_beta_n2,
                        rule, xbeta_0[rule.doc_support[0]], sum_best_beta_opt, rule.penalize);

	} // end binary_line_search().


	// Find intercept.  Uses binary line search.  Possibly this function should
	// be merged into other one.
	// return the change in the intercept
	LONG_DOUBLE intercept_search(vector<LONG_DOUBLE>* sum_best_beta_opt) {
		vector<LONG_DOUBLE> sum_best_beta_n0(*sum_best_beta_opt);
		vector<LONG_DOUBLE> sum_best_beta_n1(*sum_best_beta_opt);
		vector<LONG_DOUBLE> sum_best_beta_n2(*sum_best_beta_opt);

		//Rcout << "intercept search" << endl;

		// adjust intercept up or down?  Calc gradient
		int_rule.gradient = calc_int_gradient(*sum_best_beta_opt);

		// find range, no regularization since this is intercept.
		find_bracketing_range( int_rule, sum_best_beta_n0, sum_best_beta_n1, sum_best_beta_n2, false, false );

		vector<LONG_DOUBLE> sum_best_beta_mid_n0_n1(corpus.size());
		vector<LONG_DOUBLE> sum_best_beta_mid_n1_n2(corpus.size());

		LONG_DOUBLE xbeta_00 = (*sum_best_beta_opt)[int_rule.doc_support[0]];
		if ( step_verbosity > 3 ) {
			Rcout << "Original xbeta_00: " << xbeta_00 << endl;
		}
		find_best_range(sum_best_beta_n0, sum_best_beta_n1, sum_best_beta_n2,
                        int_rule, xbeta_00, sum_best_beta_opt, false );

		LONG_DOUBLE int_step = (*sum_best_beta_opt)[int_rule.doc_support[0]] - xbeta_00;
		if ( step_verbosity > 3 ) {
			Rcout << "Intercept step found to be " << int_step << endl;
		}
		return int_step;
	} // end intercept_search().



public:


	void make_unigram_list( ) {

		// A map from unigrams to search_space.
		string unigram;
		bool at_space = false;

		// Prepare the locations for unigrams.
		if (verbosity >= 1) {
			Rcout << "preparing inverted index for unigrams\n";
			Rcout.flush();
		}
		unsigned int l   = corpus.size();

		for (unsigned int docid = 0; docid < l; ++docid) {
			//Rcout << "\nscanning docid: " << docid << " class y: " << y[docid];
			if ( verbosity >= 9 ) {
				Rcout << "doc #" << docid << ": " << ( (corpus[docid].size() < 30) ? corpus[docid] : corpus[docid].substr(0,30)) << endl;
			}

			for (unsigned int pos = 0; pos < corpus[docid].size(); ++pos) {
				// Skip white spaces. They are not considered to be unigrams.
				if (isspace(corpus[docid][pos])) {
					at_space = true;
					continue;
				}
				// If word level tokens.
				if (token_type==WORD_TOKEN) {
					if (at_space || pos == 0) {
						at_space = false;
						if (!unigram.empty()) {
							space_t & tmp = unigrams[unigram];
							tmp.set_ne( unigram );
							if ( verbosity > 10 ) {
								Rcout << "Adding unigram.  At doc #" << docid << ", char " << pos << " doclen: " << corpus[docid].size() << " un sz: " << unigram.size() << endl;
							}
							tmp.add (docid, pos - unigram.size() - 1);
							//						if (unigram.size() >= 1)
							//    						Rcout << "\nunigram:" << unigram << "*";
							unigram.clear();
						}
						unigram.push_back(corpus[docid][pos]);
					} else {
						unigram.push_back(corpus[docid][pos]);
					}
				} else {  // Char (i.e. byte) level token.
					unigram = corpus[docid][pos];
					space_t & tmp = unigrams[unigram];
					tmp.add (docid, pos);
					tmp.set_ne( unigram );
					unigram.clear();
				}
			}

			// final unigram of document
			if (token_type==WORD_TOKEN) {
				if (!unigram.empty()) {
					space_t & tmp = unigrams[unigram];
					tmp.set_ne( unigram );
					if ( verbosity > 10 ) {
						Rcout << "Adding final unigram.  doc #" << docid << " doclen: " << corpus[docid].size() << " un sz: " << unigram.size() << endl;
					}
					tmp.add (docid, corpus[docid].size() - unigram.size());
					unigram.clear();
				}
			} //end for transaction.
		} //end for docid.

		gettimeofday(&t, NULL);
		if (verbosity >= 1) {
			Rcout << " ( " << (t.tv_sec - t_origin.tv_sec) << " seconds; " << (t.tv_sec - t_origin.tv_sec) / 60.0 << " minutes )\n";
			Rcout.flush();
		}

		if (verbosity >= 1) {
			Rcout << "\n# distinct unigrams before culling: " << unigrams.size() << endl;
			Rcout.flush();
		}
    }

    // Drop any unigram with insufficient support.
    // Also build cache of possible phrase extensions.
    void cull_unigram_list() {


		// Keep only unigrams above minsup threshold.
		for (std::map <string, space_t>::iterator it = unigrams.begin (); it != unigrams.end(); ) {

            if ( !it->second.isConverted() ) {
                it->second.calc_support_weights( Lp, binary_features, no_regularization );
            }
			if ( it->second.support() < minsup ) {
				//Rcout << "erasing unigram " << it->second.ne << endl;
				if ( verbosity >= 5 ) {
					Rcout << "killing " << it->first << "/" << it->second.support() << endl;
				}
				std::map <string, space_t>::iterator killptr = it;
				++it;
				unigrams.erase (killptr);
			} else {
				single_node_minsup_cache.insert (it->second.ne);
				if (verbosity >= 5) {
					Rcout << "distinct unigram: " << it->first << "/" << it->second.support() << endl;
				}
				++it;
			}
		}

		gettimeofday(&t, NULL);
		if (verbosity >= 1) {
			//			Rcout << "done erasing unigrams without minimal support" << endl;
			Rcout << "\n# distinct unigrams: " << single_node_minsup_cache.size();
			Rcout << " ( " << (t.tv_sec - t_origin.tv_sec) << " seconds; " << (t.tv_sec - t_origin.tv_sec) / 60.0 << " minutes )\n";
			Rcout.flush();
		}
	}



	bool setup_output_file( std::ofstream &os, const char*out_file) {
		// setup output file
		os.open (out_file);
		if (! os) {
			Rcerr << "FATAL: Cannot open output file: " << out_file << std::endl;
			return false;
		}
		if ( verbosity > 0 ) {
			Rcout << "opened outfile '" << out_file << "'\n";
		}

		os.setf(std::ios::fixed,std::ios::floatfield);
		os.precision(12);

		return true;
	}


	/**
	 * Add ngram step to model
	 *
	 * 1) Write step out to file.
	 * 2) Update the stored model in memory.
	 * 3) Add ngram to rule cache if necessary.
	 */
	void add_to_model( space_t *rule, LONG_DOUBLE step, std::ofstream &os ) {
		string &ngram = rule->ngram;
		used_rule_cache[rule->ngram] = rule;
		if ( os.is_open() ) {
			os << step << '\t' << ngram << std::endl;
		} else {
			history_steps.push_back( step );
			history_ngrams.push_back( ngram );
		}
		map<string, LONG_DOUBLE>::iterator features_it = final_model_cache.find(ngram);
		if (features_it == final_model_cache.end()) {
			final_model_cache[ ngram ] = step;
			if ( verbosity >=4 ) {
				Rcout << "Adding feature " << ngram << " with first value of " << step << endl;
			}
		} else {
			final_model_cache[ ngram ] += step;
			if ( verbosity >=4 ) {
				Rcout << "Adding " << step << " to feature " << ngram << endl;
			}
		}
	}


    /*
     * Find best ngram and stash in 'rule'
     */
	void find_best_ngram() {
		tau = 0;
		rule = NULL;
		//int_rule.gradient = 0;
		//rule = &int_rule;

		std::vector <space_t*> old_space;
		std::vector <space_t*> new_space;

		//gettimeofday(&t_start, NULL);
		pruned = total = rewritten = 0;
		old_space.clear ();
		new_space.clear ();

		// start with intercept
		//int_rule.gradient = calc_int_gradient( sum_best_xbeta );
		//update_best_rule_to_intercept( int_rule.gradient );
		//if ( verbosity > 3 ) {
		//			Rcout << "Starting Point" << endl;
		//			print_rule( int_rule );
		//}


		// check all features in cache
		check_model_features();

		if ( verbosity > 2 ) {
			Rcout << "Searching " << unigrams.size() << " unigrams" << endl;
		}

		// If BFS traversal.
		if (traversal_strategy == BFS) {

			// Iterate through unigrams to make first list.
			for (std::map <string, space_t>::iterator it = unigrams.begin (); it != unigrams.end(); ++it) {

				// TODO: should not need---banned words are already dumped, right?   (but can't due to batch
				// mode for multiple labelings and different bans.  Think this through later.)
				if ( banned_words.count(it->second.ne) > 0 ) {
					continue;
				}

				space_t *uni = &(it->second);

				map<string, LONG_DOUBLE>::iterator features_it = features_cache.find(uni->ngram);
				if (features_it == features_cache.end()) {
					calc_gradient_and_update( uni );
				}

				old_space.push_back ( uni );
			}

			// Search for best n-gram. Extend in a bfs fashion,
			// level per level, e.g., first extend unigrams to bigrams, then bigrams to trigrams, etc.
			for (unsigned int size = 2; size <= maxpat; ++size) {
				if ( verbosity > 2 ) {
					Rcout << "Searching " << old_space.size() << " clusters of " << size << "-grams" << endl;
				}

				for (unsigned int i = 0; i < old_space.size(); ++i) {
					span_bfs (old_space[i], new_space, size);
				}
				if (new_space.empty()) {
					break;
				}
				old_space = new_space;
				new_space.clear ();
			} // End search for best n-gram.

		} else { // DFS traversal.
			// Iterate through unigrams and go deep
			for (std::map <string, space_t>::iterator it = unigrams.begin (); it != unigrams.end(); ++it) {
				if ( verbosity > 3 ) {
					Rcout << "Descending from unigram " << it->second.ne << endl;
				}

				// TODO: should not need---they are already dumped, right?
				if ( banned_words.count(it->second.ne) > 0 ) {
					//Rcout << "\tSkipping" << it->second.ngram << endl;
					continue;
				}

				space_t *uni = &(it->second);
				// first check to see if we are the new best bet.
				map<string, LONG_DOUBLE>::iterator features_it = features_cache.find(uni->ngram);
				if (features_it == features_cache.end()) {
					calc_gradient_and_update( &(it->second) );
				}

				if (!can_prune ( uni )) {
					span_dfs (uni, 2);
				}
			}
		}

		// Keep best ngram rule.
		//rule_cache.insert (rule);
		if (verbosity >= 2) {
			if ( rule == NULL ) {
				Rcout << "\nFound no ngram with nonzero gradient." << endl;
			} else {
				Rcout << "\nfound best ngram! (" << rule->ngram << ") ";
				rule->print_rule( );
			}
			gettimeofday(&t, NULL);
			Rcout << " ( " << t.tv_sec - t_origin.tv_sec << " seconds; " << (t.tv_sec - t_origin.tv_sec) / 60.0 << " minutes )\n";
		}
	}




	/**
	 * Once a step-length is found, add feature to the features cache if it is not there
	 * and update the sum squared betas and abs betas by the extra step length added to the
	 * current beta.
	 */
	void update_cumulative_weights( LONG_DOUBLE step_length_opt ) {
		// Insert or update new feature.
		if ( rule->penalize ) {
			map<string, LONG_DOUBLE>::iterator features_it = features_cache.find(rule->ngram);
			if (features_it == features_cache.end()) {
				// If feature not there, insert it.
				//Rcout << "adding ngram '" << rule->ngram << "' to feature cache" << endl;
				features_cache[rule->ngram] = step_length_opt;
			} else {
				// Adjust coeficient and the sums of coeficients: subtract out old value
				//	Rcout << "updating ngram '" << rule->ngram << "' to feature cache" << endl;
				sum_squared_betas = sum_squared_betas - pow(features_it->second, 2);
				sum_abs_betas = sum_abs_betas - abs(features_it->second);
				features_it->second += step_length_opt;
			}
			sum_squared_betas += pow(features_cache[rule->ngram], 2);
			sum_abs_betas += abs(features_cache[rule->ngram]);
		}
		//} else {
		//	Rcout << "not updating non-penalized feature '" << rule->ngram << "'" << endl;
		//}
	}


	/**
	 * Do a single gradient-descent step which consists of:
	 *    a) Finding the best ngram to update (i.e., find ngram with largest magnitude gradient)
	 *    b) Doing a line search and changing the coefficient for that ngram to minimize loss.
	 *
	 *    The results gets stored in best_beta_opt
	 */
	void descend_one_step( unsigned int itr, std::ofstream &os, std::vector < LONG_DOUBLE > &sum_best_beta_opt ) {

		if ( verbosity > 1 ) {
			Rcout << "\n\n** Descending one step.  Iteration #" << itr << "\n";
		}

		find_best_ngram();

		if ( verbosity > 3 ) {
			Rcout << "\nNow Optimizing beta for ngram.\n";
		}
		if ( rule == NULL || tau == 0 || rule->gradient == 0 ) {
			if ( verbosity > 0 ) {
				Rcout << "Best ngram has 0 gradient.  At minimum.  Returning without loss calculation.  Copying over sum_xbeta to opt." << endl;
			}
			sum_best_beta_opt.clear();
			sum_best_beta_opt.assign(sum_best_xbeta.begin(), sum_best_xbeta.end());
			return;
		}
		if ( verbosity > 3 ) {
			Rcout << "Starting point for sum_best_xbeta:\n";
			print_estimates( sum_best_xbeta, 1 );
		}

		// Use line search to detect the best step_length/learning_rate.
		// The method does binary search for the step length, by using a start parameter epsilon.
		// It works by first bracketing a promising value, followed by focusing on the smallest interval containing this value.
		sum_best_beta_opt.clear();
		binary_line_search( *rule, sum_best_xbeta, &sum_best_beta_opt);

		// The optimal step length as obtained from the line search.
		LONG_DOUBLE step_length_opt = extract_delta_beta( *rule, sum_best_beta_opt, sum_best_xbeta[rule->doc_support[0]] );
		if ( verbosity > 3 ) {
			Rcout << "\nOptimal step length for feature '" << rule->ngram << "' found: " << step_length_opt << endl;
		}

		// Update the weight of the best n-gram and calc elements of regularization
		update_cumulative_weights( step_length_opt );

		// calc current loss
		map<string, LONG_DOUBLE>::iterator features_it = features_cache.end();

		//Rcout << "time now: " << gettimeofday(t, NULL);

		// The optimal step length as obtained from the line search.
		add_to_model( rule, step_length_opt, os );

		if (verbosity >= 1) {
			LONG_DOUBLE loss = calc_loss( sum_best_beta_opt, features_it, 0, true );
			Rcout <<  "\nItr " << itr << " results: size model: " << features_cache.size () << "    rewrite/prune/total: "
            << rewritten << "/" << pruned << "/" << total << " "
            << "   # nodes: " << total_nodes << "\n\tgradient: " << rule->gradient
            << "\n\tstep len: "<< step_length_opt << "\n\tngram: " << rule->ngram;

			LONG_DOUBLE pen = calc_penalty(sum_abs_betas, sum_squared_betas, 0, 0 );
			Rcout << "\n\tloss + penalty term = " << (loss - pen) << " + " << pen << " = " << loss << endl;
			Rcout.flush();
		}

		if (verbosity >= 4 ) {
			print_out_model( final_model_cache, Rcout );
			print_estimates( sum_best_beta_opt, 1 );
		}
	} // end descend_one_step


	/**
	 * A slashed down version of descend_one_step to just adjust intercept
	 */
	void adjust_intercept( unsigned int itr, std::ofstream &os ) {

		if ( verbosity > 1 ) {
			Rcout << "\n\n** Adjusting intercept.  Iteration #" << itr << "\n";
		}

        //calc intercept
		LONG_DOUBLE int_step = intercept_search( &sum_best_xbeta );
		add_to_model( &int_rule, int_step, os );

		map<string, LONG_DOUBLE>::iterator features_it = features_cache.end();
		if (verbosity > 1) {
			LONG_DOUBLE loss = calc_loss( sum_best_xbeta, features_it, 0, false );
			Rcout <<  "\nIntercept Adjust #" << itr << ": " << "\n\tgradient: " << int_rule.gradient << "\n\tstep len: "<< int_step;

			LONG_DOUBLE pen = calc_penalty(sum_abs_betas, sum_squared_betas, 0, 0 );
			Rcout << "\n\tloss + penalty term = " << (loss - pen) << " + " << pen << " = " << loss << endl;
			Rcout.flush();

			if (verbosity >= 4 ) {
				Rcout << "printing final cache" << endl;
				print_out_model( final_model_cache, Rcout );
				Rcout << "printing best beta opt" << endl;
				print_estimates( sum_best_xbeta, 1 );
				Rcout << "done printing best beta opt" << endl;

			}
		}

	} // end adjust_intercept



	LONG_DOUBLE calc_convergence_rate( vector<LONG_DOUBLE>& xbeta_t0, vector<LONG_DOUBLE>& xbeta_t1 ) {
		// calc some stats to determine convergence
		LONG_DOUBLE sum_abs_scalar_prod_diff = 0;
		LONG_DOUBLE sum_abs_scalar_prod = 0;
		for ( unsigned int k  = 0; k < (unsigned int)corpus.size(); ++k ) {
			// Compute the sum of per document difference between the scalar product at 2 consecutive iterations.
			sum_abs_scalar_prod_diff += abs(xbeta_t1[k] - xbeta_t0[k]);
			// Compute the sum of per document scalar product at current iteration.
			sum_abs_scalar_prod += abs(xbeta_t1[k]);
		}

		if ( verbosity > 0 ) {
			Rcout << "Convergence rate: " <<  sum_abs_scalar_prod_diff << " / (1+" << sum_abs_scalar_prod << ") = " << (sum_abs_scalar_prod_diff / (1 + sum_abs_scalar_prod)) << endl;
		}

		// Set the convergence rate as in paper by Madigan et al on BBR.
		return sum_abs_scalar_prod_diff / (1 + sum_abs_scalar_prod);

	}


    void initialize_run() {
        // All beta coeficients are zero when starting.
        sum_squared_betas = 0;
        sum_abs_betas = 0;

        unsigned int l = corpus.size();

        // The starting point is beta = (0, 0, 0, 0, .....).
        sum_best_xbeta.resize(l);
        std::fill(sum_best_xbeta.begin(), sum_best_xbeta.end(), 0.0);
        exp_fraction.resize (l);
        std::fill (exp_fraction.begin(), exp_fraction.end(), 1.0 /2.0);

        // intercept is a rule with a weight of 1 on all the documents.
        int_rule.ngram = INTERCEPT_STRING;
        int_rule.doc_support.resize(l);
        int_rule.weight.resize(l);
        int_rule.penalize = false;
        for ( unsigned int i = 0; i < corpus.size(); i++ ) {
            int_rule.doc_support[i] = i;
            int_rule.weight[i] = 1;
        }
        int_rule.total_support = l;
        int_rule.total_docs = l;
        int_rule.Z = 1;

    }

	bool run (const char *out_file,
              unsigned int maxitr ) {


        if ( verbosity > 0 ) {
        	print_settings( maxitr );
        }

        std::ofstream os;
		std::ofstream final_os;
        bool write_to_file = out_file != NULL;

    	tau         = 0.0;
		total_nodes = 0;

        initialize_run();


		Rcout.setf(std::ios::fixed,std::ios::floatfield);
		Rcout.precision(8);

		if ( verbosity > 0 ) {
			Rcout << "Beginning run\n";
		}

		if ( write_to_file ) {

			string dtxt = string(".txt");
			string model_file = string(out_file);
			size_t found = model_file.find(dtxt);
			if (found!=string::npos) {
				model_file.erase (found,dtxt.length());
			}
			model_file += string("_model.txt");
			if ( ! setup_output_file(os, model_file.c_str()) ) {
				return false;
			}


			if ( ! setup_output_file(final_os, out_file) ) {
				return false;
			}
		}

		features_cache.clear();
        history_steps.clear();
        history_ngrams.clear();
        used_rule_cache.clear();
        final_model_cache.clear();
        features_cache.clear();

        gettimeofday(&t_origin, NULL);

		LONG_DOUBLE convergence_rate;


		// Compute loss with start beta vector.
		map<string, LONG_DOUBLE>::iterator features_it = features_cache.end();
		LONG_DOUBLE loss = calc_loss( sum_best_xbeta, features_it, 0, false );
		if (verbosity >= 1) {
			Rcout << "start loss: " << loss << endl;
		}


/*		// Iterate through unigrams to set them up.
		for (std::map <string, space_t>::iterator it = unigrams.begin (); it != unigrams.end(); ++it) {
            space_t *uni = &(it->second);
            uni->calc_support_weights( Lp, binary_features );
		}
  */
		std::vector < LONG_DOUBLE > 				sum_best_beta_opt;


		//Rcout << "\nstart iterations... ";
		// Loop for number of given optimization iterations.
		itr = 0;
		unsigned int flat_hits = 0;
		for (itr = 0; itr < maxitr; ++itr) {

			// adjust intercept (this modifies sum_best_xbeta)
			adjust_intercept( itr, os );

			// find rule and descend
			descend_one_step( itr, os, sum_best_beta_opt );

			convergence_rate = calc_convergence_rate(sum_best_xbeta, sum_best_beta_opt );

			//sum_best_beta_opt is the optimum just found using line search.  Save it.
			//(Even if we have converged.)
			sum_best_xbeta.assign(sum_best_beta_opt.begin(), sum_best_beta_opt.end());

			if (convergence_rate < convergence_threshold) {
				flat_hits++;
				if ( flat_hits > 4 ) {
					if (verbosity >= 1) {
						Rcout << "\nFinish iterations due to convergence test!";
						Rcout << "\n# iterations: " << itr + 1;
					}
					break;
				} else {
					if (verbosity >= 1) {
						Rcout << "\nFlatness found on iteration # " << itr + 1 << endl;
					}

				}
			}


			//gettimeofday(&t, NULL);
			//Rcout << "\ntime 1 boosting iter: " << t.tv_sec - t_start.tv_sec << " seconds; " << (t.tv_sec - t_start.tv_sec) / 60.0 << " minutes ";
		} //end optimization iterations.

		if ( itr == maxitr && verbosity > 0 ) {
			Rcout << "\nFinished iterations due to hitting max iteration cut-off of " << maxitr << ".  May not have converged.\n";
		}

		gettimeofday(&t, NULL);
		if (verbosity >= 1) {
			features_it = features_cache.end();
			loss = calc_loss( sum_best_xbeta, features_it, 0, false );
			Rcout << "\nend loss: " << loss;
			Rcout << "\ntotal time: " << t.tv_sec - t_origin.tv_sec << " seconds; " << (t.tv_sec - t_origin.tv_sec) / 60.0 << " minutes\n ";
		}

        // for debugging, we can print small corpus entirely for diagnostic purposes.
		if ( corpus.size() < 20 && (verbosity > 0 )) {
			Rcout << "Individual document estimates:" << endl;
            print_estimates( sum_best_beta_opt, 0 );
		}

		if ( verbosity > 6 ) {
			Rcout << "Writing out final model:" << endl;
			Rcout << "------------------------" << endl;
			print_full_model( final_model_cache, itr, maxitr, Rcout );
			Rcout << "------------------------" << endl;
		}

		if ( write_to_file ) {
			print_full_model( final_model_cache, itr, maxitr, final_os );
		}

		if ( corpus.size() < 20 && verbosity > 4 ) {
			Rcout << "Writing out final model:" << endl;
			print_estimates( sum_best_xbeta, 1 );
		}

		return true;
	} //end run().



	NumericVector find_C ( unsigned int num_scrambles ) {
        Rcpp::RNGScope scope;         		// needed when RNGs are drawn

        if ( verbosity > 0 ) {
        	print_settings( 0 );
        }

        NumericVector res( num_scrambles + 1 );

		unsigned int l = corpus.size();

        Rcpp::NumericVector swps;

		sum_best_xbeta.resize(l);

		if ( verbosity > 1 ) {
            Rcout << "Looking for C to get empty model: " << num_scrambles << " iterations with first not scrambled"  << endl;
        }
        int n = y.size();

        for ( unsigned int sc_itr = 0; sc_itr < num_scrambles + 1; sc_itr++ ) {
            if ( verbosity > 1 ) {
                Rcout << "\nBeginning find_C iteration\n";
                Rcout.flush();
            }

            tau         = 0.0;
            total_nodes = 0;

            // All beta coeficients are zero when starting.
            sum_squared_betas = 0;
            sum_abs_betas = 0;

            features_cache.clear();


            if ( sc_itr > 0 ) {
                swps = Rcpp::runif(n);
                for ( int i = n - 1; i > 0; i-- ) {
                    swap( y[i], y[ int( swps[i]*n ) ] );
                }
            }


            // The starting point is beta = (0, 0, 0, 0, .....).
            std::fill(sum_best_xbeta.begin(), sum_best_xbeta.end(), 0.0);
            exp_fraction.resize (l);
            std::fill (exp_fraction.begin(), exp_fraction.end(), 1.0 /2.0);

            // intercept is a rule with a weight of 1 on all the documents.
            int_rule.ngram = "*intercept*";
            int_rule.doc_support.resize(l);
            int_rule.weight.resize(l);
            int_rule.penalize = false;
            for ( unsigned int i = 0; i < corpus.size(); i++ ) {
                int_rule.doc_support[i] = i;
                int_rule.weight[i] = 1;
            }
            int_rule.total_support = l;
            int_rule.total_docs = l;
            int_rule.Z = 1;

            Rcout.setf(std::ios::fixed,std::ios::floatfield);
            Rcout.precision(8);

             // adjust intercept (this modifies sum_best_xbeta)
            LONG_DOUBLE int_step = intercept_search( &sum_best_xbeta );
            if ( verbosity > 4 ) {
                Rcout << "Took step of size " << int_step << endl;
            }
            if ( verbosity > 1 ) {
                Rcout << "Finding best ngram\n";
            }

            // find best rule
            find_best_ngram();

            double needed_C = 0;
            if ( rule != NULL ) {
                needed_C = rule->gradient;

                if (verbosity >= 1) {
                    Rcout <<  "\nItr " << sc_itr << " results: size model: " << features_cache.size () << "    rewrite/prune/total: "
                    << rewritten << "/" << pruned << "/" << total << " "
                    << "   # nodes: " << total_nodes << "\n\tgradient: " << rule->gradient
                    << "\n\tFound C: "<< needed_C << "\n\t due to ngram: " << rule->ngram << endl;
                    Rcout.flush();
                }

            } else {
                if (verbosity >= 1) {
                    Rcout << "No ngram found.  C is therefore " << needed_C << endl;
                    Rcout.flush();
                }
            }

            res[sc_itr] = std::abs( needed_C );
        }

		gettimeofday(&t, NULL);
		if (verbosity >= 1) {
			Rcout << "\ntotal time: " << t.tv_sec - t_origin.tv_sec << " seconds; " << (t.tv_sec - t_origin.tv_sec) / 60.0 << " minutes\n ";
		}

		return res;
	} //end find_C().

};



SEXP textreg(Rcpp::XPtr<SeqLearner> seql_learner, Rcpp::List rparam) {

    // Step 1: read parameters for how to run ngram from R
    int maxiter = 100;
    bool find_C = false;
    int find_C_iter = 0;

    if ( seql_learner->verbosity > 1 ) {
        Rcout << "beginning c++ ngram function call\n";
        Rcout.flush();
    }

    // get parameters
    try {
        seql_learner->C = Rcpp::as<double>(rparam["C"]);
        seql_learner->pos_only = Rcpp::as<int>(rparam["positive.only"]);
        seql_learner->maxpat = Rcpp::as<int>(rparam["max.pattern"]);
        seql_learner->minpat = Rcpp::as<int>(rparam["min.pattern"]);

        unsigned int newminsup = Rcpp::as<int>(rparam["min.support"]);
        if ( seql_learner->has_unigram_cache() && newminsup != seql_learner->minsup ) {
            if ( seql_learner->verbosity > 0 ) {
                Rcerr << "Warning: Resetting tree due to min support parameter change from " << seql_learner->minsup << " to " << newminsup << endl;
            }
            seql_learner->reset_tree();
        }
        seql_learner->minsup = newminsup;

        seql_learner->maxgap = Rcpp::as<int>(rparam["gap"]);
        //  seql_learner->token_type = Rcpp::as<int>(rparam["token.type"]);
        seql_learner->convergence_threshold = Rcpp::as<double>(rparam["convergence.threshold"]);
        seql_learner->verbosity = Rcpp::as<double>(rparam["verbosity"]);
        seql_learner->step_verbosity = Rcpp::as<double>(rparam["step.verbosity"]);

        double newLp = Rcpp::as<double>(rparam["Lq"]);
        if ( seql_learner->has_unigram_cache() && seql_learner->Lp != newLp ) {
            if ( seql_learner->verbosity > 0 ) {
                Rcerr << "Warning: Resetting tree due to Lq parameter change" << endl;
            }
            seql_learner->reset_tree();
        }
        //Rcout << "Setting Lp to " << newLp << endl;
        seql_learner->set_Lp( newLp );  // WARNING: p is q in r package.

        seql_learner->positive_weight = Rcpp::as<int>(rparam["positive.weight"] );

        seql_learner->binary_features = Rcpp::as<int>(rparam["binary.features"] );
        seql_learner->no_regularization = Rcpp::as<int>(rparam["no.regularization"] );
        maxiter       = Rcpp::as<int>(rparam["maxIter"]);

        find_C       = Rcpp::as<int>(rparam["findC"]);
        find_C_iter  = Rcpp::as<int>(rparam["findCIter"]);

    } catch (std::exception &ex) {		// or use END_RCPP macro
        Rprintf("Caught error\n" );
        Rcerr << "error diagnostic '" << ex.what() << "'" << endl;
        //Rcerr << ex << endl;

        forward_exception_to_r( ex );
    } catch(...) {
    	::Rf_error( "c++ exception (unknown reason)" );
    }

    if ( seql_learner->verbosity > 1 ) {
        Rcout << "parameters loaded\n";
        Rcout.flush();
    }
    // Rprintf("Checking Labeling Vector\n" );

    // Step 3: Find all unigrams if not present.
    if ( !seql_learner->has_unigram_cache() ) {
        seql_learner->make_unigram_list();
        seql_learner->cull_unigram_list();
    } else {
        if ( seql_learner->verbosity >= 1 ) {
            Rcout << "No need to build unigram list as it is inherited from previous call" << endl;
        }
    }

    if ( !find_C ) {
        try {
            seql_learner->run( NULL, maxiter );

            if ( seql_learner->verbosity > 1 ) {
                Rcout << "assembling results to return" << endl;
                Rcout.flush();
            }
            DataFrame df = seql_learner->make_full_model_dataframe(seql_learner->final_model_cache);
            List ruleset = seql_learner->make_rule_set( seql_learner->final_model_cache );
            List search_path = seql_learner->make_search_path();
            List notes = Rcpp::List(rparam);
            notes[ "iter" ] = seql_learner->itr;
            notes[ "n" ] = (unsigned int)seql_learner->corpus.size();
            notes[ "Lp" ] = seql_learner->Lp;
            notes[ "binary.features" ] = seql_learner->binary_features;
            notes[ "no.regularization"] = seql_learner->no_regularization;
            notes[ "positive.weight" ] = seql_learner->positive_weight;

            if ( seql_learner->verbosity > 1 ) {
                Rcout << "going to return" << endl;
                Rcout.flush();
            }

            return Rcpp::List::create(_[ "model" ] = df,
                                  _[ "rules" ] = ruleset,
                                  _[ "notes" ] = notes,
                                  _[ "labeling" ] =  seql_learner->y,
 //                                 _[ "banlist" ] = bannedV,  // This was commented out at some point?  Need to bring back somehow
 //TODO: Add a get ban list to seql_learner
                                  _[ "path"] = search_path );


        } catch( std::exception &ex ) {		// or use END_RCPP macro
            forward_exception_to_r( ex );
        } catch(...) {
            ::Rf_error( "c++ exception (unknown reason)" );
        }
    } else {
        if ( seql_learner->verbosity > 1 ) {
            Rcout << "calling find_C\n";
            Rcout.flush();
        }
        try {
            NumericVector res = seql_learner->find_C( find_C_iter );

            return res;

        } catch( std::exception &ex ) {		// or use END_RCPP macro
            forward_exception_to_r( ex );
        } catch(...) {
            ::Rf_error( "c++ exception in find C block (unknown reason)" );
        }

    }
    return R_NilValue;
}


Rcpp::XPtr<SeqLearner> build_corpus(Rcpp::StringVector corpusV, Rcpp::NumericVector labelV,
									Rcpp::StringVector bannedV, Rcpp::List rparam)
{
    // Rcpp::StringVector bannedV = Rcpp::StringVector();
	SeqLearner *seql_learner = new SeqLearner();
    std::string outfile = "none";


    if ( seql_learner->verbosity > 1 ) {
        Rcout << "beginning c++ ngram function call\n";
        Rcout.flush();
    }


    // get parameters
    try {
        seql_learner->verbosity = Rcpp::as<double>(rparam["verbosity"]);
//        seql_learner->minsup = Rcpp::as<int>(rparam["min.support"]);
//        seql_learner->set_Lp( Rcpp::as<double>(rparam["Lq"]) );  // WARNING: p is q in r package.
//        seql_learner->binary_features = Rcpp::as<int>(rparam["binary.features"] );


    } catch (std::exception &ex) {		// or use END_RCPP macro
        Rprintf("Caught error\n" );
        Rcerr << "error diagnostic '" << ex.what() << "'" << endl;
        //Rcerr << ex << endl;

        forward_exception_to_r( ex );
    } catch(...) {
    	::Rf_error( "c++ exception (unknown reason)" );
    }

    if ( seql_learner->verbosity > 1 ) {
        Rcout << "parameters loaded\n";
        Rcout.flush();
    }
    // Rprintf("Checking Labeling Vector\n" );


    try {
    	//    	Rprintf("Adding Banned Words\n" );
    	// add banned words
    	// if ( ! Rf_isNull(banned) ) {

			for ( int i = 0; i < bannedV.size(); i++ ) {
				std::string stt = std::string(bannedV[i]);
				seql_learner->add_banned_word( stt );
			}
		// }

    	// Step 2: Get the corpus vector
    	//    	Rprintf("Labeling vector is %d long\n", labelV.size() );
    	//    	Rprintf("Corpus is %d long\n", corpusV.size() );

    	// did we get a file name or the actual text from the corpus?
    	if ( corpusV.size() == 1 && labelV.size() != 1 ) {

    		bool res = seql_learner->read_in_data( std::string(corpusV[0]).c_str(), labelV );
    		if ( !res ) {
    			//Rcerr << "Failing to grab " << std::string(corpusV[0]) << endl;
    			::Rf_error( "Failed to read in file data (file not found?)" );
    		} else {
    			if ( seql_learner->verbosity > 0 ) {
    				Rcout << "Finished reading in from text file" << endl;
    			}
    		}
    	} else {
			for ( int i = 0; i < corpusV.size(); i++ ) {
				std::string stt = std::string(corpusV[i]);
				seql_learner->add_document( stt, labelV[i] );
			}
	   	}
		seql_learner->finish_initializing();


    } catch( std::exception &ex ) {		// or use END_RCPP macro
    	forward_exception_to_r( ex );
    } catch(...) {
    	::Rf_error( "c++ exception (unknown reason)" );
    }

	// return the external pointer to the R side
	return Rcpp::XPtr<SeqLearner>( seql_learner, true );
}


Rcpp::XPtr<SeqLearner> update_banned(Rcpp::XPtr<SeqLearner> seql_learner, Rcpp::StringVector bannedV )
{
    if ( seql_learner->verbosity >= 1 ) {
        Rcout << "Updating ban list\n";
        Rcout.flush();
    }

    try {
        seql_learner->banned_words.clear();

        // add banned words
        for ( int i = 0; i < bannedV.size(); i++ ) {
            std::string stt = std::string(bannedV[i]);
            seql_learner->add_banned_word( stt );
        }

    } catch( std::exception &ex ) {		// or use END_RCPP macro
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "c++ exception (unknown reason)" );
    }

    // return the external pointer to the R side
    return Rcpp::XPtr<SeqLearner>( seql_learner, true );
}




// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393


// update_banned
Rcpp::XPtr<SeqLearner> update_banned(Rcpp::XPtr<SeqLearner> seql_learner, Rcpp::StringVector bannedV);
RcppExport SEXP textreg_update_banned(SEXP seql_learnerSEXP, SEXP bannedVSEXP) {
    BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::XPtr<SeqLearner> >::type seql_learner(seql_learnerSEXP );
        Rcpp::traits::input_parameter< Rcpp::StringVector >::type bannedV(bannedVSEXP );
        Rcpp::XPtr<SeqLearner> __result = update_banned(seql_learner, bannedV);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}





// This part of the file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393


// textreg
SEXP textreg(Rcpp::XPtr<SeqLearner> seql_learner, Rcpp::List rparam);
RcppExport SEXP textreg_textreg(SEXP seql_learnerSEXP, SEXP rparamSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::XPtr<SeqLearner> >::type seql_learner(seql_learnerSEXP );
        Rcpp::traits::input_parameter< Rcpp::List >::type rparam(rparamSEXP );
        SEXP __result = textreg(seql_learner, rparam);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


// build_corpus
Rcpp::XPtr<SeqLearner> build_corpus(Rcpp::StringVector corpusV, Rcpp::NumericVector labelV, Rcpp::StringVector bannedV, Rcpp::List rparam);
RcppExport SEXP textreg_build_corpus(SEXP corpusVSEXP, SEXP labelVSEXP, SEXP bannedVSEXP, SEXP rparamSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::StringVector >::type corpusV(corpusVSEXP );
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type labelV(labelVSEXP );
        Rcpp::traits::input_parameter< Rcpp::StringVector >::type bannedV(bannedVSEXP );
        Rcpp::traits::input_parameter< Rcpp::List >::type rparam(rparamSEXP );
        Rcpp::XPtr<SeqLearner> __result = build_corpus(corpusV, labelV, bannedV, rparam);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
