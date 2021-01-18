#include <RnewMat.h>
#include <R_ext/Visibility.h>
#include <combinatorics.h>
#include <dataStructure.h>
#include <hyperg.h>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <iostream>
#include <cassert>
#include <climits>
#include <conversions.h>
#include <cmath>
#include <types.h>

// using pretty much:
using std::map;
using std::set;
using std::vector;
using std::accumulate;
using std::find;
using std::set_difference;
using std::count;
using std::max_element;
using std::max;
using std::min;

typedef std::vector<long double>::size_type indexType;


SEXP exhaustiveGaussian(// declaration
                        SEXP R_x, // (not centered!) design matrix (with colnames)
                        SEXP R_xcentered, // centered design matrix
                        SEXP R_y, // response vector
                        SEXP R_fpmaxs, // vector of maximum fp degrees
                        SEXP R_fppos, // corresponding vector of fp column indices
                        SEXP R_fpcards, // corresponding vector of power set cardinalities
                        SEXP R_nFps, // number of fp terms
                        SEXP R_fpnames, // names of fp terms
                        SEXP R_ucIndices, // vector giving _unc_ertainty custer indices (column -> which group)
                        SEXP R_ucTermList, // list for group -> which columns mapping
                        SEXP R_nUcGroups, // number of uncertainty groups
                        SEXP R_totalNumber, // cardinality of model space,
                        SEXP R_hyperparam, // hyperparameter a for hyper-g prior
                        SEXP R_priorType, // type of model prior?
                        SEXP R_nModels, // number of best models to be returned
                        SEXP R_verbose); // should progress been displayed?

SEXP samplingGaussian(// declaration
                      SEXP R_x, // (not centered!) design matrix (with colnames)
                      SEXP R_xcentered, // centered design matrix
                      SEXP R_y, // response vector
                      SEXP R_fpmaxs, // vector of maximum fp degrees
                      SEXP R_fppos, // corresponding vector of fp column indices
                      SEXP R_fpcards, // corresponding vector of power set cardinalities
                      SEXP R_nFps, // number of fp terms
                      SEXP R_fpnames, // names of fp terms
                      SEXP R_ucIndices, // vector giving _unc_ertainty custer indices (column -> which group)
                      SEXP R_ucTermList, // list for group -> which columns mapping
                      SEXP R_nUcGroups, // number of uncertainty groups
                      SEXP R_hyperparam, // hyperparameter a for hyper-g prior
                      SEXP R_priorType, // type of model prior?
                      SEXP R_nModels, // number of best models to be returned
                      SEXP R_verbose, // should progress been displayed?
                      SEXP R_chainlength, // how many times should a jump been made?
                      SEXP R_nCache); // size of models cache (an STL map)

SEXP logMargLik( //declaration
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha, // hyperparamater for hyper-g prior
                SEXP R_sst); // total sum of squares computed from y

SEXP postExpectedg( //declaration
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha); // hyperparamater for hyper-g prior

SEXP postExpectedShrinkage( //declaration
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha); // hyperparamater for hyper-g prior


// export to C interface ##########################################################################

extern "C"
{
  
static const R_CallMethodDef callMethods[] = {
  {"exhaustiveGaussian", (DL_FUNC) &exhaustiveGaussian, 16},
  {"samplingGaussian", (DL_FUNC) &samplingGaussian, 17},
  {"logMargLik", (DL_FUNC) &logMargLik, 5},
  {"postExpectedg", (DL_FUNC) &postExpectedg, 4},
  {"postExpectedShrinkage", (DL_FUNC) &postExpectedShrinkage, 4},
  {NULL, NULL, 0}
};

  // Register these functions with R:
  
void attribute_visible R_init_bfp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

} // extern "C"



// other functions ##########################################################################
void permPars(PosInt pos, // current position in parameter vector, starting from 0
              const fpInfo &currFp,
              const int &nUcGroups,
              modelPar mod,
              set<model> &space,
              const hyperPriorPars &hyp,
              const dataValues &data,
              const vector<IntSet>& ucTermList,
              const set<int> &fixedCols,
              book&);

set<int> getFreeUcs( // compute set of free uc group indices
                   const modelPar& mod,
                   const vector<PosInt>& ucSizes,
                   const PosInt& currDim,
                   const PosInt& maxDim);

set<PosInt> getFreeCovs( // compute set of free cov indices
                             const modelPar& mod,
                             const fpInfo& currFp,
                             const set<int>& freeUcs,
                             const PosInt& currDim,
                             const PosInt& maxDim);
set<PosInt> getPresentCovs( // determine set of present cov indices
        const modelPar& mod);

template <class T> T discreteUniform( // return random element of myset; should be enclosed in getRNGstate() etc.
const set<T>& myset);

template <class T> typename T::iterator dU( // return iterator of random element of myset; should be enclosed in getRNGstate() etc.
const T& myset);

int discreteUniform( // get random int x with lower <= x < upper; should be enclosed in getRNGstate() etc.
                    const int& lower,
                    const int& upper);

void computeModel(const modelPar &mod,
                  const hyperPriorPars &hyp,
                  const dataValues &data,
                  const fpInfo &currFp,
                  const vector<IntSet>& ucTermList,
                  const int &nUcGroups,
                  const set<int> &fixedCols,
                  set<model> &space,
                  book&);

ReturnMatrix getDesignMatrix( // construct design matrix for the model
                             const modelPar &mod,
                             const dataValues &data,
                             const fpInfo &currFp,
                             const vector<IntSet>& ucTermList,
                             const int &nUcGroups,
                             const set<int> &fixedCols);

double getR2( // compute coefficient of determination for the model
             const Matrix &design,
             const dataValues &data,
             const set<int> &fixedCols,
             const hyperPriorPars &hyp);

double getVarLogMargLik( // compute varying part of log marginal likelihood for specific model
                        const double &R2,
                        const int &n,
                        const int &dim,
                        const hyperPriorPars &hyp);

double getVarLogPrior( // compute logarithm of model prior
                      const modelPar &mod,
                      const fpInfo &currFp,
                      const PosInt nUcGroups,
                      const hyperPriorPars &hyp);

ReturnMatrix getFpMatrix( // build Fp basis matrix from transformed cols and power indices
                         const vector<ColumnVector> &tcols,
                         const Powers &powerinds,
                         const dataValues &data);

void pushInclusionProbs( // push back index into covGroupWisePosteriors-Array
                        const modelPar &mod,
                        const fpInfo &currFp,
                        const int &nUcGroups,
                        book &bookkeep);




// definitions ####################################################################################

SEXP
exhaustiveGaussian(// definition
                   SEXP R_x, // design matrix
                   SEXP R_xcentered, // centered design matrix
                   SEXP R_y, // response vector
                   SEXP R_fpmaxs, // vector of maximum fp degrees
                   SEXP R_fppos, // corresponding vector of fp column indices
                   SEXP R_fpcards, // corresponding vector of power set cardinalities
                   SEXP R_nFps, // number of fp terms
                   SEXP R_fpnames, // names of fp terms
                   SEXP R_ucIndices, // vector giving _unc_ertainty custer indices (column -> which group)
                   SEXP R_ucTermList, // list for (group -> which columns) one-to-many mapping
                   SEXP R_nUcGroups, // number of uncertainty groups
                   SEXP R_totalNumber, // cardinality of model space
                   SEXP R_hyperparam, // hyperparameter a for hyper-g prior
                   SEXP R_priorType, // type of model prior?
                   SEXP R_nModels, // number of best models to be returned
                   SEXP R_verbose) // should progress been displayed?
{

    PosInt nProtect = 0;

    // unpack ###
	// data
	const Matrix x = getMatrix(R_x);
	const Matrix xcentered = getMatrix(R_xcentered);
	const ColumnVector y = vec2col(R_y);

	// prior specifications
	const double hyperparam = Rf_asReal(R_hyperparam);
	const std::string priorType = getStringVector(R_priorType).at(0);

	hyperPriorPars hyp(hyperparam,
	                   priorType);

	const double totalNumber = REAL(R_totalNumber)[0]; // cardinality of model space

	// constant information
	const dataValues data(x, xcentered, y, totalNumber);

	// fp info
	fpInfo currentFpInfo(R_nFps,
	                     R_fpcards,
	                     R_fppos,
	                     R_fpmaxs,
	                     R_fpnames,
	                     x);

	// uc info
	const int* ucIndicesArray = INTEGER(R_ucIndices);
	const vector<int> ucIndices(ucIndicesArray, ucIndicesArray + Rf_length(R_ucIndices));
	const int nUcGroups = INTEGER(R_nUcGroups)[0];

	vector<IntSet> ucTermList(nUcGroups); // Array with length nUcGroups of int-vectors
	if(nUcGroups){ // catch case with no uc groups
		for(R_len_t i = 0; i != Rf_length(R_ucTermList); i++) {
			SEXP temp = VECTOR_ELT(R_ucTermList, i);
			copy(INTEGER(temp), INTEGER(temp) + Rf_length(temp),
			     inserter(ucTermList.at(i), ucTermList.at(i).begin()));
		}
	}

	set<int> nonfixedCols = set<int>(currentFpInfo.fppos, currentFpInfo.fppos + currentFpInfo.nFps); // fixed info
	set<int> allCols, fixedCols;
	for(vector<int>::size_type i = 0; i != ucIndices.size(); i++){
		allCols.insert(i + 1);
		if (ucIndices.at(i))
			nonfixedCols.insert(i + 1);
	}
	set_difference(allCols.begin(), allCols.end(), nonfixedCols.begin(), nonfixedCols.end(),
				   inserter(fixedCols, fixedCols.begin()));
	// now fixedCols contains indices of columns that are always present in the design matrix


	// no map needed for exhaustive search, a set is the right thing:
	set<model> orderedModels;

	// compute marginal likelihood for every model ###
	if (orderedModels.max_size() < totalNumber)
		Rf_error("\nmodel space is too large - cannot compute every model\n");

	// start model
	modelPar startModel(currentFpInfo.nFps, 0, 0);
	PowersVector startFps(currentFpInfo.nFps); // allocate correct length of vector
	startModel.fpPars = startFps;

	// bookkeeping
	book bookkeep;
	bookkeep.verbose = LOGICAL(R_verbose)[0];

	// for computation of inclusion probs
	vector<indexSafeSum> cgwp(currentFpInfo.nFps + nUcGroups);
	bookkeep.covGroupWisePosteriors = cgwp;

	// vector<indexSafeSum> lfp(currentFpInfo.nFps);
	bookkeep.linearFpPosteriors = vector<indexSafeSum>(currentFpInfo.nFps);

	// how many models to return?
	bookkeep.nModels = INTEGER(R_nModels)[0];

	// start computation
	permPars(0, currentFpInfo, nUcGroups, startModel, orderedModels, hyp, data, ucTermList, fixedCols, bookkeep);

	if (bookkeep.verbose){
		Rprintf("\nActual number of possible models:  %d ", bookkeep.modelCounter);
		Rprintf("\nNumber of non-identifiable models: %d", bookkeep.nanCounter);
		Rprintf("\nNumber of saved possible models:   %d\n", orderedModels.size());
	}

	// normalize posterior probabilities and correct log marg lik and log prior of the models to return
	const long double normConst = bookkeep.modelPropToPosteriors.sum();
	const long double logNormConst = log(normConst);

	const double logMargLikConst = 	- (data.nObs - 1) / 2.0 * log(data.sumOfSquaresTotal) - log(hyp.a - 2.0);

	// inclusion probs
	SEXP inc;
	Rf_protect(inc = Rf_allocVector(REALSXP, currentFpInfo.nFps + nUcGroups));
	nProtect++;
	for (int i = 0; i != Rf_length(inc); i++)
		REAL(inc)[i] = bookkeep.covGroupWisePosteriors.at(i).sum(bookkeep.modelPropToPosteriors) / normConst;

	SEXP linearInc;
	Rf_protect(linearInc = Rf_allocVector(REALSXP, currentFpInfo.nFps));
	nProtect++;
	for (int i = 0; i != Rf_length(linearInc); i++)
	    REAL(linearInc)[i] = bookkeep.linearFpPosteriors.at(i).sum(bookkeep.modelPropToPosteriors) / normConst;

	SEXP ret;
	Rf_protect(ret = Rf_allocVector(VECSXP, orderedModels.size()));
	nProtect++;

	bookkeep.chainlength = 1; // prevent calculation with uninitialised variable in covert2list
	
	unsigned int i = 0;
	for(set<model>::reverse_iterator j = orderedModels.rbegin(); j != orderedModels.rend(); j++)
	{
		model modCopy = *j;
		SET_VECTOR_ELT(ret, i++, modCopy.convert2list(currentFpInfo, logMargLikConst, logNormConst, bookkeep));
	}
	Rf_setAttrib(ret, Rf_install("numVisited"), Rf_ScalarReal(bookkeep.modelCounter));
	Rf_setAttrib(ret, Rf_install("inclusionProbs"), inc);
	Rf_setAttrib(ret, Rf_install("linearInclusionProbs"), linearInc);
	Rf_setAttrib(ret, Rf_install("logNormConst"), Rf_ScalarReal(logNormConst));

	// return ###
	Rf_unprotect(nProtect);
	return ret;
}

// ***************************************************************************************************//

// recursion via:
void permPars(PosInt pos, // current position in parameter vector, starting from 0 - copied.
              const fpInfo& currFp,
              const int &nUcGroups,
              modelPar mod,	// is copied every time! everything else is call by reference.
              set<model> &space,
              const hyperPriorPars &hyp,
              const dataValues &data,
              const vector<IntSet>& ucTermList,
              const set<int> &fixedCols,
              book &bookkeep)
{
	if (pos != currFp.nFps){ // some fps are still left
		const int card = currFp.fpcards[pos]; // cardinality of this power set
		permPars(pos + 1, currFp, nUcGroups, mod, space, hyp, data, ucTermList, fixedCols, bookkeep); // degree 0
		for (int deg = 1; deg <= currFp.fpmaxs[pos]; deg++){ // different degrees for fp at pos
			mod.fpSize++; // increment sums of fp degrees
			IntVector part(card); // partition of deg into card parts
			bool more1 = false;
			int h(0), t(0); // internal variables for comp_next
			do {
				comp_next(deg, card, part, &more1, h, t);	// next partition of deg into card parts
				mod.fpPars[pos] = freqvec2multiset(part); // convert into multiset
				// and go on
				permPars(pos + 1, currFp, nUcGroups, mod, space, hyp, data, ucTermList, fixedCols, bookkeep);
			} while (more1);
		}
	} else { // no fps left
		computeModel(mod, hyp, data, currFp, ucTermList, nUcGroups, fixedCols, space, bookkeep);
		for (int deg = 1; deg <= nUcGroups; deg++){ // different number of uc groups
			mod.ucSize++; // increment number of uc groups present
			IntVector subset(deg); // partition of deg into card parts
			bool more2 = false;
			int m(0), m2(0); // internal variables for ksub_next
			do {
				ksub_next(nUcGroups, deg, subset, &more2, m, m2);	// next subset (positive integers)
				mod.ucPars = set<int>(subset.begin(), subset.end()); // convert into set
				computeModel(mod, hyp, data, currFp, ucTermList, nUcGroups, fixedCols, space, bookkeep);
			} while (more2);
		}
	}
}



// ***************************************************************************************************//

void computeModel(// compute (varying part of) marginal likelihood and prior of mod and insert into map
					const modelPar &mod,
					const hyperPriorPars &hyp,
					const dataValues &data,
					const fpInfo &currFp,
					const vector<IntSet>& ucTermList,
					const int &nUcGroups,
					const set<int> &fixedCols,
					set<model> &space,
					book &bookkeep
				 )
{
	static set<model>::size_type compCounter = 0;

	// design matrix
	Matrix thisDesign = getDesignMatrix(mod, data, currFp, ucTermList, nUcGroups, fixedCols);

	// R2
	double thisR2 = getR2(thisDesign, data, fixedCols, hyp);

	if (R_IsNaN(thisR2) == FALSE){
		// log marginal likelihood
		double thisVarLogMargLik = getVarLogMargLik(thisR2, data.nObs, thisDesign.Ncols(), hyp);

		// log prior
		const double thisLogPrior = getVarLogPrior(mod, currFp, nUcGroups, hyp);

		// posterior expected g
		double thisPostExpectedg = posteriorExpectedg_hyperg(thisR2, data.nObs, thisDesign.Ncols(), hyp.a, thisVarLogMargLik);

		// posterior expected shrinkage
		double thisPostExpectedShrinkage = posteriorExpectedShrinkage_hyperg(thisR2, data.nObs, thisDesign.Ncols(), hyp.a, thisVarLogMargLik);

		// put all this into the modelInfo
		modelInfo info(thisVarLogMargLik, thisLogPrior, thisPostExpectedg, thisPostExpectedShrinkage, thisR2);

		// altogether we have the model:
		model thisModel = model(mod, info);

		// and insert it into the model space
		if (space.size() >= bookkeep.nModels){
			set<model>::iterator it = space.begin();
			if (*it < thisModel){ // compare this model to the least probable model in the set
				space.erase(it);
				space.insert(thisModel); // exchange if it is better than this worst model in the set
			}
		} else {
			space.insert(thisModel);
		}

		const long double thisPropToPosterior = expl(thisVarLogMargLik + thisLogPrior);
		bookkeep.modelPropToPosteriors.add(thisPropToPosterior);

		pushInclusionProbs(mod, currFp, nUcGroups, bookkeep);
		bookkeep.modelCounter++;

	} else {
		bookkeep.nanCounter++;
	}
	// increase static vars
	if((++compCounter % max(data.totalNumber / 100, static_cast<dataValues::NumberType>(1)) == 0) && bookkeep.verbose)
		Rprintf("-"); // display computation progress at each percent
}

// ***************************************************************************************************//

SEXP
samplingGaussian(// definition
                 SEXP R_x, // design matrix (with colnames)
                 SEXP R_xcentered, // centered design matrix
                 SEXP R_y, // response vector
                 SEXP R_fpmaxs, // vector of maximum fp degrees
                 SEXP R_fppos, // corresponding vector of fp column indices
                 SEXP R_fpcards, // corresponding vector of power set cardinalities
                 SEXP R_nFps, // number of fp terms
                 SEXP R_fpnames, // names of fp terms
                 SEXP R_ucIndices, // vector giving _unc_ertainty custer indices (column -> which group)
                 SEXP R_ucTermList, // list for group -> which columns mapping
                 SEXP R_nUcGroups, // number of uncertainty groups
                 SEXP R_hyperparam, // hyperparameter a for hyper-g prior
                 SEXP R_priorType, // use sparse model prior?
                 SEXP R_nModels, // number of best models to be returned
                 SEXP R_verbose, // should progress been displayed?
                 SEXP R_chainlength, // how many times should a jump been made?
                 SEXP R_nCache) // size of models cache (an STL map)
{
	// important!!! We now assume that all elements of R_fpmaxs are identical!!!
	// It would be best to remove the option supporting different maximum FP degrees from the code,
	// to be inline with the paper.

	// unpack ###
	// data
	const Matrix x = getMatrix(R_x);
	const Matrix xcentered = getMatrix(R_xcentered);
	const ColumnVector y = vec2col(R_y);

	dataValues data(x, xcentered, y, 0); // totalNumber is not needed

	// fp info
	fpInfo currentFpInfo(R_nFps,
	                     R_fpcards,
	                     R_fppos,
	                     R_fpmaxs,
	                     R_fpnames,
	                     x);

	// the FP range
	const std::set<unsigned int> fpRange = constructSequence(currentFpInfo.nFps);


	// uc info
	const int* ucIndicesArray = INTEGER(R_ucIndices);
	const vector<int> ucIndices(ucIndicesArray, ucIndicesArray + Rf_length(R_ucIndices));
	const int nUcGroups = INTEGER(R_nUcGroups)[0];

	vector<unsigned int> ucSizes;
	for (int i = 1; i <= nUcGroups; i++){
		ucSizes.push_back(count(ucIndices.begin(), ucIndices.end(), i));
	}
	int maxUcDim = accumulate(ucSizes.begin(), ucSizes.end(), 0);


	vector<IntSet> ucTermList(nUcGroups); // Array with length nUcGroups of int-vectors
	if(nUcGroups){ // catch case with no uc groups
		for(R_len_t i = 0; i != Rf_length(R_ucTermList); i++){
			SEXP temp = VECTOR_ELT(R_ucTermList, i);
			copy(INTEGER(temp), INTEGER(temp) + Rf_length(temp),
			     inserter(ucTermList[i], ucTermList[i].begin()));
		}
	}

	// determine columns that are always in the design matrix
	set<int> nonfixedCols = set<int>(currentFpInfo.fppos, currentFpInfo.fppos + currentFpInfo.nFps); // fixed info
	set<int> allCols, fixedCols;
	for(vector<int>::size_type i = 0; i != ucIndices.size(); i++){
		allCols.insert(i + 1);
		if (ucIndices.at(i))
			nonfixedCols.insert(i + 1);
	}
	set_difference(allCols.begin(), allCols.end(), nonfixedCols.begin(), nonfixedCols.end(),
				   inserter(fixedCols, fixedCols.begin()));
	unsigned int fixedDim = fixedCols.size();

	// now fixedCols contains indices of columns that are always present in the design matrix

	// prior specifications
	const double hyperparam = Rf_asReal(R_hyperparam);
        const std::string priorType = getStringVector(R_priorType).at(0);

        hyperPriorPars hyp(hyperparam,
                           priorType);

	// models which can be found during chain run can be cached in here:
	ModelCache modelCache(Rf_asInteger(R_nCache));

	// bookkeeping:
	book bookkeep; // 0) initializes empty sum of prop to posteriors and modelCounter 0

	// a) length of chain
	double chainlength = REAL(R_chainlength)[0];
	if (ULONG_MAX < chainlength){
		Rf_warning("\nchainlength too high - reducing to %d \n", ULONG_MAX);
		bookkeep.chainlength = ULONG_MAX;
	} else {
		bookkeep.chainlength = static_cast<PosLargeInt>(chainlength);
	}


	// b) verbose?
	bookkeep.verbose = LOGICAL(R_verbose)[0];

	// how many models to return?
	bookkeep.nModels = INTEGER(R_nModels)[0];

	// upper limit for num of columns
	unsigned int maxDim = min(static_cast<unsigned int>(data.nObs), fixedDim + currentFpInfo.maxFpDim + maxUcDim);

	// theta and theta':
	modelmcmc old, now;

	// start model
	modelPar startModel(currentFpInfo.nFps, 0, 0);
	PowersVector startFps(currentFpInfo.nFps); // initialize empty vector of correct length
	startModel.fpPars = startFps;
	old.modPar = startModel;
	old.dim = fixedDim;
	old.freeUcs = getFreeUcs(old.modPar, ucSizes, old.dim, maxDim);
	old.freeCovs = getFreeCovs(old.modPar, currentFpInfo, old.freeUcs, old.dim, maxDim);
	old.presentCovs = getPresentCovs(old.modPar);
	old.birthprob = 1; old.deathprob = old.moveprob = 0;

	Matrix oldDesign = getDesignMatrix(old.modPar, data, currentFpInfo, ucTermList, nUcGroups, fixedCols);
	double oldR2 = getR2(oldDesign, data, fixedCols, hyp);

	// log marginal likelihood
	old.logMargLik = getVarLogMargLik(oldR2, data.nObs, oldDesign.Ncols(), hyp);

	// log prior
        old.logPrior = getVarLogPrior(old.modPar, currentFpInfo, nUcGroups, hyp);

	// posterior expected g
	double oldPostExpectedg = posteriorExpectedg_hyperg(oldR2, data.nObs, oldDesign.Ncols(), hyp.a, old.logMargLik);

	// posterior expected shrinkage
	double oldPostExpectedShrinkage = posteriorExpectedShrinkage_hyperg(oldR2, data.nObs, oldDesign.Ncols(), hyp.a, old.logMargLik);

	// insert this model into map container
	modelInfo startInfo(old.logMargLik, old.logPrior, oldPostExpectedg, oldPostExpectedShrinkage, oldR2, 1);
	modelCache.insert(old.modPar, startInfo);

	// start with this model config
	now = old;

	// Start MCMC sampler***********************************************************//
	GetRNGstate(); // use R's random number generator
	for(PosLargeInt t = 0; t != bookkeep.chainlength; /* ++t explicitly at the end */){
		double logPropRatio; // log proposal ratio
		// randomly select move type
		double u1 = unif_rand();
		if (u1 < old.birthprob){											// BIRTH
			unsigned int newCovInd = discreteUniform<unsigned int>(old.freeCovs);
			if (newCovInd <= currentFpInfo.nFps){ 					// some fp index
				int powerIndex = discreteUniform(0, currentFpInfo.fpcards[newCovInd-1]);
				now.modPar.fpPars.at(newCovInd-1).insert(powerIndex);
				now.modPar.fpSize++; // correct invariants
				now.dim++;
				unsigned int newPowersEqualPowerIndex = count(now.modPar.fpPars.at(newCovInd-1).begin(), now.modPar.fpPars.at(newCovInd-1).end(), powerIndex);
				unsigned int m = old.modPar.fpPars.at(newCovInd-1).size();
				logPropRatio = log(static_cast<double>(newPowersEqualPowerIndex)) +
				              log(static_cast<double>(currentFpInfo.fpcards[newCovInd-1])) -
				              log1p(static_cast<double>(m));
			} else { 													// uc index
				int index = discreteUniform<int>(old.freeUcs);
				now.modPar.ucPars.insert(index);
				now.modPar.ucSize++;
				now.dim += ucSizes.at(index - 1);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				logPropRatio = log(static_cast<double>(old.freeUcs.size())) -
				        log(static_cast<double>(now.modPar.ucSize));
			}
			now.presentCovs.insert(newCovInd);
			now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
			if (now.dim == maxDim){
				now.birthprob = 0; now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 1.0 / 3 : 0.5;
			} else {
				now.birthprob = now.deathprob =	now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3;
			}
			logPropRatio += log(now.deathprob) - log(old.birthprob) +
			            log(static_cast<double>(old.freeCovs.size())) -
			                log(static_cast<double>(now.presentCovs.size()));
		} else if (u1 < old.birthprob + old.deathprob){					// DEATH
			unsigned int oldCovInd = discreteUniform<unsigned int>(old.presentCovs);
			if (oldCovInd <= currentFpInfo.nFps){ 					// some fp index
				Powers::iterator powerIterator = dU<Powers >(now.modPar.fpPars.at(oldCovInd-1));
				unsigned int oldPowersEqualPowerIndex = count(old.modPar.fpPars.at(oldCovInd-1).begin(), old.modPar.fpPars.at(oldCovInd-1).end(), *powerIterator);
				now.modPar.fpPars.at(oldCovInd-1).erase(powerIterator);
				now.modPar.fpSize--; // correct invariants
				now.dim--;
				logPropRatio = - log(static_cast<double>(oldPowersEqualPowerIndex)) -
				              log(static_cast<double>(currentFpInfo.fpcards[oldCovInd-1])) +
				                  log(static_cast<double>(old.modPar.fpPars.at(oldCovInd-1).size()));
			} else { 													// uc index
				set<int>::iterator IndIterator = dU<set<int> >(now.modPar.ucPars);
				now.modPar.ucSize--;
				now.dim -= ucSizes.at(*IndIterator - 1);
				now.modPar.ucPars.erase(IndIterator);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				logPropRatio = log(static_cast<double>(old.modPar.ucSize)) -
				        log(static_cast<double>(now.freeUcs.size()));
			}
			now.presentCovs = getPresentCovs(now.modPar);
			now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
			if (now.dim == fixedDim){
				now.birthprob = 1; now.deathprob = now.moveprob = 0;
			} else {
				now.birthprob = now.deathprob =	now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3;
			}
			logPropRatio += log(now.birthprob) - log(old.deathprob) +
			            log(static_cast<double>(old.presentCovs.size())) -
			                log(static_cast<double>(now.freeCovs.size()));

		} else if (u1 < old.birthprob + old.deathprob + old.moveprob){	 // MOVE
			unsigned int CovInd = discreteUniform<unsigned int>(old.presentCovs);
			if (CovInd <= currentFpInfo.nFps){ 						// some fp index
				Powers::iterator powerIterator = dU<Powers >(now.modPar.fpPars.at(CovInd-1));
				unsigned int oldPowersEqualPowerIndex = count(old.modPar.fpPars.at(CovInd-1).begin(), old.modPar.fpPars.at(CovInd-1).end(), *powerIterator);
				now.modPar.fpPars.at(CovInd-1).erase(powerIterator);
				int powerIndex = discreteUniform(0, currentFpInfo.fpcards[CovInd-1]);
				now.modPar.fpPars.at(CovInd-1).insert(powerIndex);
				unsigned int newPowersEqualPowerIndex = count(now.modPar.fpPars.at(CovInd-1).begin(), now.modPar.fpPars.at(CovInd-1).end(), powerIndex);
				// free, present Covs and move type probs are unchanged
				logPropRatio = log(static_cast<double>(newPowersEqualPowerIndex)) -
				        log(static_cast<double>(oldPowersEqualPowerIndex));

			} else { 													// uc index
				set<int>::iterator IndIterator = dU<set<int> >(now.modPar.ucPars);
				now.modPar.ucSize--;
				now.dim -= ucSizes.at(*IndIterator - 1);
				now.modPar.ucPars.erase(IndIterator);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				int index = discreteUniform<int>(now.freeUcs);
				now.modPar.ucPars.insert(index);
				now.modPar.ucSize++;
				now.dim += ucSizes.at(index - 1);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				// here something may change, therefore:
				now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
				if (now.dim == maxDim){
					now.birthprob = 0; now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 1.0 / 3 : 0.5;
				} else {
					now.birthprob = now.deathprob =	now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3;
				}
				logPropRatio = 0.0;
			}
		} else {													// SWITCH (of FP vectors)
			// select only the FP present covs
			std::set<unsigned int> presentFps = removeElement(old.presentCovs, currentFpInfo.nFps + 1);

			// so we have the first power vector:
			unsigned int firstFpInd = discreteUniform<unsigned int>(presentFps);
			Powers first = now.modPar.fpPars.at(firstFpInd - 1);

			// the second power vector from all other FPs
			std::set<unsigned int> otherFps = removeElement(fpRange, firstFpInd);
			unsigned int secondFpInd = discreteUniform<unsigned int>(otherFps);
			Powers second = now.modPar.fpPars.at(secondFpInd - 1);

			// save the first
			Powers saveFirst = first;

			// copy second to first
			now.modPar.fpPars.at(firstFpInd - 1) = second;

			// and save to second
			now.modPar.fpPars.at(secondFpInd - 1) = saveFirst;

			// so now we have switched the power vectors.

			// move type probs are not changed, because the number of present FPs is unchanged,
			// as well as the dimension of the model.

			// but carefully update the information which covariates are free and which are present
			now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
			now.presentCovs = getPresentCovs(now.modPar);

			// and the proposal ratio is 1, thus the log proposal ratio is 0:
			logPropRatio = 0;
		}

		// search for log marg lik of proposed model
		modelInfo nowInfo = modelCache.getModelInfo(now.modPar);

		if (R_IsNA(nowInfo.logMargLik))
		{ // "now" is a new model

		    // construct design matrix and compute R^2
		    Matrix nowDesign =
		            getDesignMatrix(now.modPar, data, currentFpInfo,
		                            ucTermList, nUcGroups, fixedCols);
		    double nowR2 = getR2(nowDesign, data, fixedCols, hyp);

		    if (R_IsNaN(nowR2))
		    { // check if new model is OK, if not then nan
		        now.logMargLik = R_NaN;

		        // we do not save this model in the model cache
		        bookkeep.nanCounter++;
		    }
		    else
		    { // OK: then compute the rest, and insert into model cache

		        // log marginal likelihood and log Bayes factor
		        now.logMargLik = getVarLogMargLik(nowR2, data.nObs,
		                                          nowDesign.Ncols(), hyp);

		        now.logPrior = getVarLogPrior(now.modPar, currentFpInfo, nUcGroups, hyp);

		        // posterior expected g
		        double nowPostExpectedg =
		                posteriorExpectedg_hyperg(nowR2, data.nObs,
		                                          nowDesign.Ncols(), hyp.a,
		                                          now.logMargLik);

		        // posterior expected shrinkage
		        double nowPostExpectedShrinkage =
		                posteriorExpectedShrinkage_hyperg(nowR2,
		                                                  data.nObs,
		                                                  nowDesign.Ncols(),
		                                                  hyp.a,
		                                                  now.logMargLik);

		        // insert the model parameter/info into the model cache

		        // problem: this could erase the old model from the model cache,
		        // and invalidate the iterator old.mapPos!
		        // ==> so we cannot work with the iterators here.
		        modelCache.insert(now.modPar,
		                          modelInfo(now.logMargLik, now.logPrior,
		                                    nowPostExpectedg,
		                                    nowPostExpectedShrinkage,
		                                    nowR2,
		                                    0));
		    }
		}
		else // "now" is an old model
		{
		    // extract log marg lik and prior from the modelInfo object
		    now.logMargLik = nowInfo.logMargLik;
		    now.logPrior = nowInfo.logPrior;
		}


		// decide acceptance:
		// for acceptance, the new model must be valid and the acceptance must be sampled
                if ((R_IsNaN(now.logMargLik) == FALSE) &&
                    (unif_rand() <= exp(now.logMargLik - old.logMargLik + now.logPrior - old.logPrior + logPropRatio)))
                { // acceptance
                    old = now;
                }
                else
                { // rejection
                    now = old;
                }

                // so now definitely old == now, and we can
                // increment the associated sampling frequency.
                modelCache.incrementFrequency(now.modPar);

                // echo progress?
		if((++t % max(bookkeep.chainlength / 100, static_cast<PosLargeInt>(1)) == 0) &&
		    bookkeep.verbose)
		{
			Rprintf("-"); // display computation progress at each percent
		}
	}
	PutRNGstate(); // no RNs required anymore


	// normalize posterior probabilities and correct log marg lik
	const long double logNormConst = modelCache.getLogNormConstant();
	const double logMargLikConst = - (data.nObs - 1) / 2.0 * log(data.sumOfSquaresTotal)  - log(hyp.a - 2.0);

	// get the nModels best models from the cache as an R list
	SEXP ret;
	Rf_protect(ret = modelCache.getListOfBestModels(currentFpInfo,
	                                                logMargLikConst,
	                                                logNormConst,
	                                                bookkeep));

	// set the attributes
	Rf_setAttrib(ret, Rf_install("numVisited"), Rf_ScalarReal(modelCache.size()));
	Rf_setAttrib(ret, Rf_install("inclusionProbs"), putDoubleVector(modelCache.getInclusionProbs(logNormConst, currentFpInfo.nFps, nUcGroups)));
	Rf_setAttrib(ret, Rf_install("linearInclusionProbs"), putDoubleVector(modelCache.getLinearInclusionProbs(logNormConst, currentFpInfo.nFps)));
	Rf_setAttrib(ret, Rf_install("logNormConst"), Rf_ScalarReal(logNormConst));

	if (bookkeep.verbose){
	    Rprintf("\nNumber of non-identifiable model proposals:     %d", bookkeep.nanCounter);
	    Rprintf("\nNumber of total cached models:                  %d", modelCache.size());
	    Rprintf("\nNumber of returned models:                      %d\n", Rf_length(ret));
	}


	// return ###
	Rf_unprotect(1);
	return ret;
}

// ***************************************************************************************************//

set<int> getFreeUcs(	// compute set of free uc group indices
					const modelPar& mod,
					const vector<unsigned int>& ucSizes,
					const unsigned int& currDim,
					const unsigned int& maxDim
					)
{
	set<int> ret;
	for (int i = 1; i <= static_cast<int>(ucSizes.size()); i++){ // for every uc index
		if ((find(mod.ucPars.begin(), mod.ucPars.end(), i) == mod.ucPars.end()) && (ucSizes.at(i-1) <= maxDim - currDim))
			ret.insert(i); // insert if not already in model and enough space in design matrix
	}
	return ret;
}


// ***************************************************************************************************//

set<unsigned int> getFreeCovs(					// compute set of free cov indices
				const modelPar& mod,
				const fpInfo& currFp,
				const set<int>& freeUcs,
				const unsigned int& currDim,
				const unsigned int& maxDim
					)
{
	set<unsigned int> ret;

	if (currDim == maxDim)
		return ret;

	for (unsigned int i = 0; i != mod.nFps; i++){
		if (mod.fpPars.at(i).size() < static_cast<PosInt>(currFp.fpmaxs[i]))
			ret.insert(i + 1);
	}

	if (! freeUcs.empty())
		ret.insert(mod.nFps + 1);

	return ret;
}

// ***************************************************************************************************//


set<unsigned int> getPresentCovs	( // determine set of present cov indices
					const modelPar& mod
						)
{
	set<unsigned int> ret;
	for (unsigned int i = 0; i != mod.nFps; i++){
		if (! mod.fpPars.at(i).empty())
			ret.insert(i+1);
	}
	if (! mod.ucPars.empty())
		ret.insert(mod.nFps + 1);
	return ret;
}

// ***************************************************************************************************//


template <class T>
T discreteUniform (	// return random element of myset; should be enclosed in getRNGstate() etc.
				const set<T>& myset
					)
{
	if (myset.empty())
		Rf_error("\nmyset is empty!\n");

	double u = unif_rand();
	typename set<T>::size_type size = myset.size();
	typename set<T>::const_iterator i = myset.begin(); typename set<T>::size_type j = 1;
	while(u > 1.0 / size * j){
		i++; j++;
	}
	return *i;
}

// ***************************************************************************************************//

template <class T>
typename T::iterator dU (	// return iterator of random element of myset; should be enclosed in getRNGstate() etc.
				const T& container
					)
{
	if (container.empty())
		Rf_error("\ncontainer is empty!\n");

	double u = unif_rand();
	typename T::size_type size = container.size();
	typename T::iterator i = container.begin(); typename T::size_type j = 1;
	while(u > 1.0 / size * j){
		i++; j++;
	}
	return i;
}


// ***************************************************************************************************//

int discreteUniform ( // get random int x with lower <= x < upper; should be enclosed in getRNGstate() etc.
						const int& lower,
						const int& upper
					)
{
	if (lower >= upper)
		Rf_error("\nlower = %d >= %d = upper in discreteUniform call\n", lower, upper);

	int size = upper - lower;
	int ret = lower;
	double u = unif_rand();

	while(u > 1.0 / size * (ret - lower + 1)){
		ret++;
	}
	return ret;
}

// ***************************************************************************************************//

ReturnMatrix getDesignMatrix ( // construct centered design matrix including intercept for the model
					const modelPar &mod,
					const dataValues &data,
					const fpInfo &currFp,
					const vector<IntSet>& ucTermList,
					const int &nUcGroups,
					const set<int> &fixedCols
					)
{
//
//	MATRIXSTORE(data.response, responseStore)
//
	// build design matrix B

	// intercept column
	Matrix B = data.onesVector;

	// centered fp matrices
	for (unsigned int i = 0; i != currFp.nFps; i++){
		Powers powersi = mod.fpPars.at(i);
		if (! powersi.empty()){
			Matrix Fp = getFpMatrix(currFp.tcols.at(i), powersi, data); // this is centered
			B = B | Fp;
		}
	}

	// centered uc matrices
	for (int i = 0; i != nUcGroups; i++){
		set<int>::const_iterator ipos = find(mod.ucPars.begin(), mod.ucPars.end(), i + 1);
		if (ipos != mod.ucPars.end()){ // if mod.ucPars contains i
			Matrix Uc = getMultipleCols(data.centeredDesign, ucTermList.at(i));	// this is centered
			B = B | Uc;
		}
	}

	B.release();
	return B;
}




// ***************************************************************************************************//

double getR2(	// compute coefficient of determination for the model
			const Matrix &design, // must be the centered design matrix!! (including the intercept)
			const dataValues &data,
			const set<int> &fixedCols,
			const hyperPriorPars &hyp
			)
{

	int dim = design.Ncols();
	if (dim - 1 >= data.nObs - 3 - hyp.a) return R_NaN; // not a valid model
//
//	MATRIXSTORE(B, BStore)
//
	const int numFixedCols = fixedCols.size();

	if (dim == numFixedCols) { // then this is the null model
		return 0; // because SSE == SST in this case
	}
	else {	// else we need to work a bit

		// select nonfixed part of the centered design matrix, because the fixed part is not needed here
		Matrix X = design.Columns(numFixedCols + 1, dim);

//		double test;
//		MATRIXSTORE(X, XStore)

		// compute model-specific parts
		SymmetricMatrix XtX;
		XtX << X.t() * X;

		try // a cholesky decomposition of XtX
		{
		    LowerTriangularMatrix LeftRootOfXtX = Cholesky(XtX);

//		    MATRIXSTORE(LeftRootOfXtX, lroxtxStore)

		    // compute coefficient of determination R2
		    ColumnVector tmp = LeftRootOfXtX.i() * (X.t() * data.response);

//		    MATRIXSTORE(tmp, tmpStore)

		    double sumOfSquaresModel = tmp.sum_square();

			double R2 = sumOfSquaresModel / data.sumOfSquaresTotal;

			assert((R2 <= 1) & (R2 >= 0));

			return R2;
		}
		catch(NPDException) {return R_NaN;} // if XtX is not p.d. then return NAN
	}
}

// ***************************************************************************************************//

// compute varying part of log marginal likelihood for specific model
double getVarLogMargLik(const double &R2, const int &n, const int &dim, const hyperPriorPars &hyp)
{
    // check if any interrupt signals have been entered
    // (check here because this function is used both by the exhaustive and the sampling function)
    R_CheckUserInterrupt();

    // then start computing
	if(dim == 1){
		return 0;
	} else {

		double logBF = logBF_hyperg(R2, n, dim, hyp.a);
		assert(! ISNAN(logBF)); // each possible case has been treated above!

		return logBF;
	}
}

// ***************************************************************************************************//

double getVarLogPrior( // compute varying part of logarithm of model prior
                      const modelPar &mod,
                      const fpInfo &currFp,
                      const PosInt nUcGroups,
                      const hyperPriorPars &hyp)
{
    if (hyp.priorType == "sparse")
    {
        safeSum thisVarLogPrior;
        for (unsigned int i = 0; i != currFp.nFps; i++)
        { // for each fp covariate
            unsigned int degree = mod.fpPars.at(i).size();
            double thisVal = -Rf_lchoose(currFp.fpcards[i] - 1 + degree, degree)
                    - log1p(currFp.fpmaxs[i]);
            thisVarLogPrior.add(thisVal);
        }
        return thisVarLogPrior.sum() - (nUcGroups * M_LN2);
    }
    else if (hyp.priorType == "dependent")
    {
        // determine number of all covariates (covariate groups):
        int nCovs = nUcGroups + currFp.nFps;

        // determine number of included FPs and which are nonlinear:
        int nInclContinuous = 0;
        PosIntVector nonlinearFps;

        for(PosInt i = 0; i != currFp.nFps; i++)
        {
            Powers powersi = mod.fpPars.at(i);
            if (! powersi.empty())
            {
                ++nInclContinuous;

                if(mod.fpPars.at(i) != currFp.linearPowers)
                {
                    nonlinearFps.push_back(i);
                }
            }
        }

        // determine number of included discrete covariates:
        int nInclDiscrete = mod.ucSize;

        // so altogether there are
        int nIncluded = nInclContinuous + nInclDiscrete;
        // included covariates

        // and the number of possible nonlinear transformations
        // for each variable is also included in the computations:
        double sumLogNonlinearPossibilities = 0.0;
        for(PosIntVector::const_iterator
                i = nonlinearFps.begin();
                i != nonlinearFps.end();
                ++i)
        {
            sumLogNonlinearPossibilities += log(currFp.numberPossibleFps.at(*i) - 2.0);
            //                              Note: degree 0 and linear degree 1 FP are subtracted
        }

        double result = - log1p(nCovs) - Rf_lchoose(nCovs, nIncluded) -
                log1p(nInclContinuous) - Rf_lchoose(nInclContinuous, nonlinearFps.size()) -
                sumLogNonlinearPossibilities;

        return result;
    }
    else // if (hyp.priorType == "flat")
    {
        return - (nUcGroups * M_LN2);
    }
}

// ***************************************************************************************************//

ReturnMatrix getFpMatrix( // build Fp basis matrix from vector, power indices and power set for one covariate
			const vector<ColumnVector> &tcols,
			const Powers &powerinds,
			const dataValues &data
			)
{
	const int logInd = 3; // this index corresponds to power = 0, i.e. log.
	const int nrow = tcols.at(0).Nrows();

	Matrix ret(nrow, powerinds.size());

	// start recursion
	int lastInd = logInd;
	ColumnVector lastCol(nrow); lastCol = 1;

	// there is at least one power present
	Powers::size_type cols = 1; // invariant: about to process column number cols

	for (Powers::const_iterator now = powerinds.begin(); now != powerinds.end(); now++){
		if (*now == lastInd){ 	// repeated powers case
			lastCol = SP(lastCol, tcols.at(logInd));
		} else {				// normal case
			lastInd = *now;
			lastCol = tcols.at(lastInd);
		}
		// center the column
		ret.Column(cols++) = lastCol - (lastCol.sum() / data.nObs) * data.onesVector;
	}

	ret.Release(); return ret;
}

// ***************************************************************************************************//



void pushInclusionProbs(	// push back index into covGroupWisePosteriors-Array
						const modelPar &mod,
						const fpInfo &currFp,
						const int &nUcGroups,
						book &bookkeep)
{
    for (PosInt i = 0; i != currFp.nFps; i++){
        if (! mod.fpPars.at(i).empty())
        {
            bookkeep.covGroupWisePosteriors.at(i).add(bookkeep.modelCounter);

            // also record if this FP is just a linear effect
            if(mod.fpPars.at(i) == currFp.linearPowers)
            {
                bookkeep.linearFpPosteriors.at(i).add(bookkeep.modelCounter);
            }
        }
    }

    for (int i = 1; i <= nUcGroups; i++){
        set<int>::const_iterator ipos = find(mod.ucPars.begin(), mod.ucPars.end(), i);
        if (ipos != mod.ucPars.end()) // if mod.ucPars contains i
        {
            bookkeep.covGroupWisePosteriors.at(i - 1 + currFp.nFps).add(bookkeep.modelCounter);
        }
    }
}



// ***************************************************************************************************//

// this is an interface for R to the log marg lik computation
SEXP logMargLik( //definition
				SEXP R_R2,				// coefficient of determination
				SEXP R_n,				// number of observations
				SEXP R_dim,				// number of columns of the design matrix
				SEXP R_alpha, 			// hyperparamater for hyper-g prior
				SEXP R_sst				// total sum of squares computed from y
				)
{
	unsigned int nProtect = 0;

	// unpack
	const double R2 = REAL(R_R2)[0];
	const int n = INTEGER(R_n)[0];
	const int dim = INTEGER(R_dim)[0];
	const double alpha = REAL(R_alpha)[0];
	const double sst = REAL(R_sst)[0];

	// compute
	hyperPriorPars hyp(alpha, std::string("flat")); // prior type does not matter here
	double varLogMargLik = getVarLogMargLik(R2, n, dim, hyp);
	double logMargLikConst = - (n - 1) / 2.0 * log(sst)  - log(hyp.a - 2.0);

	SEXP ret;
	Rf_protect(ret = Rf_ScalarReal(varLogMargLik + logMargLikConst));
	nProtect++;

	Rf_unprotect(nProtect);
	return(ret);
}

// ***************************************************************************************************//

// this is an interface for R to the computation of the posterior expected g
SEXP postExpectedg( // definition
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha) // hyperparamater for hyper-g prior
{
    unsigned int nProtect = 0;

    // unpack
    const double R2 = REAL(R_R2)[0];
    const int n = INTEGER(R_n)[0];
    const int dim = INTEGER(R_dim)[0];
    const double alpha = REAL(R_alpha)[0];

    // compute
    hyperPriorPars hyp(alpha, std::string("flat")); // value of priorType does not matter here
    const double varLogMargLik = getVarLogMargLik(R2, n, dim, hyp);
    const double postExpectedg = posteriorExpectedg_hyperg(R2, n, dim, hyp.a, varLogMargLik);

    SEXP ret;
    Rf_protect(ret = Rf_ScalarReal(postExpectedg));
    nProtect++;

    Rf_unprotect(nProtect);
    return(ret);
}

// ***************************************************************************************************//

// this is an interface for R to the computation of the posterior expected shrinkage factor
SEXP postExpectedShrinkage( //definition
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha) // hyperparamater for hyper-g prior
{
    unsigned int nProtect = 0;

    // unpack
    const double R2 = REAL(R_R2)[0];
    const int n = INTEGER(R_n)[0];
    const int dim = INTEGER(R_dim)[0];
    const double alpha = REAL(R_alpha)[0];

    // compute
    hyperPriorPars hyp(alpha, std::string("flat")); // value of priorType does not matter here
    const double varLogMargLik = getVarLogMargLik(R2, n, dim, hyp);
    const double postExpectedShrinkage = posteriorExpectedShrinkage_hyperg(R2, n, dim, hyp.a, varLogMargLik);

    SEXP ret;
    Rf_protect(ret = Rf_ScalarReal(postExpectedShrinkage));
    nProtect++;

    Rf_unprotect(nProtect);
    return(ret);
}

// ################################################################################################

int main() {} // dummy


// attention, when building shared library with R CMD SHLIB enter ALL source files that are entangled with
// this, too!! (RnewMat.cpp e.g.)
