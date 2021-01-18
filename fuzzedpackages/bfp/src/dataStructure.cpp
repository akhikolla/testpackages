#include "dataStructure.h"
#include "RnewMat.h"
#include "sum.h"
#include <set>
#include <types.h>

#include <cassert>
#include <algorithm>
#include <numeric>
#include <rcppExport.h>

using std::lexicographical_compare;
using std::pair;
using std::map;
using std::set;
using std::min;
using std::accumulate;
using std::max;
using std::max_element;

using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::_;

// model info functions
modelInfo& modelInfo::operator=(const modelInfo& m)
{ // assignment operator 
	if (this != &m){
		logMargLik = m.logMargLik;
		logPrior = m.logPrior;
		logPost = m.logPost;
		postExpectedg = m.postExpectedg;
		postExpectedShrinkage = m.postExpectedShrinkage;		
		hits = m.hits;
		R2 = m.R2;
	}
	return *this;	
}
	
List
modelInfo::convert2list(double addLogMargLikConst,
                        long double logNormConst,
                        const book& bookkeep) const
{
    return List::create(_["logM"] = logMargLik + addLogMargLikConst,
                        _["logP"] = logPrior,
                        _["posterior"] = NumericVector::create(exp(logPost - logNormConst),
                                                               hits * 1.0 / bookkeep.chainlength),
                        _["postExpectedg"] = postExpectedg,
                        _["postExpectedShrinkage"] = postExpectedShrinkage,
                        _["R2"] = R2);
}

// model parameter functions
	
bool modelPar::operator<(const modelPar& m) const
{ 
	// return size() < m.size(); way too easy... lexicographical comparison, starting with uc indices:
	if (ucPars < m.ucPars)
		return true;
	else if (ucPars > m.ucPars)
		return false;
	else // uc indices are equal
		return lexicographical_compare(fpPars.begin(), fpPars.end(), m.fpPars.begin(), m.fpPars.end());	
}

modelPar& modelPar::operator=(const modelPar& m)
{
	if (this != &m){
		fpPars = m.fpPars;
		ucPars = m.ucPars;
		fpSize = m.fpSize;
		ucSize = m.ucSize;
		nFps = m.nFps;
	}
	return *this;	
}

int modelPar::size() const
{
	return fpSize + ucSize;
}

List
modelPar::convert2list(const fpInfo& currFp) const
{
    // powers
    List powers(nFps);
    powers.names() = currFp.fpnames;

    for (PosInt i = 0; i != nFps; ++i)
    {
        powers[i] = currFp.inds2powers(fpPars[i]);
    }

    // return with uc settings
    return List::create(_["powers"] = powers,
                        _["ucTerms"] = ucPars);
}


// internal helper function to combine the contents of
// two named R lists into one named R list
List
combineLists(List firstList, List secondList)
{
    // allocate result list
    List ret(firstList.length() + secondList.length());

    // collect the names
    StringVector names;

    StringVector firstNames = firstList.names();
    StringVector secondNames = secondList.names();

    // now fill in contents of the first list
    for(R_len_t i = 0; i < firstList.length(); ++i)
    {
        ret[i] = firstList[i];
        names.push_back(firstNames.at(i));
    }

    // now fill in contents of the second list
    for(R_len_t i = 0; i < secondList.length(); ++i)
    {
        ret[i + firstList.length()] = secondList[i];
        names.push_back(secondNames.at(i));
    }

    // assign names to list
    ret.names() = names;

    // unprotect and return
    return ret;
}


// for combination of modelPar and modelInfo: model //

model& model::operator=(const model& m) // assignment operator
{ 
	if (this != &m){
		info = m.info;
		par = m.par;		
	}
	return *this;
}

bool model::operator<(const model& m) const  // less		
{
	double thisLogPost = info.logMargLik + info.logPrior;
	double mLogPost = m.info.logMargLik + m.info.logPrior;
	if (thisLogPost < mLogPost)
		return true;
	else if (thisLogPost > mLogPost)
		return false;
	else  // posteriors are equal, then the parameter makes the decision
		return m.par < par;
}

SEXP model::convert2list(const fpInfo& currFp,
                         double addLogMargLikConst,
                         long double logNormConst,
                         const book& bookkeep) const // convert model into list for export to R
{
    return combineLists(par.convert2list(currFp),
                        info.convert2list(addLogMargLikConst,
                                          logNormConst,
                                          bookkeep));
}

// dataValues //

dataValues::dataValues(const Matrix &x, const Matrix &xcentered, const ColumnVector &y, const double &totalNum) : 
	design(x), centeredDesign(xcentered), response(y), totalNumber(static_cast<map<modelPar, modelInfo>::size_type>(totalNum)) 
{
	// number of observations
	nObs = design.Nrows();

	// nObs long vector of ones
	onesVector = ColumnVector(nObs);
	onesVector = 1.0;
	
	// and the SST
	ColumnVector centeredResponse = response - (response.sum() / nObs) * onesVector;
	sumOfSquaresTotal = centeredResponse.sum_square();   	
}	


// fpInfo //

fpInfo::fpInfo(SEXP R_nFps,
               SEXP R_fpcards,
               SEXP R_fppos,
               SEXP R_fpmaxs,
               SEXP R_fpnames,
               const Matrix& x) :
               nFps(INTEGER(R_nFps)[0]),
               fpcards(INTEGER(R_fpcards)),
               fppos(INTEGER(R_fppos)),
               fpmaxs(INTEGER(R_fpmaxs)),
               biggestMaxDegree(*max_element(fpmaxs,
                                             fpmaxs + Rf_length(R_fpmaxs))),
               maxFpDim(accumulate(fpmaxs, fpmaxs + nFps, 0)),
               powerset(max(8, 5 + biggestMaxDegree)),
               fpnames(R_fpnames),
               numberPossibleFps(),
               linearPowers(),
               tcols(nFps)
{
    // corresponding indices        0   1     2  3    4  5  6  7
    const double fixedpowers[] = { -2, -1, -0.5, 0, 0.5, 1, 2, 3 }; // always in powerset

    copy(fixedpowers, fixedpowers + 8,
         inserter(powerset, powerset.begin()));
    for(int more = 3; more < biggestMaxDegree; more++){ // additional powers
        powerset.at(8 + (3 - more)) = more + 1;
    }

    // numbers of possible univariate fps?
    for(PosInt i=0; i != nFps; ++i)
    {
        int thisNumber = 0;
        for(int deg = 0; deg <= fpmaxs[i]; ++deg)
        {
            thisNumber += Rf_choose(fpcards[i] - 1 + deg, deg);
        }
        numberPossibleFps.push_back(thisNumber);
    }

    // insert the index 5 for linear power 1
    linearPowers.insert(5);

    // build array of vectors of ColumnVectors holding the required
    // transformed values for the design matrices
    for (PosInt i = 0; i != nFps; i++){ // for every fp term
        const int nCols = fpcards[i];
        const ColumnVector thisCol = x.Column(fppos[i]);
        vector<ColumnVector> thisFp;
        for (int j = 0; j != nCols; j++){ // for every possible power
            ColumnVector thisTransform = thisCol;
            double thisPower = powerset.at(j);
            if(thisPower){ // not 0
                for (int k = 0; k != thisTransform.Nrows(); k++){ // transform each element
                    assert(thisTransform.element(k) > 0);
                    thisTransform.element(k) = pow(thisTransform.element(k), thisPower);
                    assert(! ISNAN(thisTransform.element(k)));
                }
            } else { // 0
                for (int k = 0; k != thisTransform.Nrows(); k++){
                    assert(thisTransform.element(k) > 0);
                    thisTransform.element(k) = log(thisTransform.element(k));
                    assert(! ISNAN(thisTransform.element(k)));
                }
            }
            // do not! center the column. This is done inside getFpMatrix, because
            // the repeated powers case cannot be treated here!!

            // and put it into vector of columns
            thisFp.push_back(thisTransform);
        }
        tcols.at(i) = thisFp;
    }
}


DoubleVector
fpInfo::inds2powers(const Powers& m) const // convert inds m (a "Powers" object) into powers vector p (a DoubleVector)
{
    DoubleVector ret;

    for (Powers::const_iterator j = m.begin();
            j != m.end();
            j++)
    {
        ret.push_back(powerset[*j]);
    }

    return ret;
}


// safeSum //

void safeSum::add(const long double &val) 
{
	vals.push_back(val);
}

long double safeSum::sum()
{
	long double ret	= modified_deflation(vals);
	return ret;
}

long double safeSum::logSumExp()
{
    // the maximum of the log contributions is:
    long double maxLogContrib = *std::max_element(vals.begin(), vals.end());

    // now compute the constant which is added to all log contributions,
    // in order to avoid infinite contributions and at the same time use
    // the whole number space (i.e. possibly avoid zero contributions)
    long double constant = log(LDBL_MAX) - 100.0L - maxLogContrib;
    // 100 is for safety.

    // so now the contributions, offset by the constant
    LongDoubleVector expVals;
    for(LongDoubleVector::const_iterator
            l = vals.begin();
            l != vals.end();
            ++l)
    {
        expVals.push_back(exp(*l + constant));
    }

    // the result is the log of the sum, corrected with the constant:
    long double ret = log(modified_deflation(expVals)) - constant;
    return ret;
}

long double safeSum::simpleSum()
{
        long double ret = 0.0;
        for(LongDoubleVector::const_iterator
                v = vals.begin();
                v != vals.end();
                ++v)
        {
            ret += *v;
        }
        return ret;
}


// indexSafeSum //

void indexSafeSum::add(const std::vector<long double>::size_type& ind)
{
	indices.insert(ind);	
}

long double indexSafeSum::sum(const safeSum& s) const
{
	vector<long double> tempVec;
	for(set<indexType>::const_iterator i = indices.begin(); i != indices.end(); i++){
		tempVec.push_back(s.vals.at(*i));	
	}
	return modified_deflation(tempVec);
}


// ModelCache //

// insert model parameter and corresponding info into cache,
// with caring about the maximum number of elements in the map.
bool
ModelCache::insert(const modelPar& par, const modelInfo& info)
{
    // first check size of cache
    if(isFull())
    {
        // if we are full, then check if this log posterior is better than
        // the worst cached model, which is pointed to by
        MapType::iterator worstModelIter = *(modelIterSet.begin());

        // the comparison
        if((worstModelIter->second.logPost) < info.logPost)
        {
            // new model is better than worst model cached.
            // so we delete the worst model from the cache.

            // first from the map
            modelMap.erase(worstModelIter);
            // and then from the set
            modelIterSet.erase(modelIterSet.begin());
        }
        else
        {
            // the new model is not better than the worst model cached,
            // so we do not cache it.
            return false;
        }
    }

    // so now we know that we want to insert the model into the cache,
    // either because the cache was not full or because the new model was better
    // than the worst model cached.

    // -> try inserting into the map:
    pair<MapType::iterator, bool> ret = modelMap.insert(MapType::value_type(par, info));

    // if we were successful:
    if(ret.second)
    {
        // then also insert the iterator pointing to the map element into the set.
        modelIterSet.insert(ret.first);

        // return success
        return true;
    }
    else
    {
        return false;
        Rf_error("Should not happen: model already contained in model cache!");
    }
}

// search for the model info of a model config in the map,
// and return an information with NA for log marg lik if not found
modelInfo
ModelCache::getModelInfo(const modelPar& par) const
{
    // search for the config in the map
    MapType::const_iterator ret = modelMap.find(par);

    // if found, return the log marg lik
    if(ret != modelMap.end())
        return ret->second;
    else
        return modelInfo(R_NaReal, R_NaReal, 0.0, 0.0, 0.0);
}

// increment the sampling frequency for a model configuration
// (of course, if this config is not cached nothing is done)
void
ModelCache::incrementFrequency(const modelPar& par)
{
    // search for the config in the map
    MapType::iterator ret = modelMap.find(par);

    // if found, increment the hits
    if(ret != modelMap.end())
        ret->second.hits++;
}

// compute the log normalising constant from all cached models
long double
ModelCache::getLogNormConstant() const
{
    // use safe summation
    safeSum vec;

    // traverse the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // and add all unnormalized log posteriors
        vec.add(m->second.logPost);
    }

    // return the log of the sum of the exp'ed saved elements
    return vec.logSumExp();
}

// compute the inclusion probabilities from all cached models,
// taking the log normalising constant, the number of FPs and the number of UC groups
DoubleVector
ModelCache::getInclusionProbs(long double logNormConstant, PosInt nFps, PosInt nUcs) const
{
    // abbreviation
    typedef std::vector<safeSum> SafeSumVector;
    // allocate vector of safeSum objects for all FPs
    SafeSumVector fps(nFps);

    // and all UC groups
    SafeSumVector ucs(nUcs);

    // now process each model in the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // abbrevs
        const modelPar& thisPar = m->first;
        const modelInfo& thisInfo = m->second;

        // first process the FPs
        {
        SafeSumVector::iterator s = fps.begin();
        for (PowersVector::const_iterator
                p = thisPar.fpPars.begin();
                p != thisPar.fpPars.end();
                ++p, ++s)
        {
            // is this FP in the model m?
            if (! p->empty())
            {
                // then add the normalized model probability onto his FP stack
                s->add(exp(thisInfo.logPost - logNormConstant));
            }
        }
        }

        // then process the UC groups
        {
        SafeSumVector::iterator s = ucs.begin();
        for (PosInt i = 1; i <= nUcs; ++i, ++s)
        {
            // is this UC group in the model m?
            if (thisPar.ucPars.find(i) != thisPar.ucPars.end())
            {
                // then add the normalized model probability onto his UC stack
                s->add(exp(thisInfo.logPost - logNormConstant));
            }
        }
        }
    } // end processing all models in the cache

    // so now we can sum up safesum-wise to the return double vector
    DoubleVector ret;

    for(SafeSumVector::iterator
            s = fps.begin();
            s != fps.end();
            ++s)
    {
        ret.push_back(s->sum());
    }

    for(SafeSumVector::iterator
            s = ucs.begin();
            s != ucs.end();
            ++s)
    {
        ret.push_back(s->sum());
    }

    return ret;
}

// compute the linear inclusion probabilities from all cached models,
// taking the log normalising constant and the number of FPs
DoubleVector
ModelCache::getLinearInclusionProbs(long double logNormConstant, PosInt nFps) const
{
    // abbreviation
    typedef std::vector<safeSum> SafeSumVector;
    // allocate vector of safeSum objects for all FPs
    SafeSumVector fps(nFps);

    // what is a "linear inclusion"?
    Powers linear;
    linear.insert(5);

    // now process each model in the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // abbrevs
        const modelPar& thisPar = m->first;
        const modelInfo& thisInfo = m->second;

        // first process the FPs
        {
        SafeSumVector::iterator s = fps.begin();
        for (PowersVector::const_iterator
                p = thisPar.fpPars.begin();
                p != thisPar.fpPars.end();
                ++p, ++s)
        {
            // is this FP linear?
            if (*p == linear)
            {
                // then add the normalized model probability onto his stack
                s->add(exp(thisInfo.logPost - logNormConstant));
            }
        }
        }
    } // end processing all models in the cache

    // so now we can sum up safesum-wise to the return double vector
    DoubleVector ret;

    for(SafeSumVector::iterator
            s = fps.begin();
            s != fps.end();
            ++s)
    {
        ret.push_back(s->sum());
    }

    return ret;
}



// convert the best nModels from the cache into an R list
List
ModelCache::getListOfBestModels(const fpInfo& currFp,
                                double addLogMargLikConst,
                                long double logNormConst,
                                const book& bookkeep) const
{
    // allocate the return list
    List ret(std::min(bookkeep.nModels,
                      static_cast<PosInt>(modelIterSet.size())));
    // cast is necessary for gcc-4.2 on Mac on R-forge.

    // process the ordered list of best models from the end (because the set is ordered increasingly)
    PosInt i = 0;
    for(SetType::const_reverse_iterator
            s = modelIterSet.rbegin();
            (i < bookkeep.nModels) && (s != modelIterSet.rend());  // so the return list has min(nModels, modelIterSet.size()) elements.
            ++s, ++i)
    {
        // and for this model, combine the config and info lists to one list and
        // put that in the i-th slot of the return list.
        ret[i] = combineLists((**s).first.convert2list(currFp),
                              (**s).second.convert2list(addLogMargLikConst,
                                                        logNormConst,
                                                        bookkeep));
    }

    // return
    return ret;
}
