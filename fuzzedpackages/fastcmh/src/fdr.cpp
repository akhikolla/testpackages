//based on Gilbert (2005) A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics

//Dean Bodenham June 2016

#include<iostream>
#include<vector>
//for sort
# include<algorithm>
//for log
#include<math.h>


//NOTE:
//'difficult' part was sorting the three vectors pvalue, tau and l according the the order in pvalue. Although it would have been relatively straightforward to use the sort function and a lambda function, for some reason the g++ compiler on MAc OSX (gcc 4.2) would not accept this. Then, after installing g++ 4.5, it created other errors. In the end, I decided to make a pair of pvalue and a perm vector (could have used iota), sort the pair, and then extract the permutation half of the pair. Later, access elements of the three vectors according to the permutation order.
//Based on method from: 
// http://stackoverflow.com/questions/236172/how-do-i-sort-a-stdvector-by-the-values-of-a-different-stdvector

//A struct for ordering by the second component of the pair
//Used for sorting in increasing order
//Later will run through list in REVERSE, from largest p-value to smallest
struct orderBySecond {
    bool operator ()(std::pair<long long, double> const& a, std::pair<long long, double> const& b) {
        return (a.second) < (b.second);
    }
};


//This function takes as input a vector of p-values and, 
//instead of returning the order of the p-values, returns the permutation
//of the indices that would sort the p-values into increasing order.
//
//This permutation index set is returned, because it is then used to extract
//pvalue, tau and l from three separate vectors (rather than sorting three 
//vectors).
std::vector<long long> extractPermutation(std::vector<double>& pvalue, std::vector<long long>& tau, std::vector<long long>& l){

    //create pairs
    std::vector<std::pair<long long, double> > perm(pvalue.size());
    long long n = 0;
    for (std::vector<double>::iterator it = pvalue.begin(); it != pvalue.end(); ++it, ++n){
        perm[n] = std::make_pair(n, *it);
    }

    //sort pairs
    sort(perm.begin(), perm.end(), orderBySecond());

    //extract perm
    std::vector<long long> onlyperm(perm.size());
    n = 0;
    for (std::vector<std::pair<long long, double> >::iterator pit = perm.begin(); pit != perm.end(); ++pit, ++n){
        onlyperm[n] = ((*pit).first);
    }
    return(onlyperm);
}


//This function compute the nth harmonic number
//this is needed when computing FDR for dependent p-values, 
//under the special case of dependence known as positive regression dependence
double computeHarmonic (long long n){
    //initialise to 1.0 to take care of case when n = 1
    double h = 1.0;
    //increase by 1 in order to use <=
    n += 1;
    for (long long i = 2; i < n; ++i){
        h += (double) 1/i;
    }
    return h;
}


//computes the lower bound for the nth harmonic number 
//using the approximation with the Euler-Mascheroni constant
//computes lower, because later will take inverse and we want approximation to be larger
double computeApproxHarmonicLower(long long n){
    const double GAMMA = 0.577215664901532;
    //to get rid of ambiguity error, casting long long to double
    //previously:
    //double lower = (double) 1/ (2*n + 1) + log(n) + GAMMA;
    //now:
    double nDouble = (double) n ;
    double lower = 1 / (2*nDouble + 1) + log(nDouble) + GAMMA;
    return lower;
}


//This function computes the nth harmonic number
//count is the number of terms in the p-value set
//either exactly (if count is less than 100) or using the approximation
//involving the Euler-Mascheroni constant; the approximation is very
//good for n >= 100
double computeHarmonicFast (long long count){
    double h = 1.0;
    const long long MAX_FOR_BRUTE_FORCE = 100;
    count += 1;

    if (count < MAX_FOR_BRUTE_FORCE){
        for (long long i = 2; i < count; ++i){
            h += (double) 1/i;
        }
    } else {
        h = computeApproxHarmonicLower(count);
    }
    return h;
}


//This function compute the FDR limit at time i.
//Here, q will either be alpha, or alpha divided by the harmonic mean.
double computeFDRLimit(double q, long long i, long long m){
    return ( ( (q * i) / m) );
}


//A function to compute the adjusted alpha, taking dependence into account
double computeAdjustedFDRAlpha(double alpha, long long m_K, bool useDependence){
    double alphaAdjusted = alpha;
    if (useDependence){
        alphaAdjusted = alpha / computeHarmonicFast(m_K);
    }
    return alphaAdjusted;
}

//This function applies the FDR-Tarone method described in Gilbert. 
//At the moment, it only implements the Benjamini-Yekutieli method
//which takes dependence into account.
std::vector<long long> gilbertFDR(std::vector<double>& pvalue, std::vector<long long>& tau, std::vector<long long>& l, double alpha, bool useDependence){

    //m_K is simply the size of the pvalue set, which has already passed
    //the minimum attainable p-value threshold
    long long m_K = pvalue.size();

    //compute necessary factor to take dependence into account
    //simply alpha divided by the harmonic factor 
//     double dependenceAlpha = alpha / computeHarmonicFast(m_K);
    double dependenceAlpha = computeAdjustedFDRAlpha(alpha, m_K, useDependence);

    std::vector<long long> perm = extractPermutation(pvalue, tau, l);

    //will keep iterating until run out of values or until k is found
    bool kNotFound = true;
    long long index = 0;
    long long maxCount = perm.size();
    long long count = maxCount-1;
    //this will be the k value
    long long kval = 0;
    double fdr_limit = 0;
    //this is the p-value
    double thisPval = 0;

//     std::cout << "depAlpha: " << dependenceAlpha << std::endl;

    //run through all pvalues in permutation order
    while ((kNotFound) && (count > 0) ){

        //decrease count: starts at maxCount, but index must be maxCount-1
        //Also, final value outside loop is count==1, then reduced to count=0
        count--;
        //obtain the pvalue via the permutation
        index = perm[count];
        thisPval = pvalue[index];

        //compute the fdr limit
        //VERY IMPORTANT: need to pass count + 1, because cpp counts from 0
        fdr_limit = computeFDRLimit(dependenceAlpha, count+1, m_K);
        if (thisPval <= fdr_limit){
            kNotFound = false;
            //since this is the first reject, the number of accepts is count
            //e.g. check case for count = 0, count = 1
            kval = count;
        }

//         std::cout << "count: " << count << ", thisPval: " << thisPval << ", this_fdr_limit: " << fdr_limit << std::endl;

    }

//     std::cout << "total number of values in perm: " << m_K << std::endl;

    //even if kval = 0, this is fine
    //now copy perm to fdr_perm
//     std::cout << " kval= " << kval << std::endl;
    std::vector<long long>::const_iterator permFirst = perm.begin();
    std::vector<long long>::const_iterator permLast = perm.begin() + kval + 1;
    std::vector<long long> fdr_perm(permFirst, permLast);
    return(fdr_perm);
}

//TODO: not very good code - should use more general extraction...
//extract for pvalue from perm
std::vector<double> extractFdrPvalue(const std::vector<double>& pvalue, const std::vector<long long>& perm){
    std::vector<double> fdrPvalue(perm.size());
    unsigned long long correctIndex = 0;
    for (unsigned long long index = 0; index < perm.size(); ++index){
        correctIndex = perm[index];
        //now checking that pvalue has correctIndex-many elements...
        if (correctIndex < pvalue.size()){
            fdrPvalue[index] = pvalue[ correctIndex ];
        }
    }
    return fdrPvalue;
}


//extractFdr for tau from perm
std::vector<long long> extractFdrTau(const std::vector<long long>& tau, const std::vector<long long>& perm){
    std::vector<long long> fdrTau(perm.size());
    unsigned long long correctIndex = 0;
    for (unsigned long long index = 0; index < perm.size(); ++index){
        correctIndex = perm[index];
        //now checking that pvalue has correctIndex-many elements...
        if (correctIndex < tau.size()){
            fdrTau[index] = tau[ correctIndex ];
        }
    }
    return fdrTau;
}


//extract for l from perm
std::vector<long long> extractFdrL(const std::vector<long long>& l, const std::vector<long long>& perm){
    std::vector<long long> fdrL(perm.size());
    unsigned long long correctIndex = 0;
    for (unsigned long long index = 0; index < perm.size(); ++index){
        correctIndex = perm[index];
        //now checking that pvalue has correctIndex-many elements...
        if (correctIndex < l.size()){
            fdrL[index] = l[ correctIndex ];
        }
    }
    return fdrL;
}

