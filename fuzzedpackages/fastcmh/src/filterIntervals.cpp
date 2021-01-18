//filterIntervals.cpp
//      read csv file containing all significant intervals, group these into
//      overlapping clusters, and return the most significant interval per cluster
//
//Dean Bodenham March 2016
#include<iostream>
#include<string>
#include<vector>
#include<iomanip>
#include<numeric>
#include<algorithm>
#include<Rcpp.h>
#include "filterIntervals.h"

//add extra sort step, so result of filtering is the same, no matter the initial order of the intervals?

//--------------------------------------------------------------------------//

size_t computeEnd(size_t tau_, size_t l_){
    return (tau_ + l_ - 1);
}


size_t Interval::getStart() const { return start;}
size_t Interval::getEnd() const { return end;}
double Interval::getPvalue() const {return pvalue;}
size_t Interval::getLength() const { return end - start + 1;}
void Interval::setStart(size_t m_start){ start = m_start;}
void Interval::setEnd(size_t m_end){ end = m_end;}
void Interval::setEnd(size_t tau_, size_t l_){ end = computeEnd(tau_, l_);}
void Interval::setPvalue(double m_pvalue){ pvalue = m_pvalue;}

//check if current interval overlaps interval [a,b]
bool Interval::overlaps(size_t a, size_t b) const{
    bool doesOverlap = false;
    //very simple check: if start <= a <= end, then overlap is true 
    if (  (start <= a)  && (a <= end)  )
        doesOverlap = true;
    return doesOverlap;

}

//--------------------------------------------------------------------------//

//for sorting, need to create an operator for comparing vectors
struct less_than_Interval{
    inline bool operator() (const Interval& interval1, const Interval& interval2){
        return (interval1.getStart() < interval2.getStart());
    }
};

//sort a vector of intervals
void sortIntervals(vector<Interval>& v){
    std::sort(v.begin(), v.end(), less_than_Interval());
}
//--------------------------------------------------------------------------//
//--------------------------------------------------------------------------//
//algo:
//done NOT needed (so not done) 1) run through all intervals, creating "start" "end" and "pvalue" vectors
//done 2) find max point (so we know how long to make the binary array
//done 3) run through binary array, and label start/end of clusters
//4) run through intervals, and create label vector which assigns to cluster
//5) count clusters, and create vector that will save min-pvalue-interval per cluster
//6) run through (now labelled) intervals, for cluster 1, 2, ..., and for each cluster, find the interval with the smallest pvalue. Save the "start", "end", "pvalue"


//a fundamental method:
//get the largest index of the interval
//e.g. a=5, tau=3 -> interval [5,7] -> endpoint=7
//BUT also need to makesure this does not exceed max number of features
size_t getMaxIntervalEndpoint(vector<size_t> v_tau, vector<size_t> v_l){
    size_t maxEndpoint = 0;
//     size_t start, end;
    size_t end;
    for (size_t index = 0; index != v_tau.size(); ++index){
        end = computeEnd(v_tau[index], v_l[index]);
        maxEndpoint = (maxEndpoint < end) ? end : maxEndpoint;
    }
    return maxEndpoint;
}


//change interval [tau, tau+1, ..., tau+l - 1] to true
void makeIntervalTrue(vector<bool>& v, const size_t tau, const size_t l){
    //this was a problem, not checking that index < v.size()
    size_t count = 0;
    const size_t maxIndexTrue = tau + l;

    //as soon as iterator reaches end OR enough values set to true, will stop loop. Need to use non-constant iterator here, because we change the value in the vector
    const size_t interval_start = tau;
    const size_t interval_end = computeEnd(tau, l);

    for(vector<bool>::iterator it = v.begin() + interval_start; 
                it != v.begin() + interval_end && count < maxIndexTrue; ++it){
        //remember - C++ counts from 0, but we have corrected for this when
        //constructing the vector
        *it = true;
        count++;
    }
}


//get vector which indicates the cluster
vector<bool> getClusterIndicatorVector(vector<size_t>& v_tau, vector<size_t>& v_l){
    //get maxEndpoint, so can create vector of false's of that size
    size_t n = getMaxIntervalEndpoint(v_tau, v_l);
    //not going to mess around with 0th index... will keep "true" indices
    n++;

    //instead of construction, use push_back
    vector<bool> clusterIndicator(n, false);

    //create two iterators, and move over both iterators
    vector<size_t>::const_iterator it_tau = v_tau.begin();
    vector<size_t>::const_iterator it_l = v_l.begin();

    //just get rid of zeroth index
    //This was an error
    for (; it_tau != v_tau.end() and it_l != v_l.end(); ++it_tau, ++it_l){
        makeIntervalTrue(clusterIndicator, *it_tau, *it_l);
    }
    return clusterIndicator;

}


//run through binary array, and label start/end of clusters
vector<Interval> getClusters(vector<size_t>& v_tau, vector<size_t>& v_l){

    vector<Interval> cluster;
    vector<bool> clusterIndicator = getClusterIndicatorVector(v_tau, v_l);
    Interval thisCluster;
    bool inCluster = false;
    size_t index = 0;
    for (vector<bool>::const_iterator it = clusterIndicator.begin(); it != clusterIndicator.end(); ++it, ++index){

        //starts off false, as soon as switches to true...
        if (*it){
            //if not in cluster, this must be the start of a cluster
            if (!inCluster){
                //switch cluster state
                inCluster = true;
                //save start of cluster
                thisCluster.setStart(index);
                //not necessary, but worth setting pvalue
                thisCluster.setPvalue(DEFAULT_PVALUE);
            }
        } else{
            //okay, we've seen a false, but are in a cluster...must be the end
            if (inCluster){
                //save the end of the cluster
                thisCluster.setEnd(index);
                //and push_back the cluster to the vector
                cluster.push_back(thisCluster);
                //switch state back to false
                inCluster = false;
            }
        }
    

    }//end of for loop
    //if loop ends, but still in cluster, must be end of cluster
    if (inCluster){
            //need to subtract from index...
            index--;
            thisCluster.setEnd(index);
            //and push_back the cluster to the vector
            cluster.push_back(thisCluster);
            //switch state back to false
            inCluster = false;
    }

    return cluster;
}

//--------------------------------------------------------------------------//
//4) run through intervals, and create label vector which assigns to cluster
//labels will start counting from 0

//for(vector<bool>::iterator it = v.begin() + tau; it != v.end() && count < maxIndexTrue; ++it){
//we use consts
vector<int> getClusterLabelsForIntervals(const vector<size_t>& tau, const vector<size_t>& l, const vector<Interval>& cluster){
    //get tau and l
    //the label vector will be returned
    vector<int> label(tau.size());

    //these ints are the temp label variables
    int thisIntervalLabel;
    size_t thisIntervalEnd;

    //create a vector of size as the same as clusters
    vector<int> clusterLabels(cluster.size());
    //fill it with 0, 1, ..., cluster.size()
    for (unsigned int i=0; i < clusterLabels.size(); ++i){
        clusterLabels[i] = i;
    }
    
    //run through all intervals/tau's, iterating over both tau and l
    //although l is not needed...
    vector<size_t>::const_iterator it_tau = tau.begin();
    vector<size_t>::const_iterator it_l = l.begin();
    
    for(size_t intervalIndex = 0; it_tau != tau.end() && it_l != l.end(); ++it_tau, ++it_l, ++intervalIndex){
        //in case no cluster assigned - which should not happen - interval will get label "0"
        //setting default value to 0, although could be -1...but worried about any error would kill the process
//         thisIntervalLabel = -1;
        thisIntervalLabel = 0;
        thisIntervalEnd = computeEnd(*it_tau, *it_l);
         
        //now check which cluster it belongs to:
        //this for-loop is iterating over clusters, which are intervals.
        //we are checking for a particular tau (original interval) if it overlaps with THIS cluster
        //if so, we give it that cluster label
        vector<int>::const_iterator it_thisClusterLabel = clusterLabels.begin();

        for (vector<Interval>::const_iterator it_cluster = cluster.begin(); 
            it_cluster != cluster.end() && it_thisClusterLabel != clusterLabels.end(); 
            ++it_cluster, ++it_thisClusterLabel){
            if (  (*it_cluster).overlaps(*it_tau, thisIntervalEnd)  ){
                //this interval label is this cluster
                thisIntervalLabel = *it_thisClusterLabel;
                //push_back
                label[intervalIndex] = thisIntervalLabel;
            }
        
        }//for it_cluster, over clusters

    }//end for it_tau, over intervals
    return label;
}



//--------------------------------------------------------------------------//
//run through all clusters (in label) and through all interval...
//i) find number of clusters
//ii) run through label, for each cluster, 
//iii) when label==thisLabel, check p.value, and save interval if smaller pvalue 
vector<Interval> getMinPvalueIntervalPerCluster(vector<size_t>& tau, vector<size_t>& l, vector<double>& pvalue, const vector<int>& label){
    
    //need to quickly get the max and min cluster values
    //use 0 and 2 as max and min, because these should definitely change
    //probably can use 1 instead of 2?
    int maxCluster = 0;
    int minCluster = 2;
    for (vector<int>::const_iterator it = label.begin(); it != label.end(); ++it){
        if (*it > maxCluster)
            maxCluster = *it;
        if (*it < minCluster)
            minCluster = *it;
    }
    int numClusters = maxCluster - minCluster +1;

    //create a vector of Intervals, which will have min pvalues...
    vector<Interval> sigInts(numClusters);
    //set all pvalues to 1, and starting/end points to 0, for initialisation
    //TODO: should I just have a constructor/method?
    for (vector<Interval>::iterator it_sigInts = sigInts.begin(); it_sigInts != sigInts.end(); ++it_sigInts){
        (*it_sigInts).setStart(DEFAULT_START);
        (*it_sigInts).setEnd(DEFAULT_END);
        (*it_sigInts).setPvalue(DEFAULT_PVALUE);
    }
    
    //iterate over tau, l, pvalue and label all at once
    vector<double>::const_iterator it_pvalue = pvalue.begin();
    vector<size_t>::const_iterator it_tau = tau.begin();
    vector<size_t>::const_iterator it_l = l.begin();
    vector<int>::const_iterator it_label = label.begin();
    double currentMinPvalueForCluster;
    size_t currentMinPvalueForClusterTau;
    size_t currentMinPvalueForClusterL;
    for (; it_pvalue != pvalue.end() && it_tau != tau.end() && it_l != l.end() && it_label != label.end();
            ++it_pvalue, ++it_tau, ++it_l, ++it_label){

        //see what the current min pvalue is for this cluster
        currentMinPvalueForCluster = sigInts[*it_label].getPvalue();
        currentMinPvalueForClusterTau = sigInts[*it_label].getStart();
        //need to compute l
        currentMinPvalueForClusterL = sigInts[*it_label].getLength();
        if (*it_pvalue < currentMinPvalueForCluster){
            //if pvalue for this interval is smaller than the min pvalue for this cluster, update the interval
            sigInts[*it_label].setStart(*it_tau);
            sigInts[*it_label].setEnd(*it_tau, *it_l);
            sigInts[*it_label].setPvalue(*it_pvalue);
        }

        //ADDED ADDITIONAL CASE: what if p-values are equal?
        //Then we keep interval with smaller tau (if interval lengths the SAME)
        //This decision is made in order to ensure consistent results
        //If the p-values are equal, we keep the LARGEST interval
        //This is motivated by an example using sample data.
        //The true  interval is [100, 103], but 
        //[100, 100], [100, 101], [101, 102], [100, 103], etc
        //all have the same p-value
        //So, we keep the interval with the LARGEST length
        //
        //It is possible that [100, 101], for example, will have a smaller
        //p-value that [100, 103], the "true" interval, but then [100, 101]
        //will be returned, because it has a smaller p-value
        //This argument is just used when p-values are EQUAL
        if (*it_pvalue == currentMinPvalueForCluster){
            //First decision: if interval length is the SAME
            //              AND tau is smaller - REPLACE
            if (*it_l == currentMinPvalueForClusterL){
                if (*it_tau < currentMinPvalueForClusterTau){
                    sigInts[*it_label].setStart(*it_tau);
                    sigInts[*it_label].setEnd(*it_tau, *it_l);
                    sigInts[*it_label].setPvalue(*it_pvalue);
                }
            } else {
                //Second decision: if interval length is larger, REPLACE
                if (*it_l > currentMinPvalueForClusterL){
                    sigInts[*it_label].setStart(*it_tau);
                    sigInts[*it_label].setEnd(*it_tau, *it_l);
                    sigInts[*it_label].setPvalue(*it_pvalue);
                }
            }
        } //end of additional case when p-values are equal

    }//end of for

    return sigInts;
}


//--------------------------------------------------------------------------//
//create an empty Interval object
vector<Interval> createEmptyInterval(){
    vector<Interval> emptyInterval;
    Interval thisInterval;
    thisInterval.setStart(DEFAULT_START);
    thisInterval.setEnd(DEFAULT_END);
    thisInterval.setPvalue(DEFAULT_PVALUE);

    return(emptyInterval);
} 







//--------------------------------------------------------------------------//
//MAIN FUNCTION
vector<Interval> cpp_filterIntervalsFromMemory(vector<long long> ll_tau,
                                              vector<long long> ll_l,
                                              vector<double> pvalue){

    if (pvalue.size() > 0){
        //create data frame
        //need to case long long to size_t; not necessary for double, but do it anyway
        vector<size_t> tau(ll_tau.begin(), ll_tau.end());
        vector<size_t> l(ll_l.begin(), ll_l.end());

        //get vector of clusters (cluster labels)
        vector<Interval> cluster = getClusters(tau, l);

        //now get labels for cluster vectors
        vector<int> label = getClusterLabelsForIntervals(tau, l, cluster);

        //get minimum p-values
        vector<Interval> sigInts = getMinPvalueIntervalPerCluster(tau, l, pvalue, label);

        //sort the significant intervals
        sortIntervals(sigInts);

        //create data frame
        return(sigInts);
    } 
        
    vector<Interval> sigInts = createEmptyInterval();
    return(sigInts);
}


