/*
 * FilterIntervals.cpp
 *
 *  Created on: Sep 12, 2016
 *
 */

#include "FilterIntervals.h"

//filterIntervals.cpp
//      read csv file containing all significant intervals, group these into
//      overlapping clusters, and return the most significant interval per cluster
//
//Laetitia Papaxanthos, Dean Bodenham
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<iomanip>
#include<numeric>
#include<iterator>
#include<algorithm> //sort
#include<cmath>
//#include<Rcpp.h>

#include "Exception.h"

using namespace std;
std::ostringstream oss;


namespace SignificantPattern
{
    // Constructor
    FilterIntervals::FilterIntervals ()
    {
    }
    // Copy constructor
    FilterIntervals::FilterIntervals (const FilterIntervals& other)
    {
        // call copy assignment operator
        *this = other;
    }
    // Copy assignment operator
    FilterIntervals& FilterIntervals::operator=(const FilterIntervals& other)
    {
       if (this != &other) {
           sigInts.clear(); // should be unnecessary
           sigInts = other.sigInts;
           sigClusters = other.sigClusters;
       }
       return *this;
    }
    // Destructor
    FilterIntervals::~FilterIntervals ()
    {
    }

    longint Interval::computeEnd(longint tau_, longint l_) {
        return (tau_ + l_ - 1);
    }

    longint Interval::getStart() const { return start;}
    longint Interval::getEnd() const { return end;}
    double Interval::getScore() const {return score;}
    double Interval::getOddsRatio() const {return odds_ratio;}
    double Interval::getPvalue() const {return pvalue;}
    longint Interval::getLength() const { return end - start + 1;}
    void Interval::setStart(longint m_start){ start = m_start;}
    void Interval::setEnd(longint m_end){ end = m_end;}
    void Interval::setEnd(longint tau_, longint l_){ end = computeEnd(tau_, l_);}
    void Interval::setScore(double m_score){ score = m_score;}
    void Interval::setOddsRatio(double m_odds_ratio){ odds_ratio = m_odds_ratio;}
    void Interval::setPvalue(double m_pvalue){ pvalue = m_pvalue;}

    //check if current interval overlaps interval [a,b]
    bool Interval::overlaps(longint a, longint b) const{
        bool doesOverlap = false;
        //very simple check: if start <= a <= end, then overlap is true
        if  (  (start <= a)  && (a <= end ))
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
    longint getMaxIntervalEndpoint(vector<longint> v_tau, vector<longint> v_l){
        longint maxEndpoint = 0;
        //longint start;
        longint end;
        for (vector<longint>::size_type index = 0; index != v_tau.size(); ++index){
            //start = v_tau[index];
            //end = start + v_l[index] - 1;
            end = Interval::computeEnd(v_tau[index], v_l[index]);
            maxEndpoint = (maxEndpoint < end) ? end : maxEndpoint;
        }
        return maxEndpoint;
    }


    //change interval [tau, tau+1, ..., tau+l - 1] to true
    void FilterIntervals::makeIntervalTrue(vector<bool>& v, const longint tau, const longint l){
        //this was a problem, not checking that index < v.size()
        longint count = 0;
        const longint maxIndexTrue = tau + l;

        //as soon as iterator reaches end OR enough values set to true, will stop loop. Need to use non-constant iterator here, because we change the value in the vector
        const longint interval_start = tau;
        const longint interval_end = Interval::computeEnd((longint)tau, (longint)l);

        for(vector<bool>::iterator it = v.begin() + interval_start - 1;
                    it != v.begin() + interval_end  && count < maxIndexTrue; ++it){
            //remember - C++ counts from 0, but we have corrected for this when
            //constructing the vector
            *it = true;
            count++;
        }
    }


    //get vector which indicates the cluster
    vector<bool> FilterIntervals::getClusterIndicatorVector(vector<longint>& v_tau, vector<longint>& v_l){




        //get maxEndpoint, so can create vector of false's of that size
        longint n = getMaxIntervalEndpoint(v_tau, v_l);
        //not going to mess around with 0th index... will keep "true" indices
        n++;

        //instead of construction, use push_back
        vector<bool> clusterIndicator(n, false);

        //create two iterators, and move over both iterators
        vector<longint>::const_iterator it_tau = v_tau.begin();
        vector<longint>::const_iterator it_l = v_l.begin();

        //just get rid of zeroth index
        //This was an error
        for (; it_tau != v_tau.end() and it_l != v_l.end(); ++it_tau, ++it_l){
            makeIntervalTrue(clusterIndicator, *it_tau, *it_l);
        }

        return clusterIndicator;

    }


    //run through binary array, and label start/end of clusters
    vector<Interval> FilterIntervals::getClusters(vector<longint>& v_tau, vector<longint>& v_l){

        vector<Interval> cluster;
        vector<bool> clusterIndicator = getClusterIndicatorVector(v_tau, v_l);
        Interval thisCluster;
        bool inCluster = false;
        longint index = 0;
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
    vector<int> FilterIntervals::getClusterLabelsForIntervals(const vector<longint>& tau, const vector<longint>& l, const vector<Interval>& cluster){
        //get tau and l
        //the label vector will be returned
        vector<int> label(tau.size());

        //these ints are the temp label variables
        int thisIntervalLabel;
        longint thisIntervalEnd;

        //create a vector of the same size as the clusters
        vector<int> clusterLabels(cluster.size());
        //fill it with 0, 1, ..., cluster.size()
        for (unsigned int i=0; i < clusterLabels.size(); ++i){
            clusterLabels[i] = i;
        }

        //run through all intervals/tau's, iterating over both tau and l
        //although l is not needed...
        vector<longint>::const_iterator it_tau = tau.begin();
        vector<longint>::const_iterator it_l = l.begin();

        for(longint intervalIndex = 0; it_tau != tau.end() && it_l != l.end(); ++it_tau, ++it_l, ++intervalIndex){
            //in case no cluster assigned - which should not happen - interval will get label "0"
            //setting default value to 0, although could be -1...but worried about any error would kill the process
    //         thisIntervalLabel = -1;
            thisIntervalLabel = 0;
            thisIntervalEnd = Interval::computeEnd((longint)*it_tau, (longint)*it_l);

            //now check which cluster it belongs to:
            //this for-loop is iterating over clusters, which are intervals.
            //we are checking for a particular tau (original interval) if it overlaps with THIS cluster
            //if so, we give it that cluster label
            vector<int>::const_iterator it_thisClusterLabel = clusterLabels.begin();

            for (vector<Interval>::const_iterator it_cluster = cluster.begin();
                it_cluster != cluster.end() && it_thisClusterLabel != clusterLabels.end();
                ++it_cluster, ++it_thisClusterLabel){
                if (  (*it_cluster).overlaps((longint)*it_tau, (longint)thisIntervalEnd)  ){
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
    vector<Interval> FilterIntervals::getMinPvalueIntervalPerCluster(vector<longint>& tau, vector<longint>& l, vector<double>& score, vector<double>& odds_ratio, vector<double>& pvalue, const vector<int>& label){

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
            (*it_sigInts).setScore(DEFAULT_SCORE);
            (*it_sigInts).setOddsRatio(DEFAULT_ODDS_RATIO);
            (*it_sigInts).setPvalue(DEFAULT_PVALUE);
        }

        //iterate over tau, l, pvalue and label all at once
        vector<double>::const_iterator it_score = score.begin();
        vector<double>::const_iterator it_odds_ratio = odds_ratio.begin();
        vector<double>::const_iterator it_pvalue = pvalue.begin();
        vector<longint>::const_iterator it_tau = tau.begin();
        vector<longint>::const_iterator it_l = l.begin();
        vector<int>::const_iterator it_label = label.begin();
        double currentMinPvalueForCluster;
        longint currentMinPvalueForClusterTau;
        longint currentMinPvalueForClusterL;

        for (; it_score != score.end() && it_odds_ratio != odds_ratio.end() && it_pvalue != pvalue.end() && it_tau != tau.end() && it_l != l.end() && it_label != label.end();
               ++it_score, ++it_odds_ratio, ++it_pvalue, ++it_tau, ++it_l, ++it_label){

            //see what the current min pvalue is for this cluster
            currentMinPvalueForCluster = sigInts[*it_label].getPvalue();
            currentMinPvalueForClusterTau = sigInts[*it_label].getStart();
            //need to compute l
            currentMinPvalueForClusterL = sigInts[*it_label].getLength();
            if (*it_pvalue < currentMinPvalueForCluster){
                //if pvalue for this interval is smaller than the min pvalue for this cluster, update the interval
                sigInts[*it_label].setStart(*it_tau);
                sigInts[*it_label].setEnd(*it_tau, *it_l);
                sigInts[*it_label].setScore(*it_score);
                sigInts[*it_label].setOddsRatio(*it_odds_ratio);
                sigInts[*it_label].setPvalue(*it_pvalue);
            }

            //what if p-values are equal?
            //If the p-values are equal, we keep the SMALLEST interval.
            //It follows the objective of the clustering which is to remove the redundancy.
            //It is consistent with the FastCMH publication.
            //Then we keep interval with smaller tau (if interval lengths the SAME)
            //This decision is made in order to ensure consistent results.

            //An alternative is to keep the LARGEST interval.
            //To do so, change the sign in the line commented by " sign < for SMALLEST interval".

            if (*it_pvalue == currentMinPvalueForCluster){
                //First decision: if interval length is the SAME
                //              AND tau is smaller - REPLACE
                if (*it_l == currentMinPvalueForClusterL){
                    if (*it_tau < currentMinPvalueForClusterTau){
                        sigInts[*it_label].setStart(*it_tau);
                        sigInts[*it_label].setEnd(*it_tau, *it_l);
                        sigInts[*it_label].setScore(*it_score);
                        sigInts[*it_label].setOddsRatio(*it_odds_ratio);
                        sigInts[*it_label].setPvalue(*it_pvalue);
                    }
                } else {
                    //Second decision: if interval length is larger, REPLACE
                    if (*it_l < currentMinPvalueForClusterL){ //sign < for SMALLEST interval
                        sigInts[*it_label].setStart(*it_tau);
                        sigInts[*it_label].setEnd(*it_tau, *it_l);
                        sigInts[*it_label].setScore(*it_score);
                        sigInts[*it_label].setOddsRatio(*it_odds_ratio);
                        sigInts[*it_label].setPvalue(*it_pvalue);
                    }
                }
            } //end of additional case when p-values are equal

        }//end of for



        return sigInts;
    }

    //--------------------------------------------------------------------------//
    //Get number of regions per cluster and leftmost and rightmost SNPs
    vector<Interval> FilterIntervals::getCLusterBoundaries(vector<longint>& tau, vector<longint>& l, const vector<int>& label){

        //need to quickly get the max and min cluster values
        int maxCluster = 0;
        int minCluster = 2;
        for (vector<int>::const_iterator it = label.begin(); it != label.end(); ++it){
            if (*it > maxCluster)
                maxCluster = *it;
            if (*it < minCluster)
                minCluster = *it;
        }
        int numClusters = maxCluster - minCluster +1;

        vector<Interval> sigClusters(numClusters);
        for (vector<Interval>::iterator it_sigClusters = sigClusters.begin(); it_sigClusters != sigClusters.end(); ++it_sigClusters){
            (*it_sigClusters).setStart(-1);
            (*it_sigClusters).setEnd(-1);
            (*it_sigClusters).setPvalue(0);
        }


        vector<int>::const_iterator it_label = label.begin();
        vector<longint>::const_iterator it_tau = tau.begin();
        vector<longint>::const_iterator it_l = l.begin();
        int currentStartCluster;
        int currentEndCluster;
        int currentNumCluster;
    for (; it_label != label.end() ; ++it_tau, ++it_l, ++it_label){
    	currentStartCluster = sigClusters[*it_label].getStart();
    	currentEndCluster = sigClusters[*it_label].getEnd();
    	currentNumCluster = sigClusters[*it_label].getPvalue();
    	if (currentStartCluster > *it_tau or currentStartCluster == -1){
    		sigClusters[*it_label].setStart(*it_tau);
    	}
    	if (currentEndCluster < *it_tau + *it_l - 1){
    		sigClusters[*it_label].setEnd(*it_tau, *it_l);
    	}
    	sigClusters[*it_label].setPvalue(currentNumCluster+1);
    }

    	return sigClusters;
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
    void FilterIntervals::cpp_filterIntervalsFromMemory(vector<longint> ll_tau,
                                                        vector<longint> ll_l,
                                                        vector<double> score,
                                                        vector<double> odds_ratio,
                                                        vector<double> pvalue){

        if (pvalue.size() > 0){
            //create data frame
            vector<longint> tau(ll_tau.begin(), ll_tau.end());
            vector<longint> l(ll_l.begin(), ll_l.end());

            //get vector of clusters (cluster labels)
            vector<Interval> cluster = getClusters(tau, l);

            //now get labels for cluster vectors
            vector<int> label = getClusterLabelsForIntervals(tau, l, cluster);

            //get minimum p-values
            sigInts = getMinPvalueIntervalPerCluster(tau, l, score, odds_ratio, pvalue, label);

            //get the boundaries
            sigClusters = getCLusterBoundaries(tau, l, label);

            //sort the significant intervals
            sortIntervals(sigInts);
        }
        else
            sigInts = createEmptyInterval();
    }

    void FilterIntervals::writeToFile(const std::string& filename)
    {
        ofstream file;
        file.exceptions ( ifstream::failbit | ifstream::badbit );
        try
        {
                file.open(filename.c_str());
        }
        catch (const std::ios_base::failure& e)
        {
                throw Exception("Failed opening " + filename + " for writing");
        }

        file << "P-value" << "\t" << "score" << "\t" << "OR" << "\t" << "index_1;...;index_N" << "\t" <<"#n regions" << "\t" <<"index L SNP" << "\t" <<"index R SNP" << endl;
		for (vector<Interval>::size_type i=0; i<sigInts.size(); ++i){
			//Compute 'index_1;index_2;index_3;...;index_N'
			vector<int> myVec;
			std::ostringstream oss;

			for(int j=sigInts[i].getStart(); j<= sigInts[i].getEnd(); j++){
				myVec.push_back(j);
			}
			if (!myVec.empty())
			  {
				// Convert all but the last element to avoid a trailing ";"
				std::copy(myVec.begin(), myVec.end()-1, std::ostream_iterator<int>(oss, ";"));

				// Now add the last element with no delimiter
				oss << myVec.back();
			  }
			file << sigInts[i].getPvalue() << "\t" << sigInts[i].getScore() << "\t" << sigInts[i].getOddsRatio() << "\t" << oss.str() << "\t" << sigClusters[i].getPvalue() << "\t" << sigClusters[i].getStart() << "\t" << sigClusters[i].getEnd() << endl;
		}
		file.close();



    }

    // Constructor
    SignificantIntervals::SignificantIntervals ()
    {
    }
    // Copy constructor
    SignificantIntervals::SignificantIntervals (const SignificantIntervals& other)
    {
        // call copy assignment operator
        *this = other;
    }
    // Copy assignment operator
    SignificantIntervals& SignificantIntervals::operator=(const SignificantIntervals& other)
    {
       if (this != &other) {
           sigInts.clear(); // should be unnecessary
           sigInts = other.sigInts;
       }
       return *this;
    }
    // Destructor
    SignificantIntervals::~SignificantIntervals ()
    {
    }

    void SignificantIntervals::cpp_intervalsFromMemory(vector<longint> ll_tau,
                                                       vector<longint> ll_l,
                                                       vector<double> score,
                                                       vector<double> odds_ratio,
                                                       vector<double> pvalue){

        if (pvalue.size() > 0){
            //create data frame
            vector<longint> tau(ll_tau.begin(), ll_tau.end());
            vector<longint> l(ll_l.begin(), ll_l.end());
            vector<longint>::size_type n = ll_tau.size();
            if (n != ll_l.size() || n != pvalue.size())
                throw Exception("given vectors ll_tau, ll_l and pvalue don't have same size");

            sigInts = vector<Interval>(n);
            for (vector<longint>::size_type i=0; i<n; ++i)
            {
                sigInts.at(i).setStart(ll_tau.at(i));
                sigInts.at(i).setEnd(ll_tau.at(i), ll_l.at(i));
                sigInts.at(i).setScore(score.at(i));
                sigInts.at(i).setOddsRatio(odds_ratio.at(i));
                sigInts.at(i).setPvalue(pvalue.at(i));
            }
        }
        else
            sigInts = createEmptyInterval();
    }

    void SignificantIntervals::writeToFile(const std::string& filename)
    {
        ofstream file;
        file.exceptions ( ifstream::failbit | ifstream::badbit );
        try
        {
                file.open(filename.c_str());
        }
        catch (const std::ios_base::failure& e)
        {
                throw Exception("Failed opening " + filename + " for writing");
        }

        file << "P-value" << "\t" << "score" << "\t" << "OR" << "\t" << "index_1;index_2;index_3;...;index_N" << endl;
        for (vector<Interval>::size_type i=0; i<sigInts.size(); ++i){
            file << sigInts[i].getStart() << " " << sigInts[i].getEnd() << " " << sigInts[i].getPvalue() << " " << sigInts[i].getLength() << endl;}
        file.close();}


} /* namespace SignificantPattern */
