/*
 * FilterIntervals.h
 *
 *  Created on: Sep 12, 2016
 *      Author: mabaker
 */

#ifndef FILTERINTERVALS_H_
#define FILTERINTERVALS_H_

//filterIntervals.h
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

#include "types.h"

using std::vector;
using std::string;
using std::stringstream;


const double DEFAULT_SCORE = 0.0;
const double DEFAULT_ODDS_RATIO = 1.0;
const double DEFAULT_PVALUE = 1.0;
const longint DEFAULT_START = 0;
const longint DEFAULT_END = 0;

namespace SignificantPattern
{
    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//
    class Interval{
        public:
            longint getStart() const;
            longint getEnd() const;
            double getScore() const;
            double getOddsRatio() const;
            double getPvalue() const;
            longint getLength() const;
            void setStart(longint);
            void setEnd(longint);
            void setEnd(longint, longint);
            void setScore(double);
            void setOddsRatio(double);
            void setPvalue(double);
            bool overlaps(longint, longint) const;
            void printInterval() const;
            static longint computeEnd(longint tau_, longint l_);
        private:
            longint start;
            longint end;
            double score;
            double odds_ratio;
            double pvalue;
    };
    //--------------------------------------------------------------------------//

    class FilterIntervals
    {
    private:
        vector<Interval> sigInts;
        vector<Interval> getMinPvalueIntervalPerCluster(vector<longint>& tau, vector<longint>& l, vector<double>& score, vector<double>& odds_ratio, vector<double>& pvalue, const vector<int>& label);
        //LP output:
        vector<Interval> sigClusters;
        vector<Interval> getCLusterBoundaries(vector<longint>& tau, vector<longint>& l, const vector<int>& label);


        vector<int> getClusterLabelsForIntervals(const vector<longint>& tau, const vector<longint>& l, const vector<Interval>& cluster);
        vector<Interval> getClusters(vector<longint>& v_tau, vector<longint>& v_l);
        vector<bool> getClusterIndicatorVector(vector<longint>& v_tau, vector<longint>& v_l);
        void makeIntervalTrue(vector<bool>& v, const longint tau, const longint l);
    public:
        FilterIntervals ();
        FilterIntervals (const FilterIntervals& other);
        virtual FilterIntervals& operator=(const FilterIntervals& other);
        virtual ~FilterIntervals ();
        void cpp_filterIntervalsFromMemory(vector<longint> ll_tau,
                                           vector<longint> ll_l,
                                           vector<double> score,
                                           vector<double> odds_ratio,
                                           vector<double> pvalue);
        void writeToFile(const std::string& filename);
        inline vector<Interval>& getSigInts() { return sigInts; }
    };

    class SignificantIntervals
    {
    private:
        vector<Interval> sigInts;

    public:
        SignificantIntervals();
        SignificantIntervals (const SignificantIntervals& other);
        virtual SignificantIntervals& operator=(const SignificantIntervals& other);
        virtual ~SignificantIntervals();
        void cpp_intervalsFromMemory(vector<longint> ll_tau,
                                     vector<longint> ll_l,
                                     vector<double> score,
                                     vector<double> odds_ratio,
                                     vector<double> pvalue);
        void writeToFile(const std::string& filename);
        inline vector<Interval>& getSigInts() { return sigInts; }
    };

} /* namespace SignificantPattern */

#endif /* FILTERINTERVALS_H_ */
