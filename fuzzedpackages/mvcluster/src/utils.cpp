/*
 * utils.cpp
 *
 *  Created on: Aug 27, 2015
 *      Author: Javon
 */

#include <vector>
#include <algorithm>
#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

uint32_t countNoneZero(const vec& v) {
	uint32_t cnt = 0;
	for (uint32_t i = 0; i < v.size(); i++) {
		if (v[i] != 0) {
			cnt++;
		}
	}
	return cnt;
}

double findMedian(vector<double> vec){
//    Find median of a vector
    double median;
    size_t size = vec.size();
    median = vec[(size/2)];
    return median;
}

double findMedianOfMedians(vector<vector<double> > values){
    vector<double> medians;

    for (uint32_t i = 0; i < values.size(); i++) {
        double m = findMedian(values[i]);
        medians.push_back(m);
    }

    return findMedian(medians);
}

double selectionByMedianOfMedians(const vector<double>& values, uint32_t k) {
//    Divide the list into n/5 lists of 5 elements each
    vector<vector<double> > vec2D;

    uint32_t count = 0;
    while (count != values.size()) {
    	uint32_t countRow = 0;
        vector<double> row;

        while ((countRow < 5) && (count < values.size())) {
            row.push_back(values[count]);
            count++;
            countRow++;
        }
        vec2D.push_back(row);
    }

//    cout<<endl<<endl<<"Printing 2D vector : "<<endl;
//    for (int i = 0; i < vec2D.size(); i++) {
//        for (int j = 0; j < vec2D[i].size(); j++) {
//            cout<<vec2D[i][j]<<" ";
//        }
//        cout<<endl;
//    }
//    cout<<endl;

//    Calculating a new pivot for making splits
    double m = findMedianOfMedians(vec2D);
//    cout<<"Median of medians is : "<<m<<endl;

//    Partition the list into unique elements larger than 'm' (call this sublist L1) and
//    those smaller them 'm' (call this sublist L2)
    vector<double> L1, L2;
    uint32_t nm = 0; // amount of number equal to m
    for (uint32_t i = 0; i < vec2D.size(); i++) {
        for (uint32_t j = 0; j < vec2D[i].size(); j++) {
            if (vec2D[i][j] > m) {
                L1.push_back(vec2D[i][j]);
            }
            else if (vec2D[i][j] < m){
                L2.push_back(vec2D[i][j]);
            }
            else {
            	nm++;
            }
        }
    }

//    Checking the splits as per the new pivot 'm'
//    cout<<endl<<"Printing L1 : "<<endl;
//    for (uint i = 0; i < L1.size(); i++) {
//        cout<<L1[i]<<" ";
//    }
//
//    cout<<endl<<endl<<"Printing L2 : "<<endl;
//    for (uint i = 0; i < L2.size(); i++) {
//        cout<<L2[i]<<" ";
//    }

//    Recursive calls
    if (k < L1.size() + 1) {
		return selectionByMedianOfMedians(L1, k);
	}
	else if (k >= L1.size() + 1 && k <= L1.size() + nm) {
		return m;
	}
	else {
		return selectionByMedianOfMedians(L2, k-((uint32_t)L1.size())-nm);
	}
}



