#ifndef __TREE_H__
#define __TREE_H__

#include <vector>

using std::vector;

class Tree
{
public:
    int num_tips;                   
    double total_time;
    vector<int> speciators;         // Order of splitting lines
    vector<double> intervals;       // Time periods between speciation events
};

#endif
