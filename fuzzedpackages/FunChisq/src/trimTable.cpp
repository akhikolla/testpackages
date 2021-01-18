//
//  trimTable.cpp
//  EFT
//
//  Created by Joe Song on 9/22/19.
//  Copyright © 2019 Joe Song. All rights reserved.
//

#include "trimTable.h"

vector<vector<int>> trimTable(const vector<vector<int>> & table)
{   // DQP::trimTable() removes all-zero rows and all-zero columns
    //   from the input table
    // MS 9/21/2019
    vector<vector<int > > trimmed;

    size_t nrows = table.size();
    size_t ncols = nrows > 0 ? table[0].size() : 0;

    vector<int> rowsums(nrows, 0);
    vector<int> colsums(ncols, 0);

    for(auto i=0u; i<nrows; ++i) {
        for(auto j=0u; j<ncols; ++j) {
            rowsums[i] += table[i][j];
            colsums[j] += table[i][j];
        }
    }

    size_t nrows_trimmed(nrows), ncols_trimmed(ncols);

    for(auto i=0u; i<nrows; ++i) {
        if(rowsums[i] == 0) nrows_trimmed --;
    }

    for(auto j=0u; j<ncols; ++j) {
        if(colsums[j] == 0) ncols_trimmed --;
    }

    if(nrows_trimmed == nrows && ncols_trimmed == ncols) {
        trimmed = table;
    } else {
        trimmed = vector<vector<int > > (nrows_trimmed, vector<int>(ncols_trimmed));
        size_t i_t = 0u;
        for(auto i=0u; i<nrows; ++i) {
            if(rowsums[i] == 0) {
                continue;
            }
            size_t j_t = 0u;
            for(auto j=0u; j<ncols; ++j) {
                if(colsums[j] == 0) {
                    continue;
                }
                trimmed[i_t][j_t] = table[i][j];
                j_t ++;
            }
            i_t ++;
        }
    }

    return trimmed;
}
