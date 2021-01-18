#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>

//[[Rcpp::export]]
Rcpp::NumericVector cget_trinucleotide_counts(std::vector<std::string> key, std::string chr, const int length = 3, long int start = 1, long int end = -1) {

    std::map<std::string, long int> count;

    //initialize
    for (size_t i = 0; i < key.size(); i++) 
        count.insert( std::pair<std::string, long int>(key[i], 0) );

    std::ifstream f(chr.c_str());

    char *tri_nucl_c = new char[length];
    std::string tri_nucl;

    if(f.is_open()) {
        
        //find the length of genome, chr_length=length_of_chr + 1(end_of_file);
        f.seekg(0, f.end);
        long int chr_length = f.tellg();

        f.seekg(start - 1, f.beg);

        long int back_move_step = length - (length / 2);

        if (end < 0 || end > chr_length) end = chr_length;
        if (start < 1) start = 1;

        for (long int i = start - 1; i <= end - length; i++) {
            f.get(tri_nucl_c, 4); //extract the next 3 charactors, file pointer moves 3 chars
            tri_nucl = tri_nucl_c;
            transform(tri_nucl.begin(), tri_nucl.end(), tri_nucl.begin(), ::toupper);
            f.seekg(- back_move_step, std::ios_base::cur); //move pointer back to 'next of beginning of this loop'
            //printf("%s\n", tri_nucl.c_str());
            ++count.at(tri_nucl);
        }

    } else {
        // Typical CRAN: "Compiled code should not call entry points which might terminate R nor
        // write to stdout/stderr instead of to the console, nor the system RNG.", so we won't printf or exit here.
        // printf("Fail to open file!");
        // exit(EXIT_FAILURE);
    }

    f.close();

    Rcpp::NumericVector count_vec(key.size());

    for (size_t i = 0; i < key.size(); ++i)
        count_vec[i] = count.at(key[i]);

    delete[] tri_nucl_c;
    return count_vec;
}
