#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//[[Rcpp::export]]
Rcpp::CharacterVector cpp_get_context(std::string chr, Rcpp::NumericVector pos, unsigned long pos_length) {

    std::ifstream f(chr.c_str()); //open file containing chromosome sequence

    char tri_nucleotide[]="NNN";
    Rcpp::CharacterVector context(pos_length, tri_nucleotide);
    
    if(f.is_open()) {
        for (int i=0; i<pos_length; i++) {
            f.seekg(pos[i] - 2, f.beg);
            f.get(tri_nucleotide, 4);
            context[i] = tri_nucleotide;
        }
    } else {
        // Typical CRAN: "Compiled code should not call entry points which might terminate R nor
        // write to stdout/stderr instead of to the console, nor the system RNG.", so we won't printf or exit here.
        // printf("Fail to open file!");
        // exit(EXIT_FAILURE);
    }

    f.close();

    return context;

}
