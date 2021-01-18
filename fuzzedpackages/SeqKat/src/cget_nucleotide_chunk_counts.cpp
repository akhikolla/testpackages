#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>

//[[Rcpp::export]]
Rcpp::NumericVector cget_nucleotide_chunk_counts(std::vector<std::string> key, std::string chr, const unsigned int up_stream = 1, 
        const unsigned int down_stream = 1, long start = 1, long end = -1) 
{
        /*
         * DESCRIPTION:
         * ============
         *      Count the occurence of successive nucleotide chunk patterns.
         *
         * INPUT:
         * ======
         *      key: a vector of all possible nucleotide chunk types in UPPER CASE(like trinecleotide), including types with 'N' if any in chr.
         *      chr: chromosome fasta file path. This file should not contain promt '>chrZZ' or newlines as in many cases
         *               One should first use shell command "tr -d '\n' < chr" to remove newlines
         *      up_stream, down_stream: Up and down stream base numbers. Single base if both = 0 and trinucleotdie if both = 1
         *      start, end: start and end position of chromosome. start >= 1, and end = -1 means till the end.
         *
         * OUTPUT:
         * =======
         *      A vector of trinucleotide counts.
         *
         * CALL FROM WITHIN R:
         * ==================
         *      library(Rcpp)
         *      sourceCpp("cget_nucleotide_chunk_counts.cpp")
         *      out <- cget_nucleotide_chunk_counts(key = R_vector_of_nucleotide_pattern, chr = 'chr1.fa')
         *
         * OTHER:
         * =====
         *      Last update: 29-Sep-2015
         *      Author: Eric Xihui Lin <xlin@oicr.on.ca>
        */

        std::map<std::string, long int> count;
        const int length = up_stream + down_stream + 1;

        //initialize
        for (size_t i = 0; i < key.size(); i++) {
                count.insert(std::pair<std::string, long int>(key[i], 0));
                }

        std::ifstream f(chr.c_str());

        char *tri_nucl_c = new char[length + 1]; // c-string (char array)
        std::string tri_nucl; // c++ type string of tri_nucl_c

        if(f.is_open()) 
        {
                //check if the last char is '\n', and decide chr length
                f.seekg(-1, f.end);
                long int chr_length = f.tellg();
                if (f.get() != '\n') ++chr_length;

                if (end <= 0 || end > chr_length) end = chr_length;
                if (start < 1) start = 1;
                start = start -1 ; // start from 0 (c index)


                f.seekg(start, f.beg);
                //if shift = r, then if will f will hit EOF exactly, if not '\n' at the end. So start with shift = r+1
                int r = chr_length % length; 

                int shift;
                for (int i = 0; i < length; i++)
                {
                        shift = (r + 1 + i) % length;
                        f.seekg(start+shift, f.beg);
                        for (long int chunk_end = start+shift+length-1; chunk_end < end; chunk_end += length)
                        {
                                //extract the next `length` charactors, file pointer moves `length` chars forward
                                f.get(tri_nucl_c, length+1); 
                                tri_nucl = tri_nucl_c;
                                transform(tri_nucl.begin(), tri_nucl.end(), tri_nucl.begin(), ::toupper);
                                ++count.at(tri_nucl);
                        }
                }
        } 
        else 
        {
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
