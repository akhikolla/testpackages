#include <Rcpp.h>
#include <iostream>
#include <map>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


//' Test DNA sequences for equality.
//' 
//' \code{seqEqual} checks if two DNA sequences are identical.
//'
//' @param    seq1    character string containing a DNA sequence.
//' @param    seq2    character string containing a DNA sequence.
//' @param    ignore  vector of characters to ignore when testing for equality.
//'                   Default is to ignore c("N",".","-","?")
//' 
//' @return   Returns \code{TRUE} if sequences are equal and \code{FALSE} if they are not.
//'           Sequences of unequal length will always return \code{FALSE} regardless of
//'           their character values.
//' 
//' @seealso  Used by \link{pairwiseEqual} within \link{collapseDuplicates}.
//'           See \link{seqDist} for calculation Hamming distances between sequences.
//' 
//' @examples
//' # Ignore gaps
//' seqEqual("ATG-C", "AT--C")
//' seqEqual("ATGGC", "ATGGN")
//' seqEqual("AT--T", "ATGGC")
//' 
//' # Ignore only Ns
//' seqEqual("ATG-C", "AT--C", ignore="N")
//' seqEqual("ATGGC", "ATGGN", ignore="N")
//' seqEqual("AT--T", "ATGGC", ignore="N")
//' 
//' @export
// [[Rcpp::export]]
bool seqEqual(std::string seq1, std::string seq2, 
              CharacterVector ignore=CharacterVector::create("N","-",".","?")) {
    
    int ig_len = ignore.length();
    
    ig_len = ignore.length();
    
    int len_seq1 = seq1.length();
    int len_seq2 = seq2.length();
    
    if (len_seq1 != len_seq2) { 
        return (FALSE);
    } else {
        for(int i = 0; i < len_seq1; i++)
        {
            char seq1_char = (char)seq1[i];
            char seq2_char = (char)seq2[i];
            
            if (seq1_char != seq2_char) {
                
                bool ignore_seq1 = FALSE;
                bool ignore_seq2 = FALSE;
                
                for(int j = 0; j < ig_len; j++) {
                    
                    char ig = *(char*)ignore[j];
                    
                    if (ig == seq1_char) {
                        ignore_seq1 = TRUE;
                    }
                    
                    if (ig == seq2_char) {
                        ignore_seq2 = TRUE;
                    }
                }
                if (!ignore_seq1 & !ignore_seq2) {
                    return FALSE;
                }
            }
        }
        return TRUE;
    }  
}


//' Calculate pairwise equivalence between sequences
//' 
//' \code{pairwiseEqual} determined pairwise equivalence between a pairs in a 
//' set of sequences, excluding ambiguous positions (Ns and gaps).
//'
//' @param    seq  character vector containing a DNA sequences.
//'
//' @return   A logical matrix of equivalence between each entry in \code{seq}. 
//'           Values are \code{TRUE} when sequences are equivalent and \code{FALSE}
//'           when they are not.
//' 
//' @seealso  Uses \link{seqEqual} for testing equivalence between pairs.
//'           See \link{pairwiseDist} for generating a sequence distance matrix.
//'           
//' @examples
//' # Gaps and Ns will match any character
//' seq <- c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C", E="NTGGG")
//' d <- pairwiseEqual(seq)
//' rownames(d) <- colnames(d) <- seq
//' d
//' 
//' @export
// [[Rcpp::export]]
LogicalMatrix pairwiseEqual(StringVector seq) {
    
    // allocate the matrix we will return
    LogicalMatrix rmat(seq.length(), seq.length());
    
    for (int i = 0; i < rmat.nrow(); i++) {
        for (int j = 0; j <= i; j++) {
            
            // check seq equal
            std::string row_seq = as<std::string>(seq[i]);
            std::string col_seq = as<std::string>(seq[j]);
            
            bool is_equal = seqEqual(row_seq, col_seq);
            
            // write to output matrix
            rmat(i,j) = is_equal;
            rmat(j,i) = is_equal;
        }
    }
    
    // Add row and column names
    Rcpp::List dimnames = Rcpp::List::create(seq.attr("names"), 
                                             seq.attr("names"));
    rmat.attr("dimnames") = dimnames;
    
    return rmat;
}


// seqDist
// [[Rcpp::export]]
double seqDistRcpp(std::string seq1, std::string seq2, 
                   NumericMatrix dist_mat) {

    // Check that seq1 and seq2 have same length
    int len_seq1 = seq1.length();
    int len_seq2 = seq2.length();
    
    if (len_seq1 != len_seq2) {
        throw std::range_error("Sequences of different length.");  
    }
    
    int len_seqs = len_seq1;
    
    List dist_mat_dims = dist_mat.attr("dimnames");
    //print (dist_mat_dims);
    CharacterVector dist_mat_rownames = dist_mat_dims[0];
    CharacterVector dist_mat_colnames = dist_mat_dims[1];
    int num_rows = dist_mat_rownames.size();
    int num_cols = dist_mat_colnames.size();
    
    List row_key_idx;
    List col_key_idx;
    
    std::map<std::string, int> rows_map;
    std::map<std::string, int> cols_map;
    
    for (int i = 0; i < num_rows; i++)
    {
        //const char *this_col = dist_mat_colnames[i].c_str();
        std::string this_row = as<std::string>(dist_mat_rownames[i]);
        rows_map[this_row] = i;
    }  
    
    for (int i = 0; i < num_cols; i++)
    {
        //const char *this_col = dist_mat_colnames[i].c_str();
        std::string this_col = as<std::string>(dist_mat_colnames[i]);
        cols_map[this_col] = i;
    } 
    
    int d_seen = 0;
    int indels = 0;
    // sum(d[d>0])
    double d_sum = 0;
    
    for (int i = 0; i < len_seqs; i++)
    {
        // find row index
        int row_idx;
        char row_char = (char)seq1[i];
        std::string row_string;
        row_string+=row_char;
        auto search_row = rows_map.find(row_string);
        if(search_row != rows_map.end()) {
            row_idx = search_row->second;
        }
        else {
            throw std::range_error("Character not found in dist_mat.");  
        }
        
        // find col index
        int col_idx;
        char col_char = (char)seq2[i];
        std::string col_string;
        col_string+=col_char;
        auto search_col = cols_map.find(col_string);
        if(search_col != cols_map.end()) {
            col_idx = search_col->second;
        }
        else {
            throw std::range_error("Character not found in dist_mat.");  
        }    
        
        // distance for current i
        double d_i = dist_mat(row_idx, col_idx);
        
        if (d_i > 0){
            // Sum distance
            d_sum = d_sum + d_i;
        } 
        else if ( (d_i == -1 ) &  (d_seen != -1) )
        {
            // Count indel
            indels++;
        }  
        d_seen = d_i;
    }
    
    double distance = d_sum + indels;
    return (distance);
}


// pairwiseDist
// [[Rcpp::export]]
NumericMatrix pairwiseDistRcpp(StringVector seq, NumericMatrix dist_mat) {
    // allocate the matrix we will return
    NumericMatrix rmat(seq.length(), seq.length());
    
    for (int i = 0; i < rmat.nrow(); i++) {
        for (int j = 0; j < i; j++) {
            
            // check seq equal
            std::string row_seq = as<std::string>(seq[i]);
            std::string col_seq = as<std::string>(seq[j]);
            
            double distance = seqDistRcpp(row_seq, col_seq, dist_mat);
            
            // write to output matrix
            rmat(i,j) = distance;
            rmat(j,i) = distance;
        }
    }
    
    // Add row and column names
    Rcpp::List dimnames = Rcpp::List::create(seq.attr("names"), 
                                             seq.attr("names"));
    rmat.attr("dimnames") = dimnames;
    return rmat;
}


// nonsquareDist
// [[Rcpp::export]]
NumericMatrix nonsquareDistRcpp(StringVector seq, NumericVector indx, NumericMatrix dist_mat)
{
    // defien variables
    int m, n, i, j;
    std::string row_seq, col_seq;
    // extract the sizes. Note: This should be satisfied (n<=m)
    m = indx.size(); //number of rows
    n = seq.size();  //number of columns
    // allocate the main matrix
    NumericMatrix rmat(m,n);
    std::fill(rmat.begin(), rmat.end(), NA_REAL);
    // sort and push indices back by 1 to match c++ indexing
    std::sort(indx.begin(), indx.end());
    indx = indx - 1;
    // find the position of the column ids in the indx vector
    NumericVector pos(n);
    for (j = 0; j < n; j++) {
        pos[j] = std::find(indx.begin(), indx.end(), j) - indx.begin();
    }
    // begin filling rmat
    for (i = 0; i < m; i++) {
        row_seq = as<std::string>(seq[indx[i]]);     //row sequence 
        for (j = 0; j < n; j++) {
            if (!R_IsNA(rmat(i,j))) continue;
            if (indx[i] == j) rmat(i,j) = 0; 
            else {
                col_seq = as<std::string>(seq[j]); //col sequence
                rmat(i,j) = seqDistRcpp(row_seq, col_seq, dist_mat);
                if (pos[j] < m) rmat(pos[j],indx[i]) = rmat(i,j);
            }
        }
    }
    // Add row and column names
    StringVector subSeq = seq[indx];
    Rcpp::List dimnames = Rcpp::List::create(subSeq.attr("names"),      //rownames
                                             seq.attr("names"));  //colnames
    rmat.attr("dimnames") = dimnames;
    // return matrix
    return rmat;
}
