#include "readcol.h"
#include "ncols.h"


//' Read one column of a whitespace delimited text file
//'
//' @param fileName Input filepath
//' @param colNum An integer. The target column number
//' @param nSkip An integer. Number of lines to skip in the beginning. 
//' @param maxRowNum An Integer. Maximum number of lines to read
//' @return A vector of strings containing the target column
//' @export
// [[Rcpp::export]]
std::vector< std::string > readcol(std::string fileName, long colNum, long nSkip, long maxRowNum) {
    long nColOfFile = ncols(fileName);
    if(colNum <  1 or colNum > nColOfFile) throw "Column number out of range!";

    std::ifstream inFileStream(fileName.c_str());
    std::vector< std::string > res;
    std::string tmpLine;
    long rowCounter = 0;

    // skip lines at the begining
    for(long i = 0; i < nSkip; i++) {
        getline(inFileStream, tmpLine);
    }

    while(getline(inFileStream, tmpLine)) {
        rowCounter++;
        if(rowCounter > maxRowNum) break;

        std::istringstream lineStream(tmpLine);
        std::string tmpWord;
        for(long i=0; i<colNum; i++) {
            lineStream >> tmpWord;
        }
        res.push_back(tmpWord);
    }

    return res;
}
