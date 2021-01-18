#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>

// [[Rcpp::export]]
std::vector<double> GetLineLocations(std::string &filename) {
  std::vector<double> indices;
  std::ifstream infile;
  std::string line;

  infile.open(filename.c_str());
  if (!infile.good()) {
    Rcpp::stop("Unable to open file");
//    Rcpp::Rcerr << "Unable to open file" << std::endl;
//    return indices;
  }

  while(!infile.fail() && indices.size() < 100) {
    indices.push_back((double)infile.tellg());
    std::getline(infile, line);
  }

  infile.close();
  return indices;
}
