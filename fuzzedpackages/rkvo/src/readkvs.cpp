#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <Rcpp.h>
using namespace Rcpp;


/**
 * Consumes a string, a delimiter and a vector as a container for
 * delimited values.
 *
 * This function is taken from the following stackoverflow thread:
 * http://stackoverflow.com/questions/21377022/
 *
 * Note that this function could be implemented with a tokenizer. A
 * tokenizer approach would be probably a lot faster and the
 * delimiters could be escaped in key/value fields. This
 * implementation here is not concerned with escape delimiters.
 */
std::vector<std::string> &split(const std::string &s,
                                char delim,
                                std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


/**
 * Reads the file at `path` line by line, splits first key/value pairs
 * by `delim_tuple` and then key/value fileds by `delim_field`, and
 * returns a list of list of key/values.
 *
 * Note that, same keys and corresponding values in the same
 * observation are returned as vectors.
 */
// [[Rcpp::export]]
List readkvs(String path, SEXP delim_tuple, SEXP delim_field) {
  // Declare and initialize the input file stream:
  std::ifstream infile(path.get_cstring());

  // Declare the line:
  std::string line;

  // Get delimiters:
  char delim_tupleC = as<std::string>(delim_tuple)[0];
  char delim_fieldC = as<std::string>(delim_field)[0];

  // Declare the return value:
  List retval;

  // Iterate over the lines (observations):
  while (std::getline(infile, line)) {
    // Declare kev/value pairs vector:
    std::vector<std::string> pairs;

    // Split key/values:
    split(line, delim_tupleC, pairs);

    // Declare R (named-)list for key/value vairs for this
    // observation:
    List ll;

    // Iterate over pairs:
    for (int i=0; i < pairs.size(); i++) {
      // Declare key/value pair vector:
      std::vector<std::string> kv;

      // Split key/value:
      split(pairs.at(i), delim_fieldC, kv);

      // Create a named list for key-value pair and push to the
      // observation list. Note that we are checking same keys and if
      // there are same keywords, values are returned in the same
      // vector:
      if (ll.containsElementNamed(kv.at(0).c_str())) {
        CharacterVector cv(as<CharacterVector>(ll[kv.at(0)]));
        cv.push_back(kv.at(1));
        ll[kv.at(0)] = cv;
      }
      else {
        ll.push_back(CharacterVector(kv.at(1)), kv.at(0));
      }
    }

    // Push the observation to the return value:
    retval.push_back(ll);
  }

  // Done, return with the list of observations as key/value pairs:
  return retval;
}
