#ifndef M2_PARSER_H
#define M2_PARSER_H

#include "algstat_parser.h"


std::vector<std::string> m2_tokenize_cpp(std::string &s);
Rcpp::List m2_parse_internal_cpp(std::vector<std::string> &tokens);


class m2_tokenizer : public algstat_tokenizer {
public:
  virtual std::vector<std::string> symbol_chars();
  virtual std::vector<std::string> operators();

protected:

};






/*
class m2_parser_factory : public algstat_parser_factory {
public:
  virtual algstat_parser *create_parser(std::vector<std::string> &tokens);
};


class m2_parser : public algstat_parser {
public:
  m2_parser(algstat_parser_factory *factory);

  virtual Rcpp::List parse(std::vector<std::string> &tokens, size_t start);

protected:

};
*/


#endif



