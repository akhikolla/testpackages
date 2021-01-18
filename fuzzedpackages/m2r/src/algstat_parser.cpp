#include <Rcpp.h>

#include "algstat_parser.h"

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//






algstat_tokenizer::algstat_tokenizer() {
}


std::vector<std::string> algstat_tokenizer::tokenize(std::string &s) {
  std::vector<std::string> m2operators = this->operators();
  std::vector<std::string> m2symbolchars = this->symbol_chars();
  std::vector<std::string> tokens;

  std::vector<std::string> operatorstartchars(m2operators.size());
  for(size_t i = 0;i < m2operators.size();i++) {
    operatorstartchars[i] = m2operators[i][0];
  }

  for(size_t i = 0;i < s.length();) {
    if (s[i] == ' ' || s[i] == '\t') {
      i++;
    } else if (std::find(m2symbolchars.begin(), m2symbolchars.end(), std::string(&s[i], 1)) != m2symbolchars.end()) {
      size_t start = i;

      i++;
      while (i < s.length() && std::find(m2symbolchars.begin(), m2symbolchars.end(), std::string(&s[i], 1)) != m2symbolchars.end()) {
        i++;
      }

      tokens.push_back(s.substr(start, i - start));
    } else if (std::find(operatorstartchars.begin(), operatorstartchars.end(), std::string(&s[i], 1)) != operatorstartchars.end()) {
      // substr() is smart enough to not index past the end of the string
      for (size_t j = 0;j < m2operators.size();j++) {
        std::string op = m2operators[j];
        if (i + op.length() <= s.length() && op == s.substr(i, op.length())) {
          tokens.push_back(op);
          i = i + op.length();
          break;
        }
      }
    } else if (s[i] == '"') {
      i++;
      size_t start = i;
      for(;i < s.length() && s[i] != '"';i++) {
        if (s[i] == '\\')
          i++;
      }

      tokens.push_back("\"");
      tokens.push_back(s.substr(start, i - start));
      tokens.push_back("\"");
      i++;
    } else if (s[i] == '\n') {
      tokens.push_back("");
    }
    else
      i++;
  }

  return tokens;
}









algstat_parser::algstat_parser(algstat_parser_factory *factory)
  : factory(factory) {
}

algstat_parser::~algstat_parser() {
}

size_t algstat_parser::get_next_index() {
  return next_index;
}



algstat_list_parser::algstat_list_parser(algstat_parser_factory *factory, const std::string &open_char, const std::string &close_char, bool unlist_if_length_one, const std::string &type_name)
  : algstat_parser(factory) {
  this->open_char = open_char;
  this->close_char = close_char;
  this->unlist_if_length_one = unlist_if_length_one;
  this->type_name = type_name;
}

List algstat_list_parser::parse(std::vector<std::string> &tokens, size_t start) {

  List ret = List::create();
  size_t i = start + 1;

  error_on_fail(tokens[start] == this->open_char, "Parsing error: malformed " + type_name);

  if (tokens[i] == this->close_char) {
    i++;
  } else {
    while(true) {

      algstat_parser *parser = factory->create_parser(tokens);
      List elem = parser->parse(tokens, i);
      ret.push_back(elem);
      i = parser->get_next_index() + 1;
      delete parser;

      if (tokens[i] == this->close_char) {
        break;
      }

      error_on_fail(tokens[i] == ",", "Parsing error: malformed " + type_name);

      i++;
      error_on_fail(i < tokens.size(), "Parsing error: malformed " + type_name);

    }
  }

  ret.attr("class") = CharacterVector::create(this->type_name, "m2");

  if (ret.size() == 1) {
    ret = ret[0];
  }

  this->next_index = i;
  return ret;
}









algstat_string_parser::algstat_string_parser(algstat_parser_factory *factory, const std::string &quote_str, const std::string &type_name)
  : algstat_parser(factory) {
  this->quote_str = quote_str;
  this->type_name = type_name;
}

List algstat_string_parser::parse(std::vector<std::string> &tokens, size_t start) {

  List ret = List::create();

  error_on_fail(start + 3 > tokens.size(), "Parsing error: malformed string.")
  error_on_fail(tokens[start+2] == "\"", "Parsing error: malformed string.")

  ret = tokens[start+1];
  ret.attr("class") = CharacterVector::create("m2_string", "m2");

  this->next_index = start + 3;
  return ret;
}



