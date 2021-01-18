#ifndef ALGSTAT_PARSER_H
#define ALGSTAT_PARSER_H


#define error_on_fail(t, e) if (!(t)) throw (e);


// splits a string containing M2 code into tokens to ease parsing
// places an empty string between each line
class algstat_tokenizer {
public:
  algstat_tokenizer();
  virtual std::vector<std::string> symbol_chars() = 0;
  virtual std::vector<std::string> operators() = 0;

  virtual std::vector<std::string> tokenize(std::string &s);

private:

};


typedef class algstat_parser algstat_parser;

class algstat_parser_factory {
public:
  virtual algstat_parser *create_parser(std::vector<std::string> &tokens) = 0;
};


class algstat_parser {
public:
  algstat_parser(algstat_parser_factory *factory);
  virtual ~algstat_parser();

  virtual Rcpp::List parse(std::vector<std::string> &tokens, size_t start) = 0;
  virtual size_t get_next_index();

protected:
  algstat_parser_factory *factory;
  size_t next_index;
};







class algstat_list_parser : public algstat_parser {
public:
  algstat_list_parser(algstat_parser_factory *factory, const std::string &open_char, const std::string &close_char, bool unlist_if_length_one, const std::string &type_name);

  virtual Rcpp::List parse(std::vector<std::string> &tokens, size_t start);

protected:
  std::string open_char;
  std::string close_char;
  bool unlist_if_length_one;
  std::string type_name;
};


class algstat_string_parser : public algstat_parser {
public:
  algstat_string_parser(algstat_parser_factory *factory, const std::string &quote_str, const std::string &type_name);

  virtual Rcpp::List parse(std::vector<std::string> &tokens, size_t start);

protected:
  std::string quote_str;
  std::string type_name;
};




#endif
