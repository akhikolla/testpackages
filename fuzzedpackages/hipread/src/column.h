// Adapted from readr's Collector.h

#ifndef HIPREAD_COLUMN_H_
#define HIPREAD_COLUMN_H_

#include "boost.h"
#include "iconv.h"
#include <Rcpp.h>

class Column;
typedef boost::shared_ptr<Column> ColumnPtr;

class Column {
protected:
  Rcpp::RObject values_;
  int n_;
  std::vector<std::string> failure_values_;
  std::vector<int> failure_rows_;
  int failure_count_;

public:
  Column(SEXP values) :
    values_(values), n_(0), failure_count_(0){}
  virtual ~Column() {}

  virtual void setValue(int i, const char* x_start, const char* x_end) = 0;

  virtual std::string getType();

  int size() {
    return n_;
  }

  virtual void resize(int n) = 0; // resize is virtual so int/dbl can update pointers
  static ColumnPtr create(std::string type, Rcpp::List var_opts, Iconv* pEncoder_);

  Rcpp::RObject vector() {
    return values_;
  }

  void add_failure(int line_number, const char* x_start, const char* x_end);

  bool has_failures() {
    return failure_count_ > 0;
  }

  std::string describe_failures(std::string var_name);


};

class ColumnCharacter : public Column {
private:
  bool trim_ws;
  Iconv* pEncoder_;
public:
  ColumnCharacter(Rcpp::List opts_, Iconv *pEncoder) : Column(Rcpp::CharacterVector()) {
    trim_ws = opts_["trim_ws"];
    pEncoder_ = pEncoder;
  }
  ~ColumnCharacter() {}
  void resize(int n);
  void setValue(int i, const char* x_start, const char* x_end);
  std::string getType() {return "character";}
};

class ColumnDouble : public Column {
private:
  int imp_dec;
  double imp_dec_base;
  double *valuepointer;
public:
  ColumnDouble(Rcpp::List opts_) : Column(Rcpp::DoubleVector()) {
    // Solaris doesn't have (int, int) method for pow
    imp_dec = opts_["imp_dec"];
    imp_dec_base = static_cast<double>(std::pow(10.0, static_cast<float>(opts_["imp_dec"])));
  }
  ~ColumnDouble() {}
  void resize(int n);
  void setValue(int i, const char* x_start, const char* x_end);
  std::string getType() {return "double";}
};


class ColumnInteger : public Column {
private:
  int *valuepointer;
public:
  ColumnInteger(Rcpp::List opts_) : Column(Rcpp::IntegerVector()) {}
  ~ColumnInteger() {}
  void resize(int n);
  void setValue(int i, const char* x_start, const char* x_end);
  std::string getType() {return "integer";}
};

std::vector<ColumnPtr> createAllColumns(Rcpp::CharacterVector types, Rcpp::List var_opts, Iconv* pEncoder_);
void resizeAllColumns(std::vector<ColumnPtr>& columns, int n);
Rcpp::RObject columnsToDf(std::vector<ColumnPtr> columns, Rcpp::CharacterVector names, int n);

#endif
