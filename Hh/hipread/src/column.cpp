// Adapted from readr's Collector.cpp

#include <Rcpp.h>
#include "column.h"
#include "string_utils.h"
#include "iconv.h"

#include <Rcpp.h>
using namespace Rcpp;

std::string Column::getType() {
  return "unknown";
}

ColumnPtr Column::create(std::string type, Rcpp::List var_opts, Iconv* pEncoder_) {
  if (type == "character") {
    return ColumnPtr(new ColumnCharacter(var_opts, pEncoder_));
  } else if (type == "double") {
    return ColumnPtr(new ColumnDouble(var_opts));
  } else if (type == "integer") {
    return ColumnPtr(new ColumnInteger(var_opts));
  }

  Rcpp::stop("Unexpected column type '%s'", type);
}


void Column::add_failure(int line_number, const char* x_start, const char* x_end) {
  if (++failure_count_ <= 5) {
    std::string value;
    size_t x_diff = static_cast<size_t>(x_end - x_start);
    value.assign(x_start, x_diff);

    failure_values_.push_back(value);
    failure_rows_.push_back(line_number + 1);
  }
}

std::string Column::describe_failures(std::string var_name) {
  std::ostringstream message;

  message << "In variable '" << var_name << "', could not convert " << failure_count_ <<
    " values to " << getType() << "; Values (and row numbers) of first " <<
      failure_rows_.size() << " failures: ";

  size_t num_failures = failure_rows_.size();
  for (size_t i = 0; i < num_failures; ++i) {
    if (i > 0) {
      message << ", ";
    }

    message << "'" << failure_values_[i] << "' (" <<
      failure_rows_[i] << ")";
  }

  return message.str();
}


void ColumnCharacter::setValue(int i, const char* x_start, const char* x_end) {
  if (trim_ws) IpStringUtils::newtrim(x_start, x_end);
  SET_STRING_ELT(values_, i, pEncoder_->makeSEXP(x_start, x_end));
}

void ColumnDouble::setValue(int i, const char* x_start, const char* x_end) {
  double value;
  IpStringUtils::newtrim(x_start, x_end);
  bool success;
  if (x_start == x_end) {
    success = true;
    value = NA_REAL;
  } else {
    success = parseDouble(x_start, x_end, value);
  }

  if (!success) {
    add_failure(i, x_start, x_end);
    value = NA_REAL;
  } else if (imp_dec != 0) {
    value = value / imp_dec_base;
  }
  valuepointer[i] = value;
}

void ColumnInteger::setValue(int i, const char* x_start, const char* x_end) {
  int value;
  IpStringUtils::newtrim(x_start, x_end);
  bool success;
  if (x_start == x_end) {
    success = true;
    value = NA_INTEGER;
  } else {
    success = parseInteger(x_start, x_end, value);
  }

  if (!success) {
    add_failure(i, x_start, x_end);
    value = NA_INTEGER;
  }
  valuepointer[i] = value;
}

void ColumnCharacter::resize(int n) {
  if (n == n_)
    return;

  if (n > 0 && n < n_) {
    SETLENGTH(values_, n);
    SET_TRUELENGTH(values_, n);
  } else {
    values_ = Rf_lengthgets(values_, n);
  }
  n_ = n;

}

void ColumnDouble::resize(int n) {
  if (n == n_)
    return;

  if (n > 0 && n < n_) {
    SETLENGTH(values_, n);
    SET_TRUELENGTH(values_, n);
  } else {
    values_ = Rf_lengthgets(values_, n);
  }
  n_ = n;
  valuepointer = REAL(values_);
}

void ColumnInteger::resize(int n) {
  if (n == n_)
    return;

  if (n > 0 && n < n_) {
    SETLENGTH(values_, n);
    SET_TRUELENGTH(values_, n);
  } else {
    values_ = Rf_lengthgets(values_, n);
  }
  n_ = n;
  valuepointer = INTEGER(values_);
}

std::vector<ColumnPtr> createAllColumns(CharacterVector types, Rcpp::List var_opts, Iconv* pEncoder_) {
  int num_cols = static_cast<int>(types.size());
  std::vector<ColumnPtr> out;

  for (int i = 0; i < num_cols; ++i) {
    out.push_back(Column::create(as<std::string>(types[i]), var_opts[i], pEncoder_));
  }

  return out;
}

void resizeAllColumns(std::vector<ColumnPtr>& columns, int n) {
  size_t num_cols = columns.size();

  for (size_t i = 0; i < num_cols; ++i) {
    columns[i]->resize(n);
  }
}

static Function as_tibble("as_tibble", Environment::namespace_env("tibble"));

RObject columnsToDf(std::vector<ColumnPtr> columns, Rcpp::CharacterVector names, int n) {
  size_t num_vars = columns.size();

  List out(num_vars);
  for (size_t i = 0; i < num_vars; ++i) {
    if (columns[i]->has_failures()) {
      std::string message = columns[i]->describe_failures(
              Rcpp::as<std::string>(names[static_cast<long>(i)]));
      Rf_warning(message.c_str());
    }
    out[static_cast<long>(i)] = columns[i]->vector();
  }
  out.attr("names") = names;
  out.attr("class") = CharacterVector::create("tbl_df", "tbl", "data.frame");
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -(n));

  return out;
}
