// Adapted from readr's source
#ifndef HIPREAD_DATASOURCE_H_
#define HIPREAD_DATASOURCE_H_

#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include "boost.h"
#include "gzstream.h"

class DataSource;
typedef boost::shared_ptr<DataSource> DataSourcePtr;


class DataSource {
  std::string filename_;
public:
  DataSource(std::string filename) : filename_(filename){}
  virtual ~DataSource(){}
  virtual void getLine(const char* &start, const char* &end) = 0;
  virtual bool isDone() = 0;
  virtual std::pair<double, size_t> progress_info() = 0;
  virtual void skipLines(int skip);
  virtual void reset() = 0;
};


class FileDataSource : public DataSource {
private:
  std::string filename_;
  size_t total_size_;
  boost::interprocess::file_mapping fm_;
  boost::interprocess::mapped_region mr_;
  char* file_begin;
  char* file_end;
  char* cur_begin;
  char* cur_end;
  void skipBOM();

public:
  FileDataSource(std::string filename) : DataSource(filename){
    try {
      fm_ = boost::interprocess::file_mapping(
        filename.c_str(), boost::interprocess::read_only);
      mr_ = boost::interprocess::mapped_region(
        fm_, boost::interprocess::read_private);
    } catch (boost::interprocess::interprocess_exception& e) {
      Rcpp::stop("Cannot read file %s: %s", filename, e.what());
    }
    total_size_ = mr_.get_size();
    file_begin = static_cast<char*>(mr_.get_address());
    file_end = file_begin + total_size_;
    cur_begin = file_begin;
    cur_end = nullptr;
    skipBOM();
  }
  ~FileDataSource() {
    file_end = nullptr;
    file_begin = nullptr;
    cur_begin = nullptr;
    cur_end = nullptr;
  }
  void getLine(const char* &start, const char* &end);
  bool isDone();
  std::pair<double, size_t> progress_info();
  void reset();
};



class GzFileDataSource : public DataSource {
private:
  std::string filename_;
  size_t total_size_;
  GzStream *data_;
  size_t get_size();
  void skipBOM();
public:
  GzFileDataSource(std::string filename) : DataSource(filename) {
    data_ = new GzStream(filename);
    total_size_ = get_size();
    skipBOM();
  }
  ~GzFileDataSource() {
    if (data_) delete data_;
  }
  void getLine(const char* &start, const char* &end);
  bool isDone();
  std::pair<double, size_t> progress_info();
  void reset();
};

DataSourcePtr newDataSource(std::string filename, bool isCompressed);

// It seems a little silly to have both the previous interface (using boost shared pointers)
// and this one (using Rcpp's external pointers), but the readr-based code uses the
// former, while the yield methods need external pointers. Not sure if there's a better
// way to reconcile
Rcpp::XPtr<DataSource> newXptrDataSource(std::string filename, bool isCompressed);

#endif
