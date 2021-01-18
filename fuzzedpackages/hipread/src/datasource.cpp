#include <Rcpp.h>
using namespace Rcpp;
#include <fstream>
#include <iostream>
#include <zlib.h>
#include <algorithm>
#include "boost.h"
#include "datasource.h"

void DataSource::skipLines(int skip) {
  const char* ignore_start;
  const char* ignore_end;
  for (int i = 0; i < skip; ++i) {
    if (isDone()) break;
    getLine(ignore_start, ignore_end);
  }
}

void FileDataSource::getLine(const char* &start, const char* &end) {
  if (cur_end != nullptr) cur_begin = cur_end + 1;
  cur_end = std::find_if(cur_begin, file_end, [](int ch) {
    return ch == '\n';
  });
  if (cur_end > file_end) cur_end = file_end;

  start = cur_begin;
  end = cur_end;
}

bool FileDataSource::isDone() {
  return cur_begin == file_end;
}

void FileDataSource::skipBOM() {
  /* Unicode Byte Order Marks
   https://en.wikipedia.org/wiki/Byte_order_mark#Representations_of_byte_order_marks_by_encoding
  00 00 FE FF: UTF-32BE
  FF FE 00 00: UTF-32LE
  FE FF:       UTF-16BE
  FF FE:       UTF-16LE
  EF BB BF:    UTF-8
  */
  switch (cur_begin[0]) {
  // UTF-32BE
  case '\x00':
    if (file_end - cur_begin >= 4 && cur_begin[1] == '\x00' && cur_begin[2] == '\xFE' &&
        cur_begin[3] == '\xFF') {
      cur_begin += 4;
    }
    break;

    // UTF-8
  case '\xEF':
    if (file_end - cur_begin >= 3 && cur_begin[1] == '\xBB' && cur_begin[2] == '\xBF') {
      cur_begin += 3;
    }
    break;

    // UTF-16BE
  case '\xfe':
    if (file_end - cur_begin >= 2 && cur_begin[1] == '\xff') {
      cur_begin += 2;
    }
    break;

  case '\xff':
    if (file_end - cur_begin >= 2 && cur_begin[1] == '\xfe') {

      // UTF-32 LE
      if (file_end - cur_begin >= 4 && cur_begin[2] == '\x00' && cur_begin[3] == '\x00') {
        cur_begin += 4;
      } else { // UTF-16 LE
        cur_begin += 2;
      }

    }
    break;
  }
}

std::pair<double, size_t> FileDataSource::progress_info() {
  if (isDone()) {
    return std::make_pair(1.0, total_size_);
  } else {
    size_t current = static_cast<size_t>(cur_begin - file_begin);
    return std::make_pair(current / static_cast<double>(total_size_), current);
  }
}

void FileDataSource::reset() {
  cur_begin = file_begin;
  cur_end = nullptr;
  skipBOM();
}

size_t GzFileDataSource::get_size() {
  return data_->getTotalSizeEstimate();
}


void GzFileDataSource::getLine(const char* &start, const char* &end) {
  data_->getLine(start, end);
}

bool GzFileDataSource::isDone() {
  return data_->isDone();
}

void GzFileDataSource::skipBOM() {
  data_->skipBOM();
}

void GzFileDataSource::reset() {
  data_->reset();
  skipBOM();
}

std::pair<double, size_t> GzFileDataSource::progress_info() {
  if (isDone()) {
    return std::make_pair(1.0, total_size_);
  } else {
    size_t current = data_->getProgress() ;
    return std::make_pair(current / static_cast<double>(total_size_), current);
  }
}

DataSourcePtr newDataSource(std::string filename, bool isCompressed) {
  if (isCompressed) {
    return DataSourcePtr(new GzFileDataSource(filename));
  } else {
    return DataSourcePtr(new FileDataSource(filename));
  }
}

Rcpp::XPtr<DataSource> newXptrDataSource(std::string filename, bool isCompressed) {
  if (isCompressed) {
    return Rcpp::XPtr<DataSource>(new GzFileDataSource(filename));
  } else {
    return Rcpp::XPtr<DataSource>(new FileDataSource(filename));
  }
}
