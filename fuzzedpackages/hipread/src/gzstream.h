#ifndef HIPREAD_GZSTREAM_H_
#define HIPREAD_GZSTREAM_H_

#include <zlib.h>
#include <iostream>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

// Adapted from https://stackoverflow.com/questions/21426427/handling-large-gzfile-in-c
// I wish I could use the boost libraries but this makes me think that's not so easy:
// https://github.com/eddelbuettel/bh/issues/49

class GzStream {
private:
  std::string filename_;
  gzFile file;
  char* buffer;
  char* cur;
  char* end;
  void fillBuffer();
  bool done;
  size_t buffer_size;

public:
  GzStream(std::string filename) : filename_(filename), done(false) {
    cur = nullptr;
    file = gzopen(filename.c_str(), "rb");
    buffer_size = 1048576; // start large so that we don't have to copy as much (1048576)
    buffer = new char[buffer_size];
    fillBuffer();
  }
  ~GzStream() {
    cur = nullptr;
    end = nullptr;
    delete[] buffer;
    if (gzclose(file) != Z_OK) Rcpp::stop("Could not close file");
  }
  bool getLine(const char* &line_start, const char* &line_end);
  bool isDone();
  size_t getTotalSizeEstimate();
  size_t getProgress();
  void skipBOM();
  void reset();
};
#endif
