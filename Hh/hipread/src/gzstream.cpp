#include <Rcpp.h>
using namespace Rcpp;
#include <zlib.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "gzstream.h"


void GzStream::fillBuffer() {
  int err;
  char* offset;
  if (cur == nullptr) {
    offset = buffer;
  } else {
    if (cur == buffer) { // Didn't make any progress in buffer, so double it's size
      char* buffer_current = buffer;
      buffer = new char[buffer_size * 2];
      std::copy_n(buffer_current, buffer_size, buffer);
      buffer_size *= 2;
      delete[] buffer_current;
    }
    offset = std::copy(cur, end, buffer);
  }
  size_t overflow = static_cast<size_t>(offset - buffer);
  if (overflow >= buffer_size) { // Shouldn't ever happen, but want to guard against
    stop("Could not create large enough buffer for gzip file.");
  }

  int len = gzread(file, offset, static_cast<unsigned int>(buffer_size - overflow));

  if (len < 0) stop(gzerror(file, &err));

  cur = buffer;
  end = offset + len;
}

bool GzStream::getLine(const char* &line_start, const char* &line_end) {
  if (isDone()) return false;
  char* eol = std::find(cur, end, '\n');

  while ((eol >= end) && (!gzeof(file))) {
    fillBuffer();
    eol = std::find(cur, end, '\n');
  }
  if (gzeof(file) && (eol >= end)) {
    done = true;
    line_start = cur;
    line_end = end;
    cur = end;
    return true;
  }

  line_start = cur;
  line_end = eol;
  cur = eol + 1;
  return true;
}

bool GzStream::isDone() {
  return done && cur == end;
}

size_t GzStream::getTotalSizeEstimate() {
  std::ifstream raw_(filename_);
  raw_.seekg(0, std::ifstream::end);
  size_t total_file_bytes = static_cast<size_t>(raw_.tellg());

  return total_file_bytes;
}

size_t GzStream::getProgress() {
  size_t out = static_cast<size_t>(gzoffset(file));
  return out;
}

void GzStream::skipBOM() {
  /* Unicode Byte Order Marks
   https://en.wikipedia.org/wiki/Byte_order_mark#Representations_of_byte_order_marks_by_encoding
  00 00 FE FF: UTF-32BE
  FF FE 00 00: UTF-32LE
  FE FF:       UTF-16BE
  FF FE:       UTF-16LE
  EF BB BF:    UTF-8
  */

  switch (cur[0]) {
  // UTF-32BE
  case '\x00':
    if (end - cur >= 4 && cur[1] == '\x00' && cur[2] == '\xFE' &&
        cur[3] == '\xFF') {
      cur += 4;
    }
    break;

    // UTF-8
  case '\xEF':
    if (end - cur >= 3 && cur[1] == '\xBB' && cur[2] == '\xBF') {
      cur += 3;
    }
    break;

    // UTF-16BE
  case '\xfe':
    if (end - cur >= 2 && cur[1] == '\xff') {
      cur += 2;
    }
    break;

  case '\xff':
    if (end - cur >= 2 && cur[1] == '\xfe') {

      if (end - cur >= 4 && cur[2] == '\x00' && cur[3] == '\x00') { // UTF-32 LE
        cur += 4;
      } else { // UTF-16 LE
        cur += 2;
      }

    }
    break;
  }
}

void GzStream::reset() {
  cur = nullptr;
  done = false;
  gzrewind(file);
  fillBuffer();
}
