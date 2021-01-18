#include <math.h>
#include <fstream>
#include <inttypes.h>

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
int serialMat(NumericMatrix x, StringVector filename, bool append) {
  // 8-bit buffer
  uint8_t buffer = 0;

  int count = 0;
  int64_t rows = x.nrow(), cols = x.ncol();

  // Open file to write the matrix
  std::ofstream ofile(filename[0], std::ios::binary | std::ios::app);
  if (!ofile)
    return -1;

  // Write header
  if (!append) {
    ofile.write("epinetr", sizeof("epinetr"));
    ofile.write(reinterpret_cast<const char*> (&rows), sizeof(rows));
    ofile.write(reinterpret_cast<const char*> (&cols), sizeof(cols));
  }

  // Write the matrix
  for (int64_t i=0; i < rows; i++)
    for (int64_t j=0; j < cols; j++) {
      buffer <<= 1;
      if (x(i, j))
        buffer |= 1;
      count++;
      if (count == 8 || j + 1 == cols) {
        ofile.write(reinterpret_cast<const char*>(&buffer), sizeof(buffer));
        buffer = 0;
        count = 0;
      }
    }

  ofile.close();
  if (!ofile)
    return -1;

  // If we just appended the data, update number of rows in the file
  if (append) {
    std::fstream file(filename[0], std::ios::binary | std::ios::in | std::ios::out);
    if (!file)
      return -1;

    // Read current number of rows
    file.seekg(8, std::ios::beg);
    int64_t crows;
    file.read(reinterpret_cast<char *>(&crows), sizeof(crows));

    // Write new number of rows
    file.seekp(8, std::ios::beg);
    rows += crows;
    file.write(reinterpret_cast<const char*> (&rows), sizeof(rows));

    file.close();
    if (!file)
      return -1;
  }

  return 0;
}


// [[Rcpp::export]]
NumericMatrix getSerialMat(StringVector filename) {
  uint8_t buffer;
  std::string head = "";

  // Open the file
  std::ifstream ifile(filename[0], std::ios::binary | std::ios::in);
  if (!ifile)
    return R_NilValue;

  // Read first eight characters and check for "epinetr"
  ifile.seekg(0, std::ios::beg);
  while ((buffer = ifile.get()) != '\0') {
    head += buffer;
  }
  if (head.compare("epinetr"))
    return R_NilValue;

  // Get number of rows
  int64_t rows;

  ifile.read(reinterpret_cast<char *>(&rows), sizeof(rows));

  // Get number of columns
  int64_t cols;

  ifile.read(reinterpret_cast<char *>(&cols), sizeof(cols));

  // Create matrix
  NumericMatrix m(rows, cols);

  int bytesPerRow = ceil(float(cols) / 8);

  // Read in matrix
  for (int64_t i = 0; i < rows; i++)
    for (int64_t j = 1; j <= bytesPerRow; j++) {
      ifile.read(reinterpret_cast<char *>(&buffer), sizeof(buffer));

      int start = 1;
      if (j == bytesPerRow)
        start = 8 - (cols % 8) + 1;

      for (int k = start; k <= 8; k++) {
        if (buffer & 1)
          m(i, j*8 - k) = 1;
        buffer >>= 1;
      }
    }

  ifile.close();
  if (!ifile)
    return R_NilValue;

  return m;
}
