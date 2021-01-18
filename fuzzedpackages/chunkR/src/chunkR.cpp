// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <chunkR.h>

using namespace Rcpp;


StringMatrix NULLstringm(0);
const StringMatrix empty_stringm = NULLstringm;

DataFrame NULLdf = DataFrame::create();
const DataFrame empty_df = NULLdf;

std::vector<std::string> autovector;
#define auto_vector autovector


namespace _chunkR {

//--------------------------------------------
// CONSTRUCTORS & DESTRUCTOR
//--------------------------------------------

//' chunker, matrix-constructor
//' @keywords internal

chunker::chunker(const std::string path, char sep, bool quoted, 
                 bool has_colnames, bool has_rownames, size_t chunksize) :
path(path), sep(sep), quoted(quoted), has_colnames(has_colnames),
has_rownames(has_rownames), chunksize(chunksize), output_format("matrix"), n_row(0), n_col(0),
rnames([&chunksize] {std::vector<std::string> out; out.reserve(chunksize); return out;}()), 
cnames([] {std::vector<std::string> out; return out;}()), 
pointer_position(0), line(new std::string), element(new std::string), 
lines_completed(0), word(auto_vector), validation_state("") {
  
  set_offset();
  data_chunk.m = empty_stringm;
  data_chunk.df = empty_df; 
  
  Rcout << "New chunker object\n";
  Rcout << "Path: " << path << std::endl;
  set_colnames();
}


//' chunker, dataframe-constructor
//' @keywords internal

chunker::chunker(const std::string path, char sep, bool quoted, 
                 bool has_colnames, bool has_rownames, size_t chunksize,
                 StringVector column_types) :
path(path), sep(sep),  quoted(quoted), has_colnames(has_colnames),
has_rownames(has_rownames), chunksize(chunksize), output_format("data.frame"), n_row(0), n_col(0),
rnames([&chunksize] {std::vector<std::string> out; out.reserve(chunksize); return out;}()), 
cnames([] {std::vector<std::string> out; return out;}()), 
pointer_position(0), line(new std::string), element(new std::string), 
lines_completed(0), word(auto_vector), validation_state("") {
  
  // set offset
  set_offset();
  // count lines in file
  count_lines();
  
  data_chunk.m = empty_stringm;
  data_chunk.df = empty_df; 
  
  // configure for data frames
  size_t df_size = column_types.size();
  if(df_size > 0) {
    std::vector<int>temp;
    temp.reserve(df_size);
    for(size_t i = 0; i < df_size; ++i) {
      if(column_types[i] == "character") {
        temp.push_back(0);
      } else  if(column_types[i] == "double" || column_types[i] == "numeric" ) {
        temp.push_back(1);
      } else if(column_types[i] == "integer") {
        temp.push_back(2);
      } else if(column_types[i] == "logical") {
        temp.push_back(3);
      }
    }
    col_types = temp; 
  }
  
  
  Rcout << "New chunker object\n";
  Rcout << "Path: " << path << std::endl;
  set_colnames();
}

//' chunker__destructor 
//' @keywords internal

chunker::~chunker() {
  delete line;
  delete element;
}

//--------------------------------------------
// NEXT_CHUNK 
//--------------------------------------------

//' chunker__next_chunk
//' @keywords internal

bool chunker::next_chunk() {
  if(output_format == "matrix") {
    return next_chunk_matrix();
  } else {
    return next_chunk_df();
  }
}


//' chunker__next_chunk_matrix
//' @keywords internal

bool chunker::next_chunk_matrix() {
  if(!is_valid_chunker()) {
    Rcout << validation_state << std::endl;
    return false;
  }
  if (!line_container.eof()) {
    
    line_container.open(path.c_str(), std::ios::binary);
    if (line_container.fail()) {
      std::ostringstream msg;
      msg << "Input file opening failed.\n";
      throw(msg.str());
    }
    
    line_container.seekg(pointer_position + offset);
    
    // for first line read before
    size_t lines_read_chunk = 0;
    size_t initial_lines  = lines_completed;
    
    while (std::getline(line_container, *line, eof)) {
      bool is_name = true;
      std::stringstream ss(*line);

      while (std::getline(ss, *element, sep)) {
	
        if(quoted) {
          element->erase(remove(element->begin(), element->end(), '\"' ), element->end());
        }
        if (is_name && has_rownames) {
          rnames.push_back(*element);
          is_name = false;
          continue;
        }
        word.push_back(*element);
      }
      
      n_row++;
      lines_completed++;
      if (++lines_read_chunk >= chunksize) {
        break;
      }
      line_container.seekg(offset, std::ios::cur);
    }
    

    pointer_position = line_container.tellg();
    
    try {
      if(word.size() != (lines_read_chunk * n_col)) throw(1);
    } catch(int e) {
      Rcout << "Number of elements are not multiple of matrix dimension" << std::endl;
      line_container.seekg(0, std::ios::end);
      word.clear();
      line_container.close();
      rnames.clear();
      data_chunk.m = empty_stringm;
      validation_state = "Invalid object: number of elements are not multiple of matrix dimension";
      return false;
    }
    
    // create output
    StringMatrix output(lines_read_chunk, n_col);
    int k = 0;
    for (size_t i = 0; i < lines_read_chunk; ++i) {
      for (size_t j = 0; j < n_col; ++j) {
        output(i, j) = word[k++];
      }
    }
    
    StringVector cnames_rcpp;
    cnames_rcpp = cnames;
    colnames(output) = cnames_rcpp;
    
    // set rownames
    StringVector rnames_rcpp;
    if (has_rownames) {
      rnames_rcpp = rnames;
      rnames.clear();
    } else {
      rnames_rcpp = set_generic_rownames("R", initial_lines + 1, lines_read_chunk);
    }
    rownames(output) = rnames_rcpp;
    
    //cnames.clear();
    word.clear();
    line_container.close();
    
    if (output.nrow() == 0) {
      data_chunk.m = empty_stringm;
      return false;
    }
    
    data_chunk.m = output;
    return true;
    
  } else {
    data_chunk.m = empty_stringm;
    return false;
  }
  
}

//' chunker__next_chunk_df
//' @keywords internal

bool chunker::next_chunk_df() {
  
  if (!line_container.eof()) {
    
    line_container.open(path.c_str(), std::ios::binary);
    if (line_container.fail()) {
      std::ostringstream msg;
      msg << "Input file opening failed.\n";
      throw(msg.str());
    }
    
    line_container.seekg(pointer_position + offset);
    
    // for first line read before
    size_t lines_read_chunk = 0;
    size_t initial_lines  = lines_completed;
    List output;
    
    if(n_lines  - initial_lines >= chunksize) {
      output = mixed_list(col_types, chunksize);
    } else {
      output = mixed_list(col_types, n_lines - initial_lines); 
    }
    
    while (std::getline(line_container, *line, eof)) {
      bool is_name = true;
      
      std::stringstream ss(*line);
      int count_elements = 0;
      while (std::getline(ss, *element, sep)) {
        if (is_name && has_rownames) {
          if(quoted) {
            element->erase(remove(element->begin(), element->end(), '\"' ), element->end());
          }
          rnames.push_back(*element);
          is_name = false;
          continue;
        }
        
        if(col_types[count_elements] == 0) {
	
          if(quoted) {
            element->erase(remove(element->begin(), element->end(), '\"' ), element->end());
          }
          std::string temp_out = *element;
          (Rcpp::as<StringVector>(output[count_elements]))[lines_read_chunk] = temp_out;
        } else if(col_types[count_elements] == 1) {
          double temp_out = std::atof((* element).c_str());
          (Rcpp::as<NumericVector>(output[count_elements]))[lines_read_chunk] = temp_out;
        } else if(col_types[count_elements] == 2 || output_format[count_elements++] == 3) {
          int temp_out = std::atoi((* element).c_str());
          (Rcpp::as<IntegerVector>(output[count_elements])[lines_read_chunk]) = temp_out;
        }
        count_elements++;
      }
      
      n_row++;
      lines_completed++;
      if (++lines_read_chunk >= chunksize) {
        break;
      }
      
      line_container.seekg(offset, std::ios::cur);
    }
    
    pointer_position = line_container.tellg();
    
    
    StringVector cnames_rcpp;
    cnames_rcpp = cnames;
    output.attr("names") = cnames_rcpp;
    
    // set rownames
    StringVector rnames_rcpp;
    if (has_rownames) {
      rnames_rcpp = rnames;
      rnames.clear();
    } else {
      rnames_rcpp = set_generic_rownames("R", initial_lines + 1, lines_read_chunk);
    }
    
    output.attr("row.names") = rnames_rcpp;
    output.attr("class") = "data.frame";
    
    line_container.close();
    
    if (lines_read_chunk == 0) {
      DataFrame output;
      data_chunk.df = empty_df;
      return false;
    }
    
    data_chunk.df = output;
    return true;
    
  } else {
    data_chunk.df = empty_df;
    return false;
  }
  
}


//--------------------------------------------
// SETTERS
//--------------------------------------------

//' chunker__set_colnames 
//' @keywords internal

void chunker::set_colnames() {
  
  line_container.open(path.c_str(), std::ios::binary);
  try {
    if (line_container.fail()) {
      std::ostringstream msg;
      msg << "Input file opening failed.\n";
      throw(msg.str());
    }
  } catch (std::string& stopmsg) {
    Rcout << stopmsg;
  }
  
  if (has_colnames) {
    std::getline(line_container, *line, eof);
    std::stringstream headss(*line);
    while (std::getline(headss, *element, sep)) {
	
      if(quoted) {
        element->erase(remove(element->begin(), element->end(), '\"' ), element->end());
      }
      cnames.push_back(*element);
    }
  }
  pointer_position = line_container.tellg();
  
  // Read first line only to set number of columns
  bool is_name = true;
  std::getline(line_container, *line, eof);
  std::stringstream ss(*line);
  while (std::getline(ss, *element, sep)) {
    if (is_name && has_rownames) {
      is_name = false;
      continue;
    }
    n_col++;
  }
  
  // check consistency between column number determined and # elements in has_colnames
  if (has_colnames) {
    try {
      size_t number_in_has_colnames = cnames.size();
      
      if (number_in_has_colnames != n_col) {
        std::ostringstream msg;
        msg << "Error: Number of strings in has_colnames ("
            << number_in_has_colnames << ") " << "has not " << n_col
            << " elements";
        throw(msg.str());
      }
    } catch (std::string& m) {
      Rcout << m;
    }
  }
  
  if (!has_colnames) {
    cnames = set_generic_colnames("C", 1, n_col);
  }
  
  line_container.close();
}


//' set_generic_rownames
//' @keywords internal

std::vector<std::string> chunker::set_generic_rownames(std::string what, size_t start_from, size_t rownumber) {
  std::ostringstream os;
  std::vector<std::string> output;
  output.reserve(rownumber);
  for (size_t i = start_from; i < start_from + rownumber; ++i) {
    os << what << "_" << i;
    output.push_back(os.str());
    os.str("");
    os.clear();
  }
  return output;
}

//' set_generic_colnames
//' @keywords internal

std::vector<std::string> chunker::set_generic_colnames(std::string what,  
                                                       size_t start_from,
                                                       size_t colnumber) {
  std::ostringstream os;
  std::vector<std::string> output;
  output.reserve(colnumber);
  for (size_t i = start_from; i < start_from + colnumber; ++i) {
    os << what << "_" << i;
    output.push_back(os.str());
    os.str("");
    os.clear();
  }
  return output;
}

//--------------------------------------------
// GETTERS
//--------------------------------------------

// GET_MATRIX & GET_DATAFRAME-----------------

//'get_matrix
//'@keywords internal

StringMatrix chunker::get_matrix() {
  if(!is_valid_chunker()) {
    Rcpp::stop(validation_state); 
  }
  return data_chunk.m;
} 

//'get_dataframe
//'@keywords internal
//'
DataFrame chunker::get_dataframe() {
  if(!is_valid_chunker()) {
    Rcpp::stop(validation_state); 
  }
  return data_chunk.df;
} 

//---------------------------------------------

//' chunker__get_colnames
//' @keywords internal

StringVector chunker::get_colnames() {
  StringVector out(cnames.size());
  out = cnames;
  return out;
}


//' chunker__get_completed
//' @keywords internal

size_t chunker::get_completed() {
  return lines_completed;
}

//' chunker__get_completed
//' @keywords internal

size_t chunker::get_total() {
  return n_lines;
}

//' chunker__get_type
//' @keywords internal

const std::string chunker::get_type() {
  return output_format;
}


//--------------------------------------------
// AUXILIARY
//--------------------------------------------

//' mixed_list
//' @keywords internal

inline List chunker::mixed_list(std::vector<int> x,  int howmuch) {
  List output;
  
  for(size_t i = 0; i < x.size(); ++i) {
    switch (x[i]) {
    case 0:
  {
    StringVector temp(howmuch);
    output.push_back(temp);
    break;
  }
    case 1:
  {
    NumericVector temp1(howmuch);
    output.push_back(temp1);
    break;
  }
    case 2:
  {
    IntegerVector temp2(howmuch);
    output.push_back(temp2);
    break;
  }
    case 3:
  {
    LogicalVector temp3(howmuch);
    output.push_back(temp3);
    break;
  }
    default:
      stop("Unknown type");
    }
  }
  return output;
}

// to create an invalid object, passa text to validation state.
// If validation_state.size() != 0, throw an error
bool chunker::is_valid_chunker() {
  if(validation_state.size() != 0) return false;
  return true;
}

void chunker::set_offset() {
  std::ifstream file(path.c_str(), std::ios::binary);
  char c;
  offset = -1;
  while (file.get(c)) {
    if (c == '\n') {
      offset = 0;
      eof = '\n';
      break;
    }
    if (c == '\r') {
      if (file.get(c)) {
        if (c == '\n') {
          offset = 1;
          eof = '\r';
          break;
        }
      }
      offset = 0;
      eof = '\r';
      break;
    }
  }
  if(offset == -1) {
    Rcpp::stop("Invalid end of line");
  }
}


void chunker::count_lines() {
Rcout << "Counting lines in file... ";
std::ifstream this_file(path.c_str(), std::ios::binary); 
n_lines = std::count(std::istreambuf_iterator<char>(this_file), 
                     std::istreambuf_iterator<char>(), eof);
this_file.close();
if(has_colnames) {
  n_lines = n_lines -1;
}
Rcout << "done.\n";
}

// bool is_space(char x) {
//   if(x == ' ' || x == '\n' || x == '\n' || x == '\v' || x == '\f' || x == '\r') {
//     return true;
//   } 
//   return false;
// }
// 
// bool test_empty(std::string buffer) {
//   for (int i = 0; i< buffer.length(); ++i){
//     if(!is_space(buffer[i])) {
//       return false;
//     } 
//     return true;
//   }
// }


} /* namespace _chunkR */
