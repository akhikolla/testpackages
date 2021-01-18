#ifndef CHUNKER_H_
#define CHUNKER_H_

#include <Rcpp.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;

namespace _chunkR {

class chunker {

public:
  // constructors & destructor-----
  // matrix constructor
  chunker(const std::string path, char sep, bool quoted, 
          bool has_colnames, bool has_rownames,
          size_t chunksize, StringVector column_types);
  // data.frame constructor
	chunker(const std::string path, char sep,  bool quoted,
         bool has_colnames, bool has_rownames, size_t chunksize);
	virtual ~chunker();
	
	// next chunk ------------------
	bool next_chunk();
	bool next_chunk_matrix();
	bool next_chunk_df();
	
  // setters ---------------------
  void set_offset();
  void count_lines();
	void set_colnames();
	std::vector<std::string> set_generic_rownames(std::string what, size_t start_from, size_t rownumber);
	std::vector<std::string> set_generic_colnames(std::string what, size_t start_from, size_t colnumber);
	
	// getters --------------------
	StringMatrix get_matrix();
	DataFrame get_dataframe();
	StringVector get_colnames();
	size_t get_total();
	size_t get_completed();
	const std::string get_type();	
	
	// auxiliary ------------------
	inline List mixed_list(std::vector<int> x,  int howmuch);
	
	// validators ----------------
	bool is_valid_chunker();

private:
	const std::string path;
	const char sep;
	bool quoted;
	bool has_colnames;
	const bool has_rownames;
	const size_t chunksize;
	const std::string output_format;
	std::vector<int> col_types;

	size_t n_row;
	size_t n_col;
	std::vector<std::string> rnames;
	std::vector<std::string> cnames;
	size_t n_lines;
	
	std::ifstream line_container;
	size_t pointer_position;

	std::string* line;
	std::string* element;
	size_t lines_completed;
	std::vector<std::string> word;
	std::string validation_state;
	
	int offset;
	char eof;
  
	struct chunkInfo {
	  StringMatrix m;
	  DataFrame df;
	} data_chunk;

};

} /* namespace _chunkR */

#endif /* CHUNKER_H_ */
