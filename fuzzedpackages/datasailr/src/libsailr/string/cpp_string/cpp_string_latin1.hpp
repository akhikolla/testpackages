#ifndef CPP_STRING_LATIN1_HPP
#define CPP_STRING_LATIN1_HPP

#include <string>

std::string* cpp_string_new_unescaped_string_latin1( std::string* ori_str );
std::string* cpp_string_subset_latin1 (std::string* ori_str, size_t from_idx , size_t to_idx );
// index starts from zero.


#endif /* CPP_STRING_LATIN1_HPP */
