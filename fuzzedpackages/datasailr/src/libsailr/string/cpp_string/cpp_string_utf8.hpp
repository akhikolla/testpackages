#ifndef CPP_STRING_UTF8_HPP
#define CPP_STRING_UTF8_HPP

#include <string>
#include "utfcpp/utf8.h"

std::string* cpp_string_subset_utf8(std::string* ori_str, int from_idx, int to_idx);
std::string* cpp_string_new_unescaped_string_utf8( std::string* ori_str );


#endif /* CPP_STRING_UTF8_HPP */
