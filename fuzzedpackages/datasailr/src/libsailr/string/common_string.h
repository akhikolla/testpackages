#ifndef COMMON_STRING_H
#define COMMON_STRING_H

#include <stdlib.h>
#include "cpp_string/cpp_string.hpp"

typedef cpp_object string_object;
typedef char c_char;

string_object* string_new(const c_char* str );
string_object* string_new_with_len( const c_char* str, int len);
string_object* string_new_unescaped_and_delete_ori(string_object* str, const char* encoding);
string_object* string_clone( string_object* );
string_object* string_new_int2str( int num); 
string_object* string_new_double2str( double num);
/*
const c_char* string_int2cstr( int num); 
const c_char* string_double2cstr( double num); 
*/
string_object* string_strip( cpp_object* ); 
string_object* string_lstrip( cpp_object* ); 
string_object* string_rstrip( cpp_object* ); 
const c_char* string_read( string_object* str);
string_object* string_concat( string_object* str1, string_object* str2 );
string_object* string_ptr_concat( string_object* str1, string_object* str2 );
void string_append_string(string_object* str1, string_object* str2 );
void string_append_cstring(string_object* str1, const char* cstr );
string_object* string_repeat( string_object* str, int rep );
string_object* string_subset( string_object* str, size_t from_idx, size_t to_idx, const char* encoding);
int string_has_char( string_object* str, char c); /* Exist: 1, Non-Exist: 0 */
int string_str2int( string_object* str );
double string_str2double( string_object* str );
int string_compare( string_object* str1, string_object* str2);
int string_move_ptr( string_object** str1, string_object** str2);
int string_free( string_object* str );


#endif /* COMMON_STRING_H */
