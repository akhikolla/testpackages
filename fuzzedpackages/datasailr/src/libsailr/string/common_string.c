#include <R_ext/Print.h>
#include "common_string.h"
#include <stdio.h>
#include <stdint.h>

string_object*
string_new(const c_char* str )
{
	return cpp_string_new( str );
}

string_object*
string_new_with_len(const c_char* str, int len)
{
	return cpp_string_new_with_len( str, len);
}

string_object*
string_new_unescaped_and_delete_ori(string_object* str, const char* encoding)
{
	string_object* new_str = cpp_string_new_unescaped_string( str , encoding);
	string_free( str );
	return new_str;
}


string_object*
string_clone(string_object* str)
{
	return (string_object*) cpp_string_clone( (cpp_object*) str );
}

string_object*
string_new_int2str(int num)
{
	return cpp_string_new_int2str(num);
}

string_object*
string_new_double2str(double num)
{
	return cpp_string_new_double2str(num);
}

/* Deprecaetd: it is confusing and difficult to track and free newly created c strings.
const char* string_int2cstr(int num){}
const char* string_double2cstr(double num){}
*/

string_object*
string_strip( string_object* str)
{
	return (string_object*) cpp_string_strip((cpp_object*) str);
}

string_object*
string_lstrip( string_object* str)
{
	return (string_object*) cpp_string_lstrip((cpp_object*) str);
}

string_object*
string_rstrip( string_object* str)
{
	return (string_object*) cpp_string_rstrip((cpp_object*) str);
}

const c_char*
string_read( string_object* str)
{
  const char* c_str = cpp_string_read( str );
//  printf("C STRING: %s", c_str);
	return c_str;
}

string_object*
string_concat( string_object* str1, string_object* str2 )
{
//	printf( "%s + %s \n", string_read(str1), string_read(str2) );
	return (string_object*) cpp_string_concat((cpp_object*) str1, (cpp_object*) str2 );
}


string_object*
string_ptr_concat( string_object* str1, string_object* str2 )
{
	string_object* str_obj = (string_object**) cpp_string_ptr_concat( str1, str2 );
	return str_obj;
}

void
string_append_string(string_object* str1, string_object* str2)
{
	cpp_string_append_string( (cpp_object*)str1, (cpp_object*)str2 );
}

void
string_append_cstring(string_object* str1, const char* cstr)
{
	cpp_string_append_cstring( (cpp_object*)str1, cstr );
}

string_object*
string_repeat( string_object* str, int rep )
{
	return cpp_string_repeat( str, rep );
}

string_object*
string_subset( string_object* str, size_t from_idx, size_t to_idx, const char* encoding)
{
	return cpp_string_subset( str, from_idx, to_idx, encoding );
}


int
string_has_char( string_object* str, char c)
{
	return cpp_string_has_char((cpp_object*) str, c);
}

int
string_str2int( string_object* str )
{
	return cpp_string_str2int((cpp_object*) str );
}

double
string_str2double( string_object* str )
{
	return cpp_string_str2double((cpp_object*) str );
}

int
string_compare( string_object* str1, string_object* str2)
{
	return cpp_string_compare( str1, str2);
}

int
string_move_ptr( string_object** str1, string_object** str2)
{
	return cpp_string_move_ptr(str1, str2);
}

int
string_free( string_object* str )
{
	return cpp_string_free( str);
}


