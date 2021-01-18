#ifndef CPP_STRING_HPP
#define CPP_STRING_HPP

typedef void cpp_object;

#ifdef __cplusplus
extern "C" {
#endif
	cpp_object* cpp_string_new (const char*);
	cpp_object* cpp_string_new_with_len (const char* , int);
	cpp_object* cpp_string_new_unescaped_string(cpp_object* ori_str, const char* encoding);
	cpp_object* cpp_string_clone( cpp_object* );
	cpp_object* cpp_string_new_int2str( int num);
	cpp_object* cpp_string_new_double2str( double num);
//	Deprecaetd: it is confusing and difficult to track and free newly created c strings.
//	const char* cpp_string_int2cstr(int num);
//	const char* cpp_string_double2cstr(double num);
	cpp_object* cpp_string_lstrip( cpp_object* );
	cpp_object* cpp_string_rstrip( cpp_object* );
	cpp_object* cpp_string_strip( cpp_object* );
	const char* cpp_string_read (cpp_object* );
	cpp_object* cpp_string_concat (cpp_object*, cpp_object*);
	cpp_object* cpp_string_ptr_concat (cpp_object*, cpp_object*);
	void cpp_string_append_string (cpp_object*, cpp_object*);
	void cpp_string_append_cstring (cpp_object*, const char*);
	cpp_object* cpp_string_repeat (cpp_object*, int);
	cpp_object* cpp_string_subset (cpp_object*, size_t, size_t , const char* encoding);
	int cpp_string_has_char (cpp_object*, char);
	int cpp_string_str2int(cpp_object* obj);
	double cpp_string_str2double(cpp_object* obj);
	// int cpp_string_copy_ptr ( cpp_object** , cpp_object** );
	int cpp_string_compare ( cpp_object* , cpp_object* );
	int cpp_string_move_ptr ( cpp_object** , cpp_object** );
	int cpp_string_free ( cpp_object* );
#ifdef __cplusplus
}
#endif


#endif /* CPP_STRING_HPP */

