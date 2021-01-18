#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstring>
#include "stdlib.h"
#include "cpp_string.hpp"
#include "cpp_string_utf8.hpp"
#include "cpp_string_latin1.hpp"

#ifdef __GNUG__
  #ifndef __clang__
    #if __GNUC__ < 4 || ( __GNUC__ == 4 && ( __GNUC_MINOR__ < 5 ))
      #define NO_LAMBDA
    #endif
  #endif
#endif

#ifdef __GNUG__
  #ifndef __clang__
    #if __GNUC__ < 4 || ( __GNUC__ == 4 && ( __GNUC_MINOR__ < 5 ))
      #define NOT_ENOUGH_TO_STRING_SUPPORT
    #endif
  #endif
#endif

cpp_object*
cpp_string_new (const char* str)
{
	std::string* new_str = new std::string(str);
	return (cpp_object*) new_str;
}

cpp_object*
cpp_string_new_with_len(const char* str, int len)
{
	std::string* new_str = new std::string(str, len);
	return (cpp_object*) new_str;
}

cpp_object*
cpp_string_new_unescaped_string( cpp_object* obj , const char* encoding)
{
	std::string* ori_str = static_cast<std::string*>(obj);
	std::string* new_str;
	if( std::strcmp( encoding , "UTF8") == 0){
		new_str = cpp_string_new_unescaped_string_utf8( ori_str );
	}else if( std::strcmp( encoding , "LATIN1") == 0){
		new_str = cpp_string_new_unescaped_string_latin1( ori_str );
	}else{
		new_str = cpp_string_new_unescaped_string_utf8( ori_str ); // Default: UTF8
	}
	return (cpp_object*) new_str;
}


cpp_object*
cpp_string_clone (cpp_object* obj)
{
	std::string* ori_str = static_cast<std::string*>(obj);
	std::string* new_str = new std::string(*ori_str);
	return (cpp_object*) new_str;
}

cpp_object*
cpp_string_new_int2str(int num)
{
	std::string* p_str;
//	std::stringstream ss;
//	ss.clear();
//	ss << num ;
//	p_str = new std::string(ss.str());

#ifdef NOT_ENOUGH_TO_STRING_SUPPORT
	p_str = new std::string( std::to_string( (long long int) num));
#else
	p_str = new std::string( std::to_string( num ));
#endif

	return (cpp_object*) p_str;
}

cpp_object*
cpp_string_new_double2str(double num)
{
	std::string* p_str;
//	std::stringstream ss;
//	ss.clear();
//	ss << num ;
//	p_str = new std::string(ss.str());

#ifdef NOT_ENOUGH_TO_STRING_SUPPORT
	p_str = new std::string( std::to_string( (long double) num));
#else
	p_str = new std::string( std::to_string( num ));
#endif
	return (cpp_object*) p_str;
}

/* Deprecaetd: it is confusing and difficult to track and free newly created c strings.
const char* cpp_string_int2cstr(int num){}
const char* cpp_string_double2cstr(double num){}
*/

cpp_object*
cpp_string_lstrip(cpp_object* obj)
{
	std::string* str = static_cast<std::string*>(obj);
	std::string* new_str = new std::string(*str);

#ifdef NO_LAMBDA
    new_str->erase(new_str->begin(), std::find_if(new_str->begin(), new_str->end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
#else
	new_str->erase(new_str->begin(), std::find_if_not(new_str->begin(), new_str->end(), [](int c){return std::isspace(c);}));
#endif

    return (cpp_object*) new_str;
}

cpp_object*
cpp_string_rstrip(cpp_object* obj)
{
	std::string* str = static_cast<std::string*>(obj);
	std::string* new_str = new std::string(*str);

#ifdef NO_LAMBDA
    new_str->erase(std::find_if(new_str->rbegin(), new_str->rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), new_str->end());
#else
	new_str->erase(std::find_if_not(new_str->rbegin(), new_str->rend(), [](int c){return std::isspace(c);}).base(), new_str->end());
#endif

    return (cpp_object*) new_str;
}

cpp_object*
cpp_string_strip(cpp_object* obj)
{
	std::string* str = static_cast<std::string*>(obj);
	std::string* lstripped_str = static_cast<std::string*>(cpp_string_lstrip(str));
    std::string* sripped_str = static_cast<std::string*>(cpp_string_rstrip(lstripped_str));
	delete(lstripped_str);
    return (cpp_object*) sripped_str;
}


const char*
cpp_string_read (cpp_object* obj )
{
	std::string* str;
	if (obj != NULL){
		str = static_cast<std::string*>(obj);
//		std::cout << *str << std::endl;
//		printf("%s\n", str->c_str());
		return str->c_str();
	}else{ 
		std::cout << "WARNING: String is NULL?\n" << std::endl;
		return "";
	}
}

cpp_object*
cpp_string_concat (cpp_object* obj1 , cpp_object* obj2 )
{
//	std::cout << "start concat" << std::endl;
	std::string* str1 = static_cast<std::string*>(obj1);
	std::string* str2 = static_cast<std::string*>(obj2);
	std::stringstream ss;
	ss << *str1 << *str2;
	std::string* new_p_str = new std::string( ss.str() );
	return (void*) new_p_str;
}


cpp_object*
cpp_string_ptr_concat (cpp_object* obj1 , cpp_object* obj2 )
{
	std::string* str1 = static_cast<std::string*>(obj1);
	// std::cout << (*str1) << std::endl;
	std::string* str2 = static_cast<std::string*>(obj2);
	// std::cout << (*str2) << std::endl;
	std::stringstream ss;
	ss << *str1 << *str2;
	std::string* new_p_str = new std::string( ss.str() );
	return (void*) new_p_str;
}

void
cpp_string_append_string(cpp_object* obj1 , cpp_object* obj2)
{
	std::string* str1 = static_cast<std::string*>(obj1);
	std::string* str2 = static_cast<std::string*>(obj2);
	str1->append(*str2);
}

void
cpp_string_append_cstring(cpp_object* obj1, const char* cstr)
{
	std::string* str1 = static_cast<std::string*>(obj1);
	str1->append(cstr);
}

cpp_object*
cpp_string_repeat(cpp_object* obj, int rep)
{
	std::string* cpp_str = static_cast<std::string*>(obj);

	std::stringstream ss;
	if (rep <= 0 ){
		std::cout << "ERROR: rep should be greater than 0. \n";
		exit(0);
	}
		
	for( ; rep > 0; rep-- ){
		ss << (*cpp_str) ;
	}
	std::string* new_str = new std::string( ss.str() );
	return (void*) new_str;	
}

cpp_object*
cpp_string_subset (cpp_object* obj, size_t from_idx , size_t to_idx , const char* encoding )  // index starts from zero.
{
	std::string* ori_str = static_cast<std::string*>(obj);
	std::string* new_str;
	if( std::strcmp( encoding , "UTF8") == 0){
		new_str = cpp_string_subset_utf8( ori_str, from_idx, to_idx );
	}else if( std::strcmp( encoding , "LATIN1") == 0){
		new_str = cpp_string_subset_latin1( ori_str, from_idx, to_idx  );
	}else{
		new_str = cpp_string_subset_utf8( ori_str, from_idx, to_idx ); // Default: UTF8
	}
	return (cpp_object*) new_str;
}

int
cpp_string_has_char (cpp_object* obj, char c)
{
	std::string* cpp_str = static_cast<std::string*>(obj);
	std::size_t pos = cpp_str->find( c );
	if( pos != std::string::npos ){
		return 1;
	} else {
		return 0 ;
	} 
}


int
cpp_string_str2int(cpp_object* obj)
{
	int tmp_int;
	std::string* str = static_cast<std::string*>(obj);
	std::istringstream ( *str ) >> tmp_int;
	return tmp_int;
}


double
cpp_string_str2double(cpp_object* obj)
{
	double tmp_dbl;
	std::string* str = static_cast<std::string*>(obj);
	std::istringstream ( *str ) >> tmp_dbl;
	return tmp_dbl;
}


int
cpp_string_compare ( cpp_object* obj1, cpp_object* obj2)
{
	std::string* cpp_str1 = static_cast<std::string*>(obj1);
	std::string* cpp_str2 = static_cast<std::string*>(obj2);
	if( cpp_str1->compare(*cpp_str2) == 0 ) { // matched
		return 1;
	} else { // not matched
		return 0;
	}

}

int
cpp_string_copy_ptr(cpp_object** ptrptr , cpp_object** obj)  // This function is dangerous.
{
	*ptrptr = *obj;
	return 1;
}

int
cpp_string_move_ptr(cpp_object** ptrptr , cpp_object** obj)
{
	cpp_object* old_object = *ptrptr;
	cpp_string_copy_ptr(ptrptr, obj);
	delete (std::string* )old_object ;
	*obj = NULL;
	return 1;
}

int
cpp_string_free ( cpp_object* obj)
{
	std::string* cpp_str = static_cast<std::string*>(obj);
	delete cpp_str ;
	return 1;
}



