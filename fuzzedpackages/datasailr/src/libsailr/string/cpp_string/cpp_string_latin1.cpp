#include <iostream>
#include <string>
#include "cpp_string_latin1.hpp"

// latin1

std::string*
cpp_string_new_unescaped_string_latin1( std::string* ori_str )
{
	std::string* new_str = new std::string();
	std::cout << (*ori_str) << "(" << ori_str->length() << ")" << std::endl;

	if (ori_str->empty()){
		std::cout << "LENGTH is zero" << std::endl;
		return new_str;
	}

	size_t ori_capacity = ori_str->capacity() ;
	if(new_str->capacity() < (ori_capacity + 1)){
		new_str->reserve(ori_capacity + 1);
	}
	
	char new_char;
	for( auto it = ori_str->begin(); it != ori_str->end(); ++it){
		if( *it == '\\'){
			++it;
			if( it == ori_str->end()){
{}//				printf("ERROR: string litereal should not end with single backslash.");
				break;
;			}
			switch (*it){
			case 't' : new_char = '\t';
			break;
			case 'n' : new_char = '\n';
			break;
			case 'r' : new_char = '\r';
			break;
			case '\\' : new_char = '\\';
			break; 
			case '\'' : new_char = '\'';
			break;
			case '"' : new_char = '"';
			break;
			case '?' : new_char = '?';
			break;
			default : new_char = *it ;
			break;
			}
		}else{
			new_char = *it;
		}
//		std::cout << "---" << new_char << "---" << std::endl;
		new_str->push_back(new_char) ;
	}
	return new_str;
}

std::string*
cpp_string_subset_latin1 (std::string* ori_str, size_t from_idx , size_t to_idx )  // index starts from zero.
{
	if(from_idx > to_idx ){
		int temp_idx = to_idx;
		to_idx = from_idx;
		from_idx = temp_idx;
	}
	if(to_idx >= ori_str->size()){
		to_idx = ori_str->size() - 1;
	}
	std::string* new_str = new std::string( ori_str->substr(from_idx, (to_idx - from_idx + 1 )));
	return new_str;	
}
