#include <iostream>
#include <string>
#include "cpp_string_utf8.hpp"

/* private */
int
utf8_length(std::string* str)
{
	return utf8::distance(str->begin(), str->end());
}

int
utf8_valid(std::string* str)
{
	auto iter_from = str->begin();
	auto iter_to = str->end();
	return utf8::is_valid(iter_from, iter_to);
}

std::string*
new_string_with_same_capacity_as( std::string* ori_str )
{
	std::string* new_str = new std::string();
//	std::cout << (*ori_str) << "(" << ori_str->length() << ")" << std::endl;

	if (ori_str->empty()){
		std::cout << "LENGTH is zero" << std::endl;
		return new_str;
	}

	size_t ori_capacity = ori_str->capacity() ;
	if(new_str->capacity() < (ori_capacity + 1)){
		new_str->reserve(ori_capacity + 1);
	}
	return new_str;
}

/* public */

std::string*
cpp_string_subset_utf8(std::string* ori_str, int from_idx, int to_idx)
{
	int str_len = utf8_length( ori_str );

	auto curr_iter = ori_str->begin();
	auto str_end = ori_str->end();
	decltype(curr_iter) from_iter;
	decltype(curr_iter) to_iter ;

	int curr_idx = 0;
	// uint32_t curr_cp;


	// Validate indexes specified.
	if(from_idx > to_idx){
		return (new std::string()); // Not allowed
	}else if(from_idx >= str_len){
		return (new std::string()); // The index is out of range.
	}else if(to_idx >= str_len){
		to_idx = str_len - 1; // to_idx should be within range.
	}

	std::string* result_str;
	if(ori_str->empty()){
		return (new std::string());
	}else{
		result_str = new_string_with_same_capacity_as( ori_str );
	}

	while(1){
		if(curr_idx == from_idx){
			from_iter = curr_iter;// Initialization
		}

		try{
			utf8::next(curr_iter, str_end);
			curr_idx = curr_idx + 1;
		}catch(const utf8::not_enough_room& e){
			break;
		}

		if(curr_idx == to_idx + 1){
			to_iter = curr_iter;
			break;
		}
	}

	std::copy( from_iter , to_iter, std::back_inserter( *result_str ));
	return result_str;
}

std::string*
cpp_string_new_unescaped_string_utf8( std::string* ori_str )
{
	std::string* new_str;
	if(ori_str->empty()){
		return (new std::string());
	}else{
		new_str = new_string_with_same_capacity_as( ori_str );
	}
	
	auto curr_iter = ori_str->begin();
	auto str_end = ori_str->end();

	utf8::uint32_t curr_cp;

	while(1){
		try{
			curr_cp = utf8::next(curr_iter, str_end);
		}catch(const utf8::not_enough_room& e){
			break;
		}
		if(curr_cp == U'\\'){
			try{
				curr_cp = utf8::next(curr_iter, str_end);
			}catch(const utf8::not_enough_room& e){
{}//				printf("ERROR: string litereal should not end with single backslash.");
				break;
			}
			switch(curr_cp){
			case 't' : curr_cp = U'\t';
			break;
			case 'n' : curr_cp = U'\n';
			break;
			case 'r' : curr_cp = U'\r';
			break;
			case '\\' : curr_cp = U'\\';
			break; 
			case '\'' : curr_cp = U'\'';
			break;
			case '"' : curr_cp = U'"';
			break;
			case '?' : curr_cp= U'?';
			break;
			default : // curr_cp = curr_cp ;
			break;
			}
		}
		try{
			utf8::append(curr_cp, std::back_inserter( *new_str ));
		}catch(const utf8::invalid_code_point& e){
			std::cout << "invalid code point" << curr_cp << std::endl; 
			std::cout << "CODE_POINT_MAX" << utf8::internal::CODE_POINT_MAX << std::endl; 
			break;
		}
	}
	return new_str;
}

