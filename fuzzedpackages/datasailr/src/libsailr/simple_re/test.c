#include <R_ext/Print.h>
#include "simple_re.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv){
	
	simple_re* re = simple_re_compile("^(\\S)(.)+",  "UTF-8");

	Rprintf("Compilation finished. \n");	

//	const char* string = "こんにちは World";
	const char* string = "Hello World";
	simple_re* ptr_to_store;
	int result = simple_re_match( re , string , &ptr_to_store);
	Rprintf("Matching finished. \n");

	int num;
	char* str;
	int idx;
	if(result == 1 ){
		Rprintf("success match found. \n");
		Rprintf("Showing ... \n");
		
		num = simple_re_matched_group_num(ptr_to_store);
		for(idx = 0; idx < num ; ++idx){
			str = simple_re_matched_str( ptr_to_store , idx );
			if(idx == 0){
				Rprintf("Whole Matched : %s \n", str );
			} else {
				Rprintf("Matched Group %d : %s \n", idx, str );
			}
			free(str);
		}
	}else{
		Rprintf("No match. %d \n", result);
	}
	simple_re_free(re);
}
