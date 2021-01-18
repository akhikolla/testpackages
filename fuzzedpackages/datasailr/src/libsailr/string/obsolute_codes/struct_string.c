#include <R_ext/Print.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include "struct_string.h"

struct_string*
new_struct_string( char* str , int byte_len )
{
	size_t basic_struct_size = sizeof(struct_string);
	size_t additional_size = sizeof(char)*(byte_len + 1 );
    #ifdef DEBUG
	Rprintf("basic_struct_size: %zu \n", basic_struct_size);
	Rprintf("additional_size: %zu \n", additional_size);
    #endif
    struct_string* temp = (struct_string *)malloc( basic_struct_size + additional_size ); 
    temp->len = byte_len; 
    strcpy( temp->buf, str);
	return temp;
}

struct_string*
cat_struct_2strings( struct_string* str1, struct_string* str2)
{
	size_t basic_struct_size = sizeof(struct_string);
	int str1_size = str1->len;
	int str2_size = str2->len;
	size_t additional_size = sizeof(char)*(str1_size + str2_size + 1 );
    struct_string* temp = (struct_string *)malloc( basic_struct_size + additional_size );
	temp->len =  str1_size + str2_size;
	strcpy(temp->buf, str1->buf);
	strcat(temp->buf, str2->buf);
	return temp;
}

struct_string*
cat_struct_strings( int num_of_strings , ... )  
{
    struct_string* tmp_struct;

	// Start va_ macroes
	va_list varg_list;

	va_start(varg_list, num_of_strings );
	int i;
	int total_size = 0;
	for(i = 0 ; i < num_of_strings; i++ ){
		tmp_struct = va_arg( varg_list , struct_string* );
		total_size = total_size + (tmp_struct->len) ;
	}
    va_end(varg_list);


	size_t basic_struct_size = sizeof(struct_string);
	size_t additional_size = sizeof(char)*( total_size + 1 );
    struct_string* temp = (struct_string *)malloc( basic_struct_size + total_size );
	temp->buf[0] = '\0';
	temp->len = 0;

	va_start(varg_list, num_of_strings );
	int j;
	for(j = 0 ; j < num_of_strings; j++ ){
		tmp_struct = va_arg( varg_list , struct_string* );
		strcat(temp->buf, tmp_struct->buf);
		temp->len = temp->len + tmp_struct->len ; 
	}
    va_end(varg_list);

	return temp;
}

void
free_struct_string(struct_string* sstr)
{
	free(sstr);
}

/* 
int main(int argc, char** argv){
	char* example_str1 = "Hello World!";
	size_t example_len1 = strlen("Hello World!");
	Rprintf("%s \n", example_str1);
	Rprintf("%zu \n", example_len1);

	struct_string* str1 = new_struct_string(example_str1, example_len1 );
	Rprintf("%s \n", str1->buf);
	Rprintf("%zu \n", strlen(str1->buf));

    // No! Do not do this!!
    // https://stackoverflow.com/questions/9504588/should-i-free-char-initialized-using-string-literals
	// printf("Free original string. \n");
	// free(example_str);


	char* example_str2 = "KONNICHIWA";
	size_t example_len2 = strlen(example_str2);
	Rprintf("%s \n", example_str2);
	Rprintf("%zu \n", example_len2);

	struct_string* str2 = new_struct_string(example_str2, example_len2 );	
	Rprintf("%s \n", str2->buf);
	Rprintf("%zu \n", strlen(str2->buf));

	struct_string* strX = cat_struct_2strings( str1, str2 );	
	Rprintf("%s \n", strX->buf);
	Rprintf("%zu \n", strlen(strX->buf));
	Rprintf("%zu \n", strX->len);


	char* example_str3 = "OHAYO";
	size_t example_len3 = strlen(example_str3);
	Rprintf("%s \n", example_str3);
	Rprintf("%zu \n", example_len3);

	struct_string* str3 = new_struct_string(example_str3, example_len3 );	
	Rprintf("%s \n", str3->buf);
	Rprintf("%zu \n", strlen(str3->buf));

	struct_string* strY = cat_struct_strings( 2, str1, str3 );	
	Rprintf("%s \n", strY->buf);
	Rprintf("%zu \n", strlen(strY->buf));
	Rprintf("%zu \n", strY->len);

}
*/

