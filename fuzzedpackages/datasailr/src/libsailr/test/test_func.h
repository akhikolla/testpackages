#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_func_test1( void );

void
test_func_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test funcs ", test_func_test1 );
}

void
test_func_test1( void )
{
	// Code 
	const char* code = " "
"space_pi = ' ' + str_subset(num_to_str(3.14),1,4) + num_to_str(1592) + num_to_str(6535) ;"
"pi_str = 'PI is' + space_pi ;"
"pi_num = str_to_num( str_strip(space_pi)) ;"
"three_str = str_subset( num_to_str(10/3) , 1 , 4);"
"seven_str = str_repeat( num_to_str(7), 7);"
"ten_str = num_to_str( 2.5 * 4 ) \n"
"\n"
"spaced_str = '    hello world    ';"
"stripped_str = str_strip(spaced_str);"
"jpn_hello = 'こんにちは';"
"new_str = str_concat( jpn_hello, ' ' , stripped_str );"
"new_str2 = str_concat( new_str, ', PI = ', pi_num);"
;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "space_pi" );
	sailr_ptr_table_create_null(&table, "pi_str" );
	sailr_ptr_table_create_null(&table, "pi_num" );
	sailr_ptr_table_create_null(&table, "three_str" );
	sailr_ptr_table_create_null(&table, "seven_str" );
	sailr_ptr_table_create_null(&table, "ten_str" );

	sailr_ptr_table_create_null(&table, "spaced_str");
	sailr_ptr_table_create_null(&table, "stripped_str");

	sailr_ptr_table_create_null(&table, "jpn_hello");
	sailr_ptr_table_create_null(&table, "new_str");
	sailr_ptr_table_create_null(&table, "new_str2");

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_space_pi = sailr_ptr_table_get_type(&table, "space_pi");
	char st_pi_str = sailr_ptr_table_get_type(&table, "pi_str");
	char st_pi_num = sailr_ptr_table_get_type(&table, "pi_num"); 
	char st_three_str = sailr_ptr_table_get_type(&table, "three_str");
	char st_seven_str = sailr_ptr_table_get_type(&table, "seven_str");
	char st_ten_str = sailr_ptr_table_get_type(&table, "ten_str");
	char st_spaced_str = sailr_ptr_table_get_type(&table, "spaced_str");
	char st_stripped_str = sailr_ptr_table_get_type(&table, "stripped_str");
	char st_new_str = sailr_ptr_table_get_type(&table, "new_str");
	char st_new_str2 = sailr_ptr_table_get_type(&table, "new_str2");

	CU_ASSERT_EQUAL( st_space_pi , 's');
	CU_ASSERT_EQUAL( st_pi_str , 's');
	CU_ASSERT_EQUAL( st_pi_num , 'd');
	CU_ASSERT_EQUAL( st_three_str , 's');
	CU_ASSERT_EQUAL( st_seven_str , 's');
	CU_ASSERT_EQUAL( st_ten_str , 's');
	CU_ASSERT_EQUAL( st_spaced_str , 's');
	CU_ASSERT_EQUAL( st_stripped_str , 's');
	CU_ASSERT_EQUAL( st_new_str , 's');
	CU_ASSERT_EQUAL( st_new_str2 , 's');


	const char* s_space_pi = sailr_ptr_table_read_string(&table, "space_pi");
	const char* s_pi_str = sailr_ptr_table_read_string(&table, "pi_str");
	double s_pi_num = *((double*) *sailr_ptr_table_get_pptr(&table, "pi_num"));
	const char* s_three_str = sailr_ptr_table_read_string(&table, "three_str");
	const char* s_seven_str = sailr_ptr_table_read_string(&table, "seven_str");
	const char* s_ten_str = sailr_ptr_table_read_string(&table, "ten_str");
	const char* s_spaced_str = sailr_ptr_table_read_string(&table, "spaced_str");
	const char* s_stripped_str = sailr_ptr_table_read_string(&table, "stripped_str");
	const char* s_new_str = sailr_ptr_table_read_string(&table, "new_str");
	const char* s_new_str2 = sailr_ptr_table_read_string(&table, "new_str2");

	CU_ASSERT_STRING_EQUAL( s_space_pi , " 3.1415926535");
	CU_ASSERT_STRING_EQUAL( s_pi_str , "PI is 3.1415926535");
	CU_ASSERT_DOUBLE_EQUAL( s_pi_num , 3.1415926535, 0.00001 );
	CU_ASSERT_STRING_EQUAL( s_three_str , "3.33");
	CU_ASSERT_STRING_EQUAL( s_seven_str , "7777777");
	CU_ASSERT_STRING_EQUAL( s_ten_str , "10.000000");
	CU_ASSERT_STRING_EQUAL( s_spaced_str , "    hello world    ");
	CU_ASSERT_STRING_EQUAL( s_stripped_str , "hello world");
	CU_ASSERT_STRING_EQUAL( s_new_str , "こんにちは hello world");
	CU_ASSERT_STRING_EQUAL( s_new_str2 , "こんにちは hello world, PI = 3.141593" );

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


