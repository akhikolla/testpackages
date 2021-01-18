#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_encoding_test1( void );

void
test_encoding_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test encoding", test_encoding_test1 );
}

void
test_encoding_test1( void )
{
	// Code 
	const char* code = " "
"hello_world_en = 'Hello World' \n"
"hello_world_jp = 'こんにちは世界' \n"
"hello_en=str_subset(hello_world_en, 1, 5) \n"
"hello_jp=str_subset(hello_world_jp, 1, 5) \n"
;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "hello_world_en" );
	sailr_ptr_table_create_null(&table, "hello_world_jp" );
	sailr_ptr_table_create_null(&table, "hello_en" );
	sailr_ptr_table_create_null(&table, "hello_jp" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	const char* ori_encoding = sailr_vm_stack_get_encoding(vmstack);
	sailr_vm_stack_set_encoding(vmstack, "LATIN1");
	const char* latin1_encoding = sailr_vm_stack_get_encoding(vmstack);
	sailr_vm_stack_set_encoding(vmstack, "UTF8");

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	CU_ASSERT_STRING_EQUAL( ori_encoding , "UTF8");
	CU_ASSERT_STRING_EQUAL( latin1_encoding , "LATIN1");

	char st_hello_world_en = sailr_ptr_table_get_type(&table, "hello_world_en");
	char st_hello_world_jp = sailr_ptr_table_get_type(&table, "hello_world_jp");
	char st_hello_en = sailr_ptr_table_get_type(&table, "hello_en"); 
	char st_hello_jp = sailr_ptr_table_get_type(&table, "hello_jp");

	CU_ASSERT_EQUAL( st_hello_world_en , 's');
	CU_ASSERT_EQUAL( st_hello_world_jp , 's');
	CU_ASSERT_EQUAL( st_hello_en , 's');
	CU_ASSERT_EQUAL( st_hello_jp , 's');

	const char* s_hello_world_en = sailr_ptr_table_read_string(&table, "hello_world_en");
	const char* s_hello_world_jp = sailr_ptr_table_read_string(&table, "hello_world_jp");
	const char* s_hello_en =  sailr_ptr_table_read_string(&table, "hello_en"); 
	const char* s_hello_jp =  sailr_ptr_table_read_string(&table, "hello_jp");

	CU_ASSERT_STRING_EQUAL( s_hello_world_en , "Hello World");
	CU_ASSERT_STRING_EQUAL( s_hello_world_jp , "こんにちは世界");
	CU_ASSERT_STRING_EQUAL( s_hello_en , "Hello");
	CU_ASSERT_STRING_EQUAL( s_hello_jp , "こんにちは");

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


