#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_string_test1( void );

void
test_string_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test strings", test_string_test1 );
}

void
test_string_test1( void )
{
	// Code 
	const char* code = " "
"hello = 'Hello' \n"
"firstname = 'Mickey' \n"
"lastname = 'Mouse' \n"
"fullname = firstname + ' ' + lastname \n"
"greeting = hello + ', ' + fullname \n"
"greeting2 = greeting \n"
"greeting = 'Good bye, ' + fullname \n" 
;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "hello" );
	sailr_ptr_table_create_null(&table, "firstname" );
	sailr_ptr_table_create_null(&table, "lastname" );
	sailr_ptr_table_create_null(&table, "fullname" );
	sailr_ptr_table_create_null(&table, "greeting" );
	sailr_ptr_table_create_null(&table, "greeting2" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_hello = sailr_ptr_table_get_type(&table, "hello");
	char st_firstname = sailr_ptr_table_get_type(&table, "firstname");
	char st_lastname = sailr_ptr_table_get_type(&table, "lastname"); 
	char st_fullname = sailr_ptr_table_get_type(&table, "fullname");
	char st_greeting = sailr_ptr_table_get_type(&table, "greeting");
	char st_greeting2 = sailr_ptr_table_get_type(&table, "greeting2");

	CU_ASSERT_EQUAL( st_hello , 's');
	CU_ASSERT_EQUAL( st_firstname , 's');
	CU_ASSERT_EQUAL( st_lastname , 's');
	CU_ASSERT_EQUAL( st_fullname , 's');
	CU_ASSERT_EQUAL( st_greeting, 's');
	CU_ASSERT_EQUAL( st_greeting2, 's');

	const char* s_hello = sailr_ptr_table_read_string(&table, "hello");
	const char* s_firstname = sailr_ptr_table_read_string(&table, "firstname");
	const char* s_lastname =  sailr_ptr_table_read_string(&table, "lastname"); 
	const char* s_fullname =  sailr_ptr_table_read_string(&table, "fullname");
	const char* s_greeting =  sailr_ptr_table_read_string(&table, "greeting");
	const char* s_greeting2 = sailr_ptr_table_read_string(&table, "greeting2");

	CU_ASSERT_STRING_EQUAL( s_hello , "Hello");
	CU_ASSERT_STRING_EQUAL( s_firstname , "Mickey");
	CU_ASSERT_STRING_EQUAL( s_lastname , "Mouse");
	CU_ASSERT_STRING_EQUAL( s_fullname , "Mickey Mouse");
	CU_ASSERT_STRING_EQUAL( s_greeting, "Good bye, Mickey Mouse");
	CU_ASSERT_STRING_EQUAL( s_greeting2, "Hello, Mickey Mouse");

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


