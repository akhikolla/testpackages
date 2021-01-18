#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_rexgexp_test1( void );

void
test_regexp_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test rexgexp ", test_rexgexp_test1 );
}

void
test_rexgexp_test1( void )
{
	// Code 
	const char* code = " "
"greeting = 'Hello World' \n"
"greeting_rexp = re/^(\\w*)/ \n"
"greeting_rexp =~ greeting  \n"
"hello = rexp_matched(1) \n"
"\n"
"email_tom = 'tom2019@gmail.com' \n"
"email_rexp = re/^([a-zA-z0-9\\.\\+\\-\\*\\/\\?]+)@([a-zA-Z0-9\\.]+)\\.([a-zA-Z]+?)$/ \n"
"email_rexp =~ email_tom \n"
"email_full = rexp_matched(0) \n"
"email_local =  rexp_matched(1) \n"
"email_domain = rexp_matched(2) + '.' + rexp_matched(3) \n"
"\n"
;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_parser_state_set_source_encoding( ps , "UTF-8"); // Optional
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "greeting" );
	sailr_ptr_table_create_null(&table, "greeting_rexp" );
	sailr_ptr_table_create_null(&table, "hello" );
	sailr_ptr_table_create_null(&table, "email_tom" );
	sailr_ptr_table_create_null(&table, "email_rexp" );
	sailr_ptr_table_create_null(&table, "email_full" );
	sailr_ptr_table_create_null(&table, "email_local" );
	sailr_ptr_table_create_null(&table, "email_domain" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_greeting = sailr_ptr_table_get_type(&table, "greeting");
	char st_greeting_rexp = sailr_ptr_table_get_type(&table, "greeting_rexp");
	char st_hello = sailr_ptr_table_get_type(&table, "hello"); 
	char st_email_tom = sailr_ptr_table_get_type(&table, "email_tom");
	char st_email_rexp = sailr_ptr_table_get_type(&table, "email_rexp");
	char st_email_full = sailr_ptr_table_get_type(&table, "email_full");
	char st_email_local = sailr_ptr_table_get_type(&table, "email_local");
	char st_email_domain = sailr_ptr_table_get_type(&table, "email_domain");

	CU_ASSERT_EQUAL( st_greeting , 's');
	CU_ASSERT_EQUAL( st_greeting_rexp , 'r');
	CU_ASSERT_EQUAL( st_hello , 's');
	CU_ASSERT_EQUAL( st_email_tom , 's');
	CU_ASSERT_EQUAL( st_email_rexp, 'r');
	CU_ASSERT_EQUAL( st_email_full, 's');
	CU_ASSERT_EQUAL( st_email_local, 's');
	CU_ASSERT_EQUAL( st_email_domain, 's');

	const char* s_greeting = sailr_ptr_table_read_string(&table, "greeting");
//	const char* s_greeting_rexp = sailr_ptr_table_read_string(&table, "greeting_rexp");
	const char* s_hello =  sailr_ptr_table_read_string(&table, "hello"); 
	const char* s_email_tom =  sailr_ptr_table_read_string(&table, "email_tom");
//	const char* s_email_rexp =  sailr_ptr_table_read_string(&table, "email_rexp");
	const char* s_email_full = sailr_ptr_table_read_string(&table, "email_full");
	const char* s_email_local = sailr_ptr_table_read_string(&table, "email_local");
	const char* s_email_domain = sailr_ptr_table_read_string(&table, "email_domain");

	CU_ASSERT_STRING_EQUAL( s_greeting , "Hello World");
	CU_ASSERT_STRING_EQUAL( s_hello , "Hello");
	CU_ASSERT_STRING_EQUAL( s_email_tom , "tom2019@gmail.com");
	CU_ASSERT_STRING_EQUAL( s_email_full , "tom2019@gmail.com");
	CU_ASSERT_STRING_EQUAL( s_email_local , "tom2019");
	CU_ASSERT_STRING_EQUAL( s_email_domain , "gmail.com");

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


