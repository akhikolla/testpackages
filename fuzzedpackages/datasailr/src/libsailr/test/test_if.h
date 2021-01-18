#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_if_test1( void );

void
test_if_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test ifs", test_if_test1 );
}

void
test_if_test1( void )
{
	// Code 
	const char* code = " "
"greeting1 = '' \n"
"message1 = ''\n"
"time = 15 \n"
"country = 'Germany' \n"
"if( country == 'Germany' ){\n"
"  if( time < 12 ){ \n"
"    greeting1 = 'Guten Morgen' "
"  }else if( 18 < time ){ \n"
"    greeting1 = 'Guten Abent' "
"  }else if( 12 <= time &&  time <= 18 ){ \n"
"    greeting1 = 'Hallo' \n"
"  }\n"
"  if( time == 15)\n"
"    message1 = 'It is snack time!'\n"
"}\n" ;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "greeting1" );
	sailr_ptr_table_create_null(&table, "time" );
	sailr_ptr_table_create_null(&table, "country" );
	sailr_ptr_table_create_null(&table, "message1" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_greeting1 = sailr_ptr_table_get_type(&table, "greeting1");
	char st_time = sailr_ptr_table_get_type(&table, "time");
	char st_country = sailr_ptr_table_get_type(&table, "country");
	char st_message1 = sailr_ptr_table_get_type(&table, "message1");

	CU_ASSERT_EQUAL( st_greeting1 , 's');
	CU_ASSERT_EQUAL( st_time , 'i');
	CU_ASSERT_EQUAL( st_country , 's');
	CU_ASSERT_EQUAL( st_message1 , 's');

	const char* s_greeting1 = sailr_ptr_table_read_string(&table, "greeting1");
	int s_time = *((int*) *sailr_ptr_table_get_pptr(&table, "time"));
	const char* s_country = sailr_ptr_table_read_string(&table, "country"); 
	const char* s_message1 = sailr_ptr_table_read_string(&table, "message1"); 

	CU_ASSERT_STRING_EQUAL( s_greeting1 , "Hallo");
	CU_ASSERT_EQUAL( s_time , 15 );
	CU_ASSERT_STRING_EQUAL( s_country , "Germany");
	CU_ASSERT_STRING_EQUAL( s_message1 , "It is snack time!");

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


