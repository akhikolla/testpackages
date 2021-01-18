#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_comment_test1( void );

void
test_comment_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test ifs", test_comment_test1 );
}

void
test_comment_test1( void )
{
	// Code 
	const char* code = " "
"old_asian = 0 \n"
"age = 65 \n"
"/* Old Asians live \n"
"all over the world. \n"
"*/ \n"
"race = 'Asian-Pac-Islander' \n"
" \n"
"// Hello \n"
"if(race == 'White'){ \n"
"  race_id = 1   // White \n"
"}else if(race == 'Black'){ \n"
"  race_id = 2    // Black \n"
"}else if(race =~ re/Asian/ ){ \n"
"  race_id = 3 ; // Asian \n"
"  if(age > 50 ){   \n"
"    old_asian = 1  // Old Asian \n"
"  } \n"
"}else{ \n"
"  race_id = 4  // Others \n"
"} \n"
"";

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "old_asian" );
	sailr_ptr_table_create_null(&table, "age" );
	sailr_ptr_table_create_null(&table, "race" );
	sailr_ptr_table_create_null(&table, "race_id" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_old_asian = sailr_ptr_table_get_type(&table, "old_asian");
	char st_age = sailr_ptr_table_get_type(&table, "age");
	char st_race = sailr_ptr_table_get_type(&table, "race");
	char st_race_id = sailr_ptr_table_get_type(&table, "race_id");

	CU_ASSERT_EQUAL( st_old_asian , 'i');
	CU_ASSERT_EQUAL( st_age , 'i');
	CU_ASSERT_EQUAL( st_race , 's');
	CU_ASSERT_EQUAL( st_race_id , 'i');

	int s_old_asian = *((int*) *sailr_ptr_table_get_pptr(&table, "old_asian"));
	int s_age = *((int*) *sailr_ptr_table_get_pptr(&table, "age"));
	const char* s_race = sailr_ptr_table_read_string(&table, "race"); 
	int s_race_id = *((int*) *sailr_ptr_table_get_pptr(&table, "race_id")); 

	CU_ASSERT_EQUAL( s_old_asian , 1 );
	CU_ASSERT_EQUAL( s_age , 65 );
	CU_ASSERT_STRING_EQUAL( s_race , "Asian-Pac-Islander");
	CU_ASSERT_EQUAL( s_race_id , 3 );

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


