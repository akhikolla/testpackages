#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_int_test1( void );

void
test_int_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test integers", test_int_test1 );
}

void
test_int_test1( void )
{
	// Code 
	const char* code = " "
"age = 37 \n"
"bw = 90 \n"
"height = 172 \n"
"bmi = bw / (height/100) / (height/100) \n"
"obesity = 0 \n"
" \n"
"if( bmi > 30){ \n"
"    obesity = 1 \n"
"} else { \n"
"    obesity = 0 \n"
"} \n" ;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "age" );
	sailr_ptr_table_create_null(&table, "bw" );
	sailr_ptr_table_create_null(&table, "height" );
	sailr_ptr_table_create_null(&table, "bmi" );
	sailr_ptr_table_create_null(&table, "obesity" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_age = sailr_ptr_table_get_type(&table, "age");
	char st_bw = sailr_ptr_table_get_type(&table, "bw");
	char st_height = sailr_ptr_table_get_type(&table, "height"); 
	char st_bmi = sailr_ptr_table_get_type(&table, "bmi");
	char st_obesity = sailr_ptr_table_get_type(&table, "obesity");

	CU_ASSERT_EQUAL( st_age , 'i');
	CU_ASSERT_EQUAL( st_bw , 'i');
	CU_ASSERT_EQUAL( st_height , 'i');
	CU_ASSERT_EQUAL( st_bmi , 'd');
	CU_ASSERT_EQUAL( st_obesity, 'i');

	int s_age = *((int*) *sailr_ptr_table_get_pptr(&table, "age"));
	int s_bw = *((int*) *sailr_ptr_table_get_pptr(&table, "bw"));
	int s_height = *((int*) *sailr_ptr_table_get_pptr(&table, "height")); 
	double s_bmi = *((double*) *sailr_ptr_table_get_pptr(&table, "bmi"));
	int s_obesity = *((int*) *sailr_ptr_table_get_pptr(&table, "obesity"));

	CU_ASSERT_EQUAL( s_age , 37);
	CU_ASSERT_EQUAL( s_bw , 90);
	CU_ASSERT_EQUAL( s_height , 172);
	double bmi_should_be = 90 / (172/100.0) / (172/100.0);
	CU_ASSERT_DOUBLE_EQUAL( s_bmi , bmi_should_be, 0.00000001 );
	CU_ASSERT_EQUAL( s_obesity, 1);

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


