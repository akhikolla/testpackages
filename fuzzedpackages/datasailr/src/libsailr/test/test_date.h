#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

void test_date_test1( void );

void
test_date_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test dates ", test_date_test1 );
}

void
test_date_test1( void )
{
	// Code 
	const char* code = " "
"epoch_date = date_ymd(1970,1,1)\n"
"today = date_ymd(2019, 5, 1);"
"birthday = date_ymd(1999, 10, 31);"
"age = (today - birthday) / 365.25;"
"future = date_add_n_months( date_add_n_years( today, 10 ), 5);"
"future_str = date_format(future, '%Y-%m-%d');"
;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "epoch_date" );
	sailr_ptr_table_create_null(&table, "today" );
	sailr_ptr_table_create_null(&table, "birthday" );
	sailr_ptr_table_create_null(&table, "age" );
	sailr_ptr_table_create_null(&table, "future" );
	sailr_ptr_table_create_null(&table, "future_str" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_epoch_date = sailr_ptr_table_get_type(&table, "epoch_date");
	char st_today = sailr_ptr_table_get_type(&table, "today");
	char st_birthday = sailr_ptr_table_get_type(&table, "birthday");
	char st_age = sailr_ptr_table_get_type(&table, "age"); 
	char st_future = sailr_ptr_table_get_type(&table, "future");
	char st_future_str = sailr_ptr_table_get_type(&table, "future_str");

	CU_ASSERT_EQUAL( st_epoch_date , 'i');
	CU_ASSERT_EQUAL( st_today , 'i');
	CU_ASSERT_EQUAL( st_birthday , 'i');
	CU_ASSERT_EQUAL( st_age , 'd');
	CU_ASSERT_EQUAL( st_future , 'i');
	CU_ASSERT_EQUAL( st_future_str , 's');

	int s_epoch_date = *((int*) *sailr_ptr_table_get_pptr(&table, "epoch_date"));
//	int s_today = *((int*) *sailr_ptr_table_get_pptr(&table, "today"));
//	int s_birthday = *((int*) *sailr_ptr_table_get_pptr(&table, "birthday"));
	double s_age = *((double*) *sailr_ptr_table_get_pptr(&table, "age"));
//	int s_future = *((int*) *sailr_ptr_table_get_pptr(&table, "future"));
	const char* s_future_str = sailr_ptr_table_read_string(&table, "future_str");

	CU_ASSERT_EQUAL( s_epoch_date , 0 );
	CU_ASSERT_DOUBLE_EQUAL( s_age , 19.4, 0.1 );
	CU_ASSERT_STRING_EQUAL( s_future_str , "2029-10-01");

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


