#include <R_ext/Print.h>
#include "sailr.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <math.h>

void test_missing_test1( void );

void
test_missing_add_tests(CU_pSuite testSuite)
{
	CU_add_test(testSuite, "test missings", test_missing_test1 );
}

void
test_missing_test1( void )
{
	// Code 
	const char* code = " "
"us_pop = 3.3 \n"
"phil_pop = 1 \n"
"mars_pop = . \n"
"mars_us_ratio = mars_pop / us_pop  // nan\n"
"mars_phil_ratio = mars_pop / phil_pop // nan \n"
"us_mars_ratio = us_pop / mars_pop // inf or nan \n"
" \n"
"message1 = '' \n"
"if( mars_us_ratio > 1 ){ \n"
"  message1 = 'Mars is bigger than US.' \n"
"}else if(mars_us_ratio <= 1){\n"
"  message1 = 'Mars is not bigger than US.'  \n"
"}else{ \n"
"  if( mars_us_ratio == . ) {\n"
"    message1 = 'The relationship data is missing.' \n"
"  }else{\n"
"    message1 = 'The relation is not unknown, but not missing.'\n"
"  }\n"
"}\n"
" \n"
"message2 = '' \n"
"if( mars_phil_ratio > 1 ){ \n"
"  message2 = 'Mars is bigger than Phillipine.'  \n"
"}else if(mars_phil_ratio <= 1){\n"
"  message2 = 'Mars is not bigger than Phillipine.'  \n"
"}else{ \n"
"  if( mars_phil_ratio == . ) {\n"
"    message2 = 'The relationship data is missing.' \n"
"  }else{\n"
"    message2 = 'The relation is not unknown, but not missing.'\n"
"  }\n"
"}\n" 
" \n"
"message.us.mars.ratio = '' \n"
"if( us_mars_ratio > 1 ){ \n"
"  message.us.mars.ratio = 'US is bigger than Mars.'  \n"
"}else if( us_mars_ratio <= 1){\n"
"  message.us.mars.ratio = 'US is not bigger than Mars.'  \n"
"}else{ \n"
"  if( us_mars_ratio == . ) {\n"
"    message.us.mars.ratio = 'The relationship data is missing.' \n"
"  }else{\n"
"    message.us.mars.ratio = 'The relation is not unknown, but not missing.'\n"
"  }\n"
"}\n" ;

	// Parser Initialization
	ptr_table_object* table = sailr_ptr_table_init() ;
	parser_state_object* ps = sailr_new_parser_state ("souce from string literal", table);
	sailr_run_parser( code, ps ); 


	// Add variables
	sailr_ptr_table_create_null(&table, "us_pop" );
	sailr_ptr_table_create_null(&table, "phil_pop" );
	sailr_ptr_table_create_null(&table, "mars_pop" );
	sailr_ptr_table_create_null(&table, "mars_us_ratio" );
	sailr_ptr_table_create_null(&table, "mars_phil_ratio" );
	sailr_ptr_table_create_null(&table, "us_mars_ratio" );
	sailr_ptr_table_create_null(&table, "message1" );
	sailr_ptr_table_create_null(&table, "message2" );
	sailr_ptr_table_create_null(&table, "message.us.mars.ratio" );

	// Creating virtual machine codes
	vm_inst_object* inst_list = sailr_gen_code( ps, table); // VM Code is generated.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list);
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	vm_stack_object* vmstack = sailr_vm_stack_init();

	// Run
	sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack);

	// Assert
	// sailr_ptr_table_show_all(&table);

	char st_us_pop = sailr_ptr_table_get_type(&table, "us_pop");
	char st_phil_pop = sailr_ptr_table_get_type(&table, "phil_pop");
	char st_mars_pop = sailr_ptr_table_get_type(&table, "mars_pop"); 
	char st_mars_us_ratio = sailr_ptr_table_get_type(&table, "mars_us_ratio");
	char st_mars_phil_ratio = sailr_ptr_table_get_type(&table, "mars_phil_ratio");
	char st_us_mars_ratio = sailr_ptr_table_get_type(&table, "us_mars_ratio");
	char st_message1 = sailr_ptr_table_get_type(&table, "message1");
	char st_message2 = sailr_ptr_table_get_type(&table, "message2");
	char st_message3 = sailr_ptr_table_get_type(&table, "message.us.mars.ratio");

	CU_ASSERT_EQUAL( st_us_pop , 'd');
	CU_ASSERT_EQUAL( st_phil_pop , 'i');
	CU_ASSERT_EQUAL( st_mars_pop , 'd');
	CU_ASSERT_EQUAL( st_mars_us_ratio , 'd');
	CU_ASSERT_EQUAL( st_mars_phil_ratio, 'd');
	CU_ASSERT_EQUAL( st_us_mars_ratio, 'd');
	CU_ASSERT_EQUAL( st_message1, 's');
	CU_ASSERT_EQUAL( st_message2, 's');
	CU_ASSERT_EQUAL( st_message3, 's');

	double s_us_pop = *((double*) *sailr_ptr_table_get_pptr(&table, "us_pop"));
	int s_phil_pop = *((int*) *sailr_ptr_table_get_pptr(&table, "phil_pop"));
	double s_mars_pop = *((double*) *sailr_ptr_table_get_pptr(&table, "mars_pop")); 
	double s_mars_us_ratio = *((double*) *sailr_ptr_table_get_pptr(&table, "mars_us_ratio"));
	double s_mars_phil_ratio = *((double*) *sailr_ptr_table_get_pptr(&table, "mars_phil_ratio"));
	double s_us_mars_ratio = *((double*) *sailr_ptr_table_get_pptr(&table, "us_mars_ratio"));
	const char* s_message1 = sailr_ptr_table_read_string(&table, "message1");
	const char* s_message2 = sailr_ptr_table_read_string(&table, "message2");
	const char* s_message3 = sailr_ptr_table_read_string(&table, "message.us.mars.ratio"); 

//	double nan_value = sqrt(-1);
//	printf("nan value is %f \n", nan_value);
	CU_ASSERT_DOUBLE_EQUAL( s_us_pop , 3.3, 0.01 );
	CU_ASSERT_EQUAL( s_phil_pop , 1 );
	CU_ASSERT_TRUE( isnan( s_mars_pop ) ); // nan
	CU_ASSERT_TRUE( isnan( s_mars_us_ratio ) ); // nan
	CU_ASSERT_TRUE( isnan( s_mars_phil_ratio) ); // nan
	CU_ASSERT_TRUE( isnan( s_us_mars_ratio ) ); // nan
	CU_ASSERT_STRING_EQUAL( s_message1, "The relationship data is missing.");
//	printf("%s \n", s_message1);
	CU_ASSERT_STRING_EQUAL( s_message2, "The relationship data is missing.");
//	printf("%s \n", s_message2);
	CU_ASSERT_STRING_EQUAL( s_message3, "The relationship data is missing.");
//	printf("%s \n", s_message3);

	// Clean up
	sailr_tree_free(ps);
	sailr_ptr_table_del_all(&table);
	sailr_parser_state_free(ps);
}


