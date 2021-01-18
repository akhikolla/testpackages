#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

// Include test functions
#include "test_int.h" 
#include "test_double.h" 
#include "test_arithmetic.h"
#include "test_string.h" 
#include "test_encoding.h" 
#include "test_if.h"
#include "test_regexp.h"
#include "test_missing.h"
#include "test_func.h"
#include "test_date.h"
#include "test_comment.h"

// Prepare add and run tests
int
main()
{
	// Prepare suite 
	CU_initialize_registry();
	CU_pSuite testSuite = CU_add_suite("", NULL, NULL);

	// Suit > Test (> Assert within test functions)
	// Add tests to the suite
	test_int_add_tests(testSuite);
	test_double_add_tests(testSuite);
//	test_arithmetic_add_tests(testSuite);
	test_string_add_tests(testSuite);
	test_encoding_add_tests(testSuite);

	test_if_add_tests(testSuite);
	test_missing_add_tests(testSuite);
	test_regexp_add_tests(testSuite);
	test_func_add_tests(testSuite);
	test_date_add_tests(testSuite);

	test_comment_add_tests(testSuite); 

	// Run CUnit
//	CU_console_run_tests();  // This is interacitve version  #include <CUnit/Console.h>
//	CU_automated_run_tests();  // This outputs XML file. #include <CUnit/Automated.h>
	CU_basic_run_tests();  // This prints out results on console.  #include <CUnit/Basic.h>

	// Finish
	CU_cleanup_registry();

	return 0;
}
