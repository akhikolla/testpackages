# HISTORY

## Project Started [Aug 2018]

* I started this project to make it easier to manipulate dataset.
* First, I tried to make a simple calculator


## Milestone Ver 0.1 [Nov 2018]

* The program starts woring.
	+ Virtual stack machine is implemented.
	+ Arithmetic calculations are implemented. 
	+ How to compile & run

```
cd vm; ./compile_vm.sh; cd .. 
./compile.sh
./myparser sample_code/code1.slr
```

## Ver 0.2 [Nov. 23 2018]

* The library now deals with C++ strings.
    + cpp_string.cpp implements that part.
    + struct_string.cpp is now not updated. (But remains there)
    + If this library might be used from C program, mechanism to switch between struct_string and cpp_string shoul be implementd
* The library is now used from C++.
    + C++ main program can calculate using this library. See main.cpp
* Makefile (GNU Make) can compile 
    + Meaning that this can be easily integrated with Rcpp project.

```
./compile_cpp.sh
```

* The project name is now changed from RCppCalc to LibsailrDev.
* To use libsailr part, use only under LibsailrDev/sailr/ directory.

```
cd <LibsailrDev_Project_Directory>/sailr/
make build
```


## Ver 0.21 [Nov. 28 2018]

* Bug fixes
    + string manipulations are improved. 
* From Rcpp this library is successfully loaded and utilized.

## Ver 0.22 [Dec. 17 2018]

* Bug fixes
    + Nested if else did not work well. Now it's fixed.
    + working fine.

## Ver 0.23 [Dec. 18 2018]

* Newly created variables are supported.
    1. Variable should be prepared as PTR_NULL on ptr_table.
    2. The variable type should be changed dynamically.


## Ver 0.24 [Dec.21 2018]

* Bug fix: double was not dealt properly, and it's fixed.

## Ver 0.25 [Dec.21 2018]

* Change directory structure.
    + Easier to be used from other projects.

## Ver 0.26 [Dec. 21 2018]

* For undefined varialbe, its type is assigned to the type of first value.
* Number (int and double) are dealt in the same way at VM instruction level.

## Ver 0.30 [Dec. 25 2018]

* Numbers are converted between double and int.
	* Add extra address field to ptr_table/ptr_record.
	* This enables number-type record to keep memory for int and double.
* This happens when assingning and LHS type is different from RHS type.

* Add unitary minus operator support.

## Ver 0.31 [Dec. 28 2018]

* Minor fixes
	* Function ptr_table_get_pptr() 
	* Some bug fixes.

## Ver 0.32

* Minor fixes for R package.
    * Code21, in which strings are manipulated, passed. 

## Ver 0.33 [Jan. 16 2019]

* Minor fixes
    + Add missing values support. (Period is dealt as missing number.)

## Ver 0.35 [Feb. 18 2019]

* Regular expression support 
* Function support

## Ver 0.40 [Mar. 5 2019]

* Functions are supported. (Just a few functions)
    + Built-in functions (written in C) are avaible.
    + Only print() function is implemented. 


## Ver 0.50 [Apr. 6 2019]

* More functions are added
    + print( str1, str2 ... )
    + num_to_str( num ) 
    + str_strip( str )
    + str_lstrip( str ) 
    + str_rstrip( str )
    + str_concat( str1, str2 ... ) 
    + str_repeat( str )
    + str_subset( str, num, num )  // index starts from zero.
    + str_to_num( str )
    + rexp_matched( num )  // back reference
    + date_ymd( year, month, day )
    + date_ym_weekday_nth( year, month, weekday, nth )  // weekday should be "Sun", "Mon" ... 
    + date_add_n_years( unix_date, years )
    + date_add_n_months( unix_date, months )
    + date_add_n_days( unix_date, days )
    + date_format( unix_date, str_format )  // This format should follow C++'s std::chrono::format. "%m/%d/%y" or "%Y-%m-%d" is popular.
* Back reference mechanism for regular expresion
    + Use function to do this. => rexp_matched(1)

## Ver 0.51 (Apr. 21 2019)

* Refactoring1
    + Stop Memory Leak
    + Free pointer table at last

## Ver 0.60 (May. 15 2019)

* Major improvements.
* Refactoring2
    + Introduce DEBUG_PRINT macro (in helper.h) for C part.
* Improve vm_calc.c (1)
    + Take care of INT_MIN and INT_MAX for integers. Performance is not efficient, but before int calculation, double calculation and integer range check is conducted.
    + deal with missing values.
        + At parsing phase, missing values are already treated as nan in double.
        + Be careful not to convert it into integer unexpectedly.
* Improve vm_calc.c (2) , vm/func/c_func , vm_assign.c , vm_rexp.c and vm_stack.c
    + Prevent from dealing with stack pointer directly. 
        + Use stack vm_stack_push_*** and vm_stack_clean_and_pop()
    + For vm_assign.c, it is refactored.
* Rethink about libsailr API. (for users to prevent memory leak easily.)
    + "Users should use pointers like PTR_INT, PTR_DBL and PTR_STR on ptr_table."
        + User should not use IVAL, DVAL. 
    + *Users should usually prepare memory for those known variables* , and they should be freed manually.
* Made vm_assign.c tidy
* Refactor pp2val.c (using vm_stack_item_is_temp() function)
* Refactor ptr_table: ptr_record_free_memory_if_gc_required() function. This should use string_free(), simple_re_free(). 

## Ver 0.61 (Jun. 1 2019)

* Unit tests are introduced. (CUnit)
* Along with the tests, bug fixes wree done.
	+ Examples:
	+ fcall was reduced to expr, but now is reduced to arg (in parse.y)
	+ Operations for missing values are corrected.
		+ For eq(==), when both sides have nan, return true. If either side has nan, return false. For neq(!=), the behavior is opposite.
		+ For other operators, like + - * =, follow C math.h implementation, meaning returning false
			+ https://stackoverflow.com/questions/38798791/nan-comparison-rule-in-c-c
	+ Dealing with comments is difficult.
		+ (Start) states are explicitly defined for most of the rules.
		+ One line comments and multiple line comments are implemented.
		+ One line comments should work like a terminator.
		+ Multiple line comments should start in a new line at current implementation. DON'T insert it wihtin a statement.
			+ (ref.) http://www.cs.man.ac.uk/~pjj/cs2121/ex2_str_comm.html 

## Ver 0.62 (Jun. 16 2019)

* Start to support 32bit/64bit windows platform using mingw-w64.
	+ See build scripts/envs under mingw_env directory.
* Period(.) is allowed to be included for variable name or column name. (iris data includes period.)

## Ver 0.63 (Aug. 8 2019)

* Move 32bit/64bit windows build scripts outside of this repository.
* Added compiler flags
    + Explicitly use "-std=c99", "std=c++11" for compilers
    + Use "-g" for debugging purpose (This does not slow donw execution)
    + "-fstack-protector-strong" for C/C++ compilers
        + This flag is disabled for mingw compiler.
        + See the conditional in makefile.
    + These are learnt from testthat library compilations.
* In Makefile
    + $^ is replaced with $<
        + $^ : All the prerequisites
        + $< : The 1st prerequisite
        + The Makefile uses -include $(DEPS) & -MD -MMP mechanism, meaning that header files are set to be prerequisites for target file.
            + Only the 1st prerequisite should be compiled. Prerequisite at 2nd and after are header files, so they should not be listed in gcc/g++ arguments.
* Under dev_env directory
    + Makefile is used for test parser (myparsercpp) compilation.
    + Onigmo source for this test parser is put under dev_env/onigmo_src
        + This onigmo_src is ignored from git. 
        + Download onigmo source and extract the files there. 


## Ver 0.64 (Sep. 10 2019)

* Parser (bison/flex) is now reentrant.
* (Note) Under simple_re directory, global variable still exists to store the regular expression that was last executed.
    + In Ver 0.74, this point is updated.

## Ver 0.65 (Sep.11 2019)

* Global variables in ptr_table.c are removed. 
    + Those counter information for anonymous string and regexp are stored in the 1st element of ptr_table object.
    + The type is defined in ptr_table_info that is stored in the seed element called "_HEAD_OF_UTHASH_"

## Ver 0.66 (Sep.23 2019)

* ptr_table now holds information whether null variable (Known but not defined) obtains type definition.
    + ptr_table_info->null_update holds this information. This does not holds which variable became defined.


## Ver 0.67 (Sep.24 2019)

* vm/func/c_func.c is updated.
    + Previous codes created anonymous string on ptr_table every time the function returns new string even within RHS of assignment operator.
    + The current code jsut creates temporary string on vm stack as return from function.
    + Minor related fixes.

## Ver 0.68 (Sep.25 2019)

* vm_stack is now able to hold string encoding information. 
    + Based on this information, vm stack is going to call appropriate string functions.
    + Add test case for this functionality.

## Ver 0.69 (Sep. 30 2019)

* UTF8 support. The following functions can deal with utf8 strings
    + cpp_string_subset() calls appropriate functions based on encodings. (Now this fuction requires encoding.)
        + cpp_string_subset_utf8()
        + cpp_string_subset_latin1()
    + cpp_string_new_unescaped_string() calls appropriate functions based on encodings. (Now this fuction requires encoding.)
        + cpp_string_new_unescaped_string_utf8()
        + cpp_string_new_unescaped_string_latin1()
* vm_stack holds information how to deal with strings.
    + calls appropriate functions.
* parser_state object now holds source file encoding, which can be used for creating new regular expression object.
    + regular expression object continues to hold this encoding information.
    + string objects do not need this information. They just hold byte sequences that come from input data or source file string literals.

## Ver 0.70 (Oct. 6 2019)

* Avoid (char*) casting
    + ptr_table functions take const char* as key.
        + ,though the key string of ptr_record continues to be char[].
    + simple_date_format()


## Ver 0.71 (Oct. 25 2019)

* Minor fixes for codes showing VM instructions. (vm_inst_list_show_all() uses printf(), not DEBUG_PRINT(). )
* Minor fixes for codes dumping parser tree. (tree_dump() uses printf(), not DEBUG_PRINT(). )
* Minor fixes for lexer. The following if-else is now allowed.
    + Previoulsy, line beginning directly with "else" was not allowed. Only "} else" was allowed.

```
if ( carname =~ re/(^Merc)/ ) { country = "Germany" ; type = rexp_matched(1) }
else if( carname =~ re/(^Cadillac|^Ford)/ ) { country = "USA" ; type = rexp_matched(1); }
else if( carname =~ re/(^Honda|^Toyota)/ ) { country = "Japan" ; type = rexp_matched(1); } 
else { carname = "other country" }
```


## Ver 0.72 (Oct. 26 2019)

* From Ver 0.66, ptr_table holds information whether null variable was updated. 
    + Now it holds which type is newly created, PTR_INT, PTR_DBL, PTR_STR or PTR_REXP
    + The default value of null_update is set to 0b0000.
        + From right, the 1st bit represents PTR_INT, the 2nd bit PTR_DBL, the 3rd bit PTR_STR and the 4th PTR_REXP. The rest of bits are not used.
        + If a bit is set to 1, the corresponding type is newly created.
* ptr_table_info_reset_null_update() is added.
    + This function is available through sailr.h


## Ver 0.72b 

* Minor fix. Additional fix for Ver 0.71.
    + Resolved reduce/reduce conflicts.


## Ver 0.73 (Oct. 29 2019)

* libsailr interface for adding string onto ptr_table is changed.
    + sailr_ptr_table_create_string_from_ptr is deprecated. 
    + Instead, sailr_ptr_table_create_string_from_cstring is introduced. 
        + Reasons for this change. 
        + Strings need to be tracked adn freed at an appropriate timing. At current implementation, the mechanism is not implemented enough.
        + To make this software available in public as soon as possible, stability should be preferred to performance. It should never happen to break user data.


* Current rules when using libsair
    1. For int and double, pass int or double pointer to libsailr. The value the pointer points to will be updated. 
        + You can obtain the calculation result by dereferencing the pointers. 
        + sailr_ptr_table_create_int_from_ptr() and sailr_ptr_table_
    2. For string, pass the initial value as cstring (constant char*). Meaning the original string objects are never modified or destroyed during libsailr calculation.
        + sailr_ptr_table_create_string_from_cstring() pass the initial string value.
        + const char* sailr_ptr_table_get_cstring() to obtain the result. 


## Ver 0.74 (Nov. 5 2019)

* The global variable was found that should have been deleted in Ver.0.64. 
    * simple_re* re_last_matched in simple_re.h is removed.
    * Instead, vm_stack now holds the last executed regular expression as vmstack->stack[0].p_vm_stack_info->last_rexp.
    * Related functions in simple_re.c/.h, vm_rexp.c, vm_stack.c, vm/func/c_func/c_func.c are updated.

## Ver 0.75 (Nov. 6 2019)

* Regular expressions should be reset every time calculation finishes from library user.
* To enable this, functions to extract specific type of records (in this case PTR_REXP) and reset regular expressions are implemented. 
    * ptr_record_get_type()
    * ptr_record_next()
    * ptr_table_first_record()
    * ptr_record_reset_rexp()
* Also corresponding sailr functions are implemented.


## Ver 0.76 (Nov. 9 2019)

* Fix parsing of if_stmt
    + lex.l and parse.y are updated.
    + Since the Ver 0.71 update, "if statement" is defined to take optional TERMIN (=opt_termin) between then_stmts and opt_else.
        + This resulted in if-statement-without-else not working. When opt_else is empty, opt_termin matches the end of if(){} statement, and TERMIN is lost between the current if-statemnt and next statement.   
    + This is resolved by removing opt_termin from if_statement in parse.y. Instead, in lex.l, else token is redefined to be [\t \n]*else[\t \n]* as follows.


* Excerpted from the output of "git diff HEAD"

```
// Main change in lex.l
-<INITIAL,IFSTATE,ELSESTATE>else
+<INITIAL,IFSTATE,ELSESTATE>[\t \n]*else[\t \n]*

// Main change in parse.y
-if_stmt        : KEY_IF condition then_stmts opt_termin opt_else
+if_stmt        : KEY_IF condition then_stmts opt_else
```

* Now the following code works.

```
if(condition){then_statement} TERMIN
next_normal_statment
```

## Ver 0.77 (Nov. 16 2019)

* str_subset()'s arguments are now one-indexed.
    + (e.g.) str_subset("Hello World", 1, 5) returns "Hello"
* print() function can now take not only string but numbers (integer + double).
    + Some new functions are added to common_string and cpp_string.


## Ver 0.78 (Nov. 17 2019)

* Division calculation is updated. It now always returns double.
    + Generating Inf from division is now properly handeled.

## Ver 0.79 (Jan. 5 2020)

* Assignment operator did not work when the stack item correspoding to RHS of assignment, i.e. top item on stack , is still PP_INT or PP_DBL. 

```
# e.g.
# Suppose age variable already exists on ptr_table as PTR_INT
# The following code did not work, because age on stack is still PP_INT
age2 = age
```

* To solve this, vm/vm_assign.c is updated to convert PP_IVAL/PP_DVAL to IVAL/DVAL for the top item of stack.
	+ Now, stack_item_pp2value(rvalue) is called every time assignment operation is conducted, which converts PP_IVAL/PP_DVAL to IVAL/DVAL for rvalue on stack.


## Ver 0.80 (Jan. 19 2020)

* API function name is renamed. Files using this function are updated.
    + sailr_construct_parser() => sailr_run_parser() 
* Copyright files are updated (Jan. 23 2020)


## Ver 0.81 

* Resolving warnings, and improving for portability.
* Makefile is updated (Feb. 4 2020) 
    + Step to generate C file from lex file is seperated.
    + "lex.o: lex.y.c y.tab.h" , y.tab.h is added as prerequisite for lex.o target.
    + "lex.yy.c : lex.l" rule is added.
        + Even when only parse.y is updated (which triggers $(YACC) commnd and generates y.tab.h and y.tab.c), and as a result lex.o is regenerated. lex.yy.c uses y.tab.h, so when y.tab.h is updated, lex should be run.
    + With this update, source codes can be distributed for machines without bison and flex.
        + Run "make y.tab.c" and "make lex.yy.c" before disribution.
* Binary files are removed from git source tree. (Feb. 5 2020)
* Strncpy is changed to memcpy  (Feb. 5 2020)
    + strncpy is changed to memcpy in lex.l. Dynamically the length of the original string is obtained, and allocate exact memory for the length + 1. (+1 if for null terminator.) In this case, strncpy and memcpy works the same, and memcpy is faster.
    + strncpy is chnaged to memcpy also in simple_re/simple_re.c.
    + strncpy is chnaged to memcpy also in gen_code.c. This is copying string into array that has enough size.
* Makefile is updated for compilation on CRAN. (Feb. 8 2020) 
* Some systems use macro function for memcpy definition. In vm_stack.c, compound literal was passed for the second argument, and macro function wrongly seperate those compound literals b/c they include comma within them. I put () parentheses for the second argument. (Feb. 9 2020)
* For g++, -D_GLIBCXX_USE_CXX11_ABI=0 option is added in Makefile to avoid errors in some environment (Mar. 8 2020)
    + About -D_GLIBCXX_USE_CXX11_ABI=0 option, see https://stackoverflow.com/questions/33394934/converting-std-cxx11string-to-stdstring
    + This option is discarded later. Instead use -std=c11 option to C compiler. (Mar. 22 2020)


## Ver 0.8.2

* Makefile updated to detect operating system. (gcc -dumpmachine is used) (Mar. 13 2020)
* The following declarations are added to prevent implicit declaration warnings. (Mar. 13 2020)
    + parse.y: yylex() and yyerror() are explicitly declared.
    + lex.l: fileno(FILE *stream) is explicitly declared.
    + sialr.c: yylex_init(), yy_scan_string(), and yylex_destroy() are explicitly declared.
* Use -std=c11 for C, as well as -std=c++11 for C++ (Mar. 22 2020)
    + Instead, discard -D_GLIBCXX_USE_CXX11_ABI=0 option.
* COPYRIGHT is updated. (Mar. 22 2020)


## Ver 0.8.3

* Solve valgrind errors (Apr. 5-11 2020)
    + Fix memory leaks.
    + Avoid multiple deallocation for the same memory.
    + Memory for string was unintentionally freed.
* Warnings by -Wreturn-type -Wparentheses -Wunused-value are solved. (Apr. 12 2020)
* Macro functions are changed to use do{}while(0). For portable codes, statement expressions are removed. (Apr.16 2020)
* Anonymous unions are changed to be used only in c11 or higher (__STDC_VERSION__ >= 201112L). (Apr.16 2020)
* Binary literals are now only used in C11 or higher.
* C++ errors are now catched by reference.
* sprintf() is removed. snprintf() is now used to avoid warnings.


## Ver 0.8.4

* Warnings are resolved (Apr.30-May.4 2020)
* Location information (not only line but also column number) becomes available in lex.l and parse.y. (May.6 2020)
    + YYLTYPE* is used for this purpose.
    + yyerror() now prints out line number and column number, when encountering syntax error.
* When VM detects runtime errors, it now stops further execution.
    + Currently, two mechanisms are used for this purpose. (May. 8 2020)
        1. Function return values (1:success, 0:fail)
        2. vm_error_raise() and vm_error_exist() in vm_error.c (available from this version)
            + The former is more strainghtforward. The latter is useful when function calls (caller/callee relationships) are complex.
* Corresponding script position (location) is reported when runtime error happens. (long-awaited feature!)
    + `struct script_loc` (defined in script_loc.h) holds script location information.
    + `struct TreeNode_` (i.e. TreeNode) and `struct _vm_inst` (i.e. vm_inst) now have location field, `struct script loc`.
        + TreeNode's loc field is set during node construction in parse.y. 
        + vm_inst's loc field is set during converting node tree to vm instructions in gen_code.c
    + When runtime error is deteted in vm_exec_code in vm.c, the location is reported.
* anonym field is added to ptr_record (=ptr_table).
    + When anonym is set 1, the record is generated as an anonymous object (e.g. strings or regular expressions generated from literals).
    + This enables library users to choose ptr_record's that are not assigned to varialbes (anonymous objects are not assigned to variables).
* gc field (and ex_gc field) on ptr_record (=ptr_table) is effectively used.
    + gc (and ex_gc field) field tells whether the object managed by the ptr_record needs to be freed (destroyed) or not when it is detached from the ptr_record.
        + Objects are usually detached when new object is assigned or each execution finishes. At this time, based on gc field, they are freed.
            + For this purpose, ptr_record_free_gc_required_memory() is used.
            + ptr_record_free_gc_required_memory() becomes available also for library users via sailr_ptr_record_free_objects().
        + Objects that come from outside of library usually need not be freed, so their gc fields are set GC_NO.
* Variable with PTR_REXP type is now allowed to be on left hand side of assignment. This allows assigning regular expression object to variable more than once.
    + Currently, assignment of regular expression creates new regualr expression objects. Therefore, note that this may affect performance.


## Ver 0.8.5

* No chnage. Only version number is updated.


## Ver 0.8.6

* Mechanism to allow externally defined functions is introduced. 
    + Currently, function pointers with the type of int (*func)( arg_list* , unsigned int, vm_stack* ) can be registered using sailr_ext_func_hash_add() function, which stores pairs of Sailr function name, number of expected arguments and corresponding C function pointer. The C function needs to manipulate items on VM stack, so more header files (such as sailr_ext.h) need to be included by library user.
* vm_exec_code() returns a new status code, 2, which means suspend. 0:fail, 1:success, 2:suspend.
* vm stack (vm_stack) now holds information about the next instruction code position when suspended, which is used when vm resumes.
* ext_func_hash holds the last executed function name, which is useful when vm suspends and you need to know which function has caused that suspension.
* Bug fix
    + When a variable with string type is assigned to a variable with the same name (i.e. itself), the string object was garbage collected. Now it's fixed, and garbage collection does not happen. (Nothing needs to happen in this case.)
    + ptr_table_update_int() and ptr_table_update_double() now update values, not pointers.


## Ver 0.8.7

* Bug fix
    + tree_free() tried to free tree nodes, even when parser_state's tree points to NULL. It does not do anything and just return in such a case.


## Plan 

* Avoid directly manipulate ptr_table's properties. Provide functions and use them.
* Consider some script language extension. BSD licensed language is best (e.g. Lua, mruby or Gauche??)
* Macro to add variables for users to ptr_table.
    + When adding value to ptr_table, missing values should be taken care of.
        + Missing value should be added as nan in double.
* Refactoring2
    + Functions in ptr_table.c. Pointer to pointer may be used wrongly; possibility for some local pointers are destroyed unintentionally.


## Abandoned Ideas

* Ways to deal with general objects.
    + At ptr_table level?? Only at function call level?? => Implement only at function call level.
    + For example. let's think about how to deal with tm structure defined in time.h in C.
        + At ptr_table level, PTR_OBJ is used.
        + At vm stack level, PP_OBJ is used.
        + PP_OBJ pointts to a wrapper for some structure. 
        + That wrapper (wrapper_obj) holds the type of real object.


