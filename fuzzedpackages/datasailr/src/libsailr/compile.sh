# (Dynamic library)
# g++ -shared -o libcpp_string.so cpp_string.o -lstdc++
# gcc -o cpp_string_test cpp_string_test.c -L. -lcpp_string
# export LD_LIBRARY_PATH=. 
# ./cpp_string_test 

# (Static library?)
# ar rvs libcpp_string.a cpp_string.o 

yacc -d -v parse.y  # yacc creates y.tab.c. -d option Creates y.tab.h. -v option creates y.output.
gcc -fPIC -o parse.o -c y.tab.c -DYYDEBUG=1 -Istring

lex -olex.yy.c --debug  lex.l
gcc -fPIC -o lex.o -c lex.yy.c -Istring

gcc -fPIC -o tree_dump.o -c tree_dump.c  -Istring
gcc -fPIC -o tree_free.o -c tree_free.c  -Istring
gcc -fPIC -o node.o -c node.c  -Istring

gcc -fPIC -o parser_state.o -c parser_state.c  -Istring
gcc -fPIC -o var_hash.o -c var_hash.c -Istring
gcc -fPIC -o ptr_table.o -c ptr_table.c  -Istring
gcc -fPIC -o gen_code.o -c gen_code.c -Ivm -Istring
gcc -fPIC -o gen_code_util.o -c gen_code_util.c -Ivm
gcc -fPIC -o vm_label.o -c vm_label.c -lm 
gcc -fPIC -o gen_code_util.o -c gen_code_util.c -Ivm

gcc -fPIC -o sailr.o -c sailr.c  -Ivm  -Istring
# gcc -shared -g -o libsailr.so sailr.o tree_dump.o tree_free.o node.o gen_code.o gen_code_util.o ptr_table.o  parser_state.o var_hash.o vm_label.o string/common_string.o string/cpp_string.o vm/vm.o vm/vm_assign.o vm/vm_calc.o vm/vm_cmd.o vm/vm_code.o vm/vm_stack.o parse.o lex.o 

ar rvs libsailr.a sailr.o tree_dump.o tree_free.o node.o gen_code.o gen_code_util.o ptr_table.o  parser_state.o var_hash.o vm_label.o string/common_string.o string/cpp_string.o vm/vm.o vm/vm_assign.o vm/vm_calc.o vm/vm_cmd.o vm/vm_code.o vm/vm_stack.o parse.o lex.o 


echo "object files are created.  Enjoy!! \n"



