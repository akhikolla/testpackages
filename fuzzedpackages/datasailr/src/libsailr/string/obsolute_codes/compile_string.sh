g++ -fPIC -std=c++11 -o cpp_string.o -c cpp_string.cpp
# gcc cpp_string_test.c cpp_string.o -lstdc++  -o cpp_string_test -I.

gcc -fPIC -o common_string.o -c common_string.c 

gcc -fPIC -o struct_string.o -c struct_string.c 

echo "String libraries are created."
