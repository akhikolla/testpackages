g++ -fPIC -std=c++11 -c cpp_date.cpp -o cpp_date.o
gcc -fPIC -c simple_date.c -o simple_date.o
gcc  test.c simple_date.o cpp_date.o -o test_run  -lstdc++
