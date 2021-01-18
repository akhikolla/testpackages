#!/bin/sh

gcc -c simple_re.c -o simple_re.o -I./ -I../ -I../dev_env/onigmo_build/include -L../dev_env/onigmo_build/lib -lonigmo

gcc test.c simple_re.o -o test_run  -I./ -I../ -I../dev_env/onigmo_build/include -L../dev_env/onigmo_build/lib -lonigmo

