YACC = bison -y
LEX = flex
CC = gcc
CXX = g++
AR = ar rvs
RM = rm -f

CFLAGS := $(CFLAGS) -std=c11 -g -O3 -I. -Ivm -Istring -Isimple_re -Isimple_date -I$(ONIG_INCLUDE_DIR) $(CC_USER_DEFINES)
CXXFLAGS := $(CXXFLAGS) -std=c++11 -g -O3 -Ivm -Istring $(CXX_USER_DEFINES) 

ARCH_TRIPLE := $(subst -, ,$(shell $(CC) -dumpmachine))

ifneq (,$(findstring mingw,$(ARCH_TRIPLE)))
# mingw
else
    ifneq (,$(findstring linux,$(ARCH_TRIPLE)))
		# linux
        CFLAGS := $(CFLAGS) -fstack-protector-strong -fPIC
        CXXFLAGS := $(CXXFLAGS) -fstack-protector-strong -fPIC 
    else
		# e.g. daarwin
        CFLAGS := $(CFLAGS) -fPIC
        CXXFLAGS := $(CXXFLAGS) -fPIC 
    endif
endif
$(info SHOW ARCH_TRIPLE)
$(info $(ARCH_TRIPLE))
$(info SHOW CFLAGS)
$(info $(CFLAGS))

TARGET = libsailr.a

SRCS=$(filter-out y.tab.c lex.yy.c, $(wildcard *.c)) $(wildcard vm/*.c) 
OBJS:=$(SRCS:.c=.o) parse.o lex.o
SRCS_CXX_STR:=$(wildcard string/cpp_string/*.cpp)
OBJS_CXX_STR:=$(SRCS_CXX_STR:.cpp=.o) 
OBJS_STR:= string/common_string.o 
OBJS_RE:= simple_re/simple_re.o
OBJS_DATE:= simple_date/simple_date.o simple_date/cpp_date.o
SRCS_CFUNC:= $(wildcard vm/func/c_func/*.c)
OBJS_CFUNC:=$(SRCS_CFUNC:.c=.o)
SRCS_EXTFUNC:= $(wildcard vm/func/ext_func/*.c)
OBJS_EXTFUNC:=$(SRCS_EXTFUNC:.c=.o)
DEPS:=$(OBJS:.o=.d) $(OBJS_STR:.o=.d) $(OBJS_RE:.o=.d) $(OBJS_DATE:.o=.d) $(OBJS_CFUNC:.o=.d) $(OBJS_EXTFUNC:.o=.d) $(OBJS_CXX_STR:.o=.d)

# List "objfile depends on source and header"
-include $(DEPS)

.PHONY: build test clean distclean

build : parse.o lex.o  $(TARGET) 

$(TARGET) : $(OBJS) $(OBJS_CXX_STR) $(OBJS_STR)  $(OBJS_RE) $(OBJS_DATE) $(OBJS_CFUNC) $(OBJS_EXTFUNC)
	$(AR) $(TARGET) $(OBJS) $(OBJS_CXX_STR) $(OBJS_STR) $(OBJS_RE) $(OBJS_DATE) $(OBJS_CFUNC) $(OBJS_EXTFUNC)

parse.o : y.tab.c
	$(CC) -o parse.o -c y.tab.c $(CFLAGS) -MMD -MP 

y.tab.c : parse.y
	$(YACC) -d -v --verbose --debug parse.y
    # yacc creates y.tab.c. -d => y.tab.h. -v => y.output

lex.o : lex.yy.c y.tab.h
	$(CC) -o lex.o -c lex.yy.c $(CFLAGS)  # This file requires y.tab.h

lex.yy.c : lex.l
	$(LEX) -olex.yy.c lex.l

%.o : %.c 
	$(CC) -c -o $@ $<  $(CFLAGS) -MMD -MP

vm/%.o : vm/%.c 
	$(CC) -c -o $@ $<  $(CFLAGS) -I. -MMD -MP

string/cpp_string/%.o : string/cpp_string/%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) -MMD -MP 

string/%.o : string/%.c 
	$(CC) -c -o $@ $<  $(CFLAGS) -I. -MMD -MP

simple_re/%.o : simple_re/%.c
	$(CC) -c -o $@ $<  $(CFLAGS) -I. -MMD -MP

simple_date/cpp_date.o : simple_date/cpp_date.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) -MMD -MP 

simple_date/%.o : simple_date/%.c
	$(CC) -c -o $@ $<  $(CFLAGS) -I. -MMD -MP

test : 
	@echo "\033[1;34m Make sure you have rebuilt $(TARGET) before running tests \033[0m"
	@echo "\033[1;34m CUnit tests for $(TARGET) \033[0m"
	$(CC) test/test_main.c $(TARGET) -o test/a.out -Wall -I. -lstdc++ -lm -L/usr/local/lib -lcunit -Ldev_env/onigmo_build/lib -lonigmo
	test/a.out

distclean :
	$(RM) *.o
	$(RM) vm/*.o
	$(RM) string/*.o
	$(RM) string/cpp_string/*.o
	$(RM) simple_re/*.o
	$(RM) simple_date/*.o
	$(RM) vm/func/c_func/*.o
	$(RM) vm/func/ext_func/*.o
	$(RM) test/*.o
	$(RM) *.d
	$(RM) vm/*.d
	$(RM) string/*.d
	$(RM) string/cpp_string/*.d
	$(RM) simple_re/*.d
	$(RM) simple_date/*.d
	$(RM) vm/func/c_func/*.d
	$(RM) vm/func/ext_func/*.d
	$(RM) test/*.d
	$(RM) *.dll
	$(RM) *.so
	$(RM) vm/*.so
	$(RM) string/*.so
	$(RM) string/cpp_string/*.so
	$(RM) simple_re/*.so
	$(RM) simple_date/*.so
	$(RM) vm/func/c_func/*.so
	$(RM) vm/func/ext_func/*.so
	$(RM) test/*.so
	$(RM) *.a
	$(RM) vm/*.a
	$(RM) string/*.a
	$(RM) string/cpp_string/*.a
	$(RM) simple_re/*.a
	$(RM) simple_date/*.a
	$(RM) vm/func/c_func/*.a
	$(RM) vm/func/ext_func/*.a
	$(RM) test/*.a
	$(RM) a.out
	$(RM) vm/a.out
	$(RM) string/a.out
	$(RM) string/cpp_string/a.out
	$(RM) simple_re/a.out
	$(RM) simple_date/a.out
	$(RM) vm/func/c_func/a.out
	$(RM) vm/func/ext_func/a.out
	$(RM) test/a.out
	$(RM) core
	$(RM) vm/core
	$(RM) string/core
	$(RM) string/cpp_string/core
	$(RM) simple_re/core
	$(RM) simple_date/core
	$(RM) vm/func/c_func/core
	$(RM) vm/func/ext_func/core
	$(RM) test/core
	$(RM) test/*.xml
	$(RM) $(TARGET)

clean: distclean
	$(RM) y.tab.c y.tab.h y.output
	$(RM) lex.yy.c



