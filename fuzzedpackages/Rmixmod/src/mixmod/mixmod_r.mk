include $(R_HOME)/etc${R_ARCH}/Makeconf

#-----------------------------------------------------------------------
# Variables
#
LIB = ../mixmod.a
SRC_DIR = ..
#EIGENDIR = "../eigen3"

#-----------------------------------------------------------------------
# Sources files
#
SRCS = $(wildcard */*.cpp) $(wildcard */*/*.cpp) $(wildcard */*/*/*.cpp)

#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS = $(SRCS:%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
# NOTE: $(ALL_CPPFLAGS) contains include path for RcppEigen; do not remove!
#
#%.o: %.cpp
#	$(CXX) -std=gnu++11 $(ALL_CPPFLAGS) -DRPACKAGE -I. -I.. -fpic -march=x86-64 -O2 -pipe $< -c -o $@
PKG_CPPFLAGS = -DRPACKAGE -I. -I..

#-----------------------------------------------------------------------
# The rule lib create the library MIXMOD
#
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rc $@ $?

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o
