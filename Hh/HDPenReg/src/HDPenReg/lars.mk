## lars makefile for R package

#include $(R_HOME)/etc${R_ARCH}/Makeconf
all:lib

#-----------------------------------------------------------------------
# Variables
# 
LIB = ./lib/libEM.a

#-----------------------------------------------------------------------
# Sources files
#
SRCS =./lars/Lars.cpp \
			./lars/Path.cpp \
			./lars/PathState.cpp \
			./lars/functions.cpp \
 			./lars/Fusion.cpp \
			./lars/Cvlars.cpp \
			./larsRmain.cpp



#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS= $(SRCS:%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
#
%.o: %.cpp	
	$(CXX) $(HD_CXXFLAGS)  $(HD_CPPFLAGS) -DSTKBASEARRAYS=1  $< -c -o $@
#	$(CXX) $(HD_CXXFLAGS)  $(HD_CPPFLAGS) -DSTKBASEARRAYS=1 -DSTK_BOUNDS_CHECK -DFUSION_DEBUG -DCVLARS_DEBUG $< -c -o $@
#	$(CXX) $(HD_CXXFLAGS)  $(HD_CPPFLAGS) -DSTKBASEARRAYS=0 -DSTK_BOUNDS_CHECK -DLARS_DEBUG -DFUSION_DEBUG -DCVLARS_DEBUG $< -c -o $@

#-----------------------------------------------------------------------
# The rule lib create the library 
#
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -r $@ $?
  
mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f $(OBJS)
