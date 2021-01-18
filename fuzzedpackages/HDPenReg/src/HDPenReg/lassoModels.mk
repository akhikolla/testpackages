## lars makefile for R package

#include $(R_HOME)/etc${R_ARCH}/Makeconf
all:lib

#-----------------------------------------------------------------------
# Variables
#
LIB= ./lib/libHD.a

#-----------------------------------------------------------------------
# Sources files
SRCS=./lassoModels/LassoSolver.cpp \
			./lassoModels/FusedLassoPenalty.cpp \
			./lassoModels/FusedLassoSolver.cpp \
			./lassoModels/LogisticLassoSolver.cpp \
			./lassoModels/LogisticFusedLassoSolver.cpp \
			./lassoModels/CV.cpp \
			./EMmain.cpp \
			./EMCVmain.cpp
#			./lassoModels/EnetPenalty.cpp \


#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
OBJS= $(SRCS:%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
%.o: %.cpp
	$(CXX)  $(HD_CXXFLAGS) $(HD_CPPFLAGS) $< -c -o $@

#-----------------------------------------------------------------------
# The rule lib create the library
#
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -r $@ $?

mostlyclean: clean

clean:
	@-rm -f $(OBJS) $(LIB)
