#------------ LARS MAKEFILE

#-----------------------------------------------------------------------
# Macro definitions. Overwrite or add flags
ifneq ($(wildcard makevars.local),)
  include makevars.local
endif

#-----------------------------------------------------------------------
# source files
CXX_SRCS = $(wildcard *.cpp)

#-----------------------------------------------------------------------
# binary files
CPP_OBJS = $(CXX_SRCS:%.cpp=$(BIN_DIR)%.o)

#-----------------------------------------------------------------------
# Compilation rules
$(BIN_DIR)%.o: %.h makefile
$(BIN_DIR)%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CXXPICFLAGS) $(HD_CXXFLAGS) $< -c -o $@

#-----------------------------------------------------------------------
# Targets
all: obj

#-----------------------------------------------------------------------
# creating objects files
obj: $(CPP_OBJS)
