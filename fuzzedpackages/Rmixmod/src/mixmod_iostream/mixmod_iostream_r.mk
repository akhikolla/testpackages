include $(R_HOME)/etc${R_ARCH}/Makeconf

#-----------------------------------------------------------------------
# Variables
#
LIB = ../mixmod.a
SRC_DIR = ..
IOSTREAM_DIR = ../mixmod_iostream/
#EIGENDIR = "../eigen3"

#-----------------------------------------------------------------------
# Sources files
#
SRCS = $(wildcard *.cpp)

#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS= $(SRCS:%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
#
%.o: %.cpp
	$(CXX) $(CXXFLAGS) ${CPICFLAGS} $(OMPSTATUS) -DXEM_RESOURCES_PATH='"${R_PACKAGE_DIR}/XML_specs"' -I${SRC_DIR} -I${IOSTREAM_DIR} $(LIBXMLXX_CFLAGS) $< -c -o $@

#-----------------------------------------------------------------------
# The rule lib create the library MIXMOD [??!]
#
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rc $@ $?

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o
