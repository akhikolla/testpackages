LIB = libcoclust.a
SRCS = $(wildcard src/Algorithms/*.cpp\
                  src/CoClustFacade/*.cpp\
                  src/Strategy/*.cpp\
                  src/InputParameters/*.cpp\
                  src/Models/*.cpp)
OBJS= $(SRCS:%.cpp=%.o)
%.o: %.cpp
	$(CXX) $(CO_CXXFLAGS) $(CO_CPPFLAGS) $< -c -o $@
	
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rc $@ $?
  
.PHONY: clean

clean:
	@-rm $(OBJS)

clean_all:
	@-rm $(OBJS) $(LIB)
