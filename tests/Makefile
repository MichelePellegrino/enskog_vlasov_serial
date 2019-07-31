CXX = g++

STANDARD = -std=c++11
WARNINGS = -Wall
OPTIMIZATION = -g
CPPFLAGS = -I../utility -I../quadrature -I../romberg
CXXFLAGS = $(WARNINGS) $(STANDARD) $(OPTIMIZATION)
LDLIBS = -L../libraries -lquadrature -lromberg

EXEC = main

SRC = test_integration.cpp

OBJS = $(SRC: .cpp = .o)

make.dep: $(SRC)
	$(RM) make.dep
	for f in $(SRCS); do \
  	$(CXX) $(CPPFLAGS) -MM $$f >> make.dep; \
  done
-include make.dep

.DEFAULT_GOAL = all

.PHONY: all clean distclean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDLIBS) $(OPTIMIZATION) $^ -o $@

clean:
	$(RM) $(EXEC)

distclean:
	$(RM) $(EXEC)
	$(RM) *.o *.dep
	$(RM) -r *.dSYM
