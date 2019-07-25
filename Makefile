CXX = g++

STANDARD = -std=c++11
WARNINGS = -Wall
OPTIMIZATION = -g
CPPFLAGS = -I./utility -I/usr/local/Cellar/eigen/3.3.7/include/eigen3 -I./quadrature
CXXFLAGS = $(WARNINGS) $(STANDARD) $(OPTIMIZATION)
LDLIBS = -L./libraries -lquadrature

EXEC = main

SRC = dsmc.cpp configuration.cpp boundary.cpp grid.cpp particles.cpp density.cpp
SRC += potential.cpp force_field.cpp collisions.cpp thermostat.cpp sampling.cpp output.cpp
SRC += $(EXEC).cpp

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
