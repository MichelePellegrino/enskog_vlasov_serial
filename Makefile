CXX = g++
# CXX = gcc

INCLUDE_EIGEN = -I/usr/local/Cellar/eigen/3.3.7/include/eigen3
# INCLUDE_EIGEN = -I/home/matematica/mpellegrino/eigen/3.3.7/include/eigen3

STANDARD = -std=c++11
WARNINGS =
# WARNINGS = -Wall
OPTIMIZATION = -g
CPPFLAGS = -I./utility $(INCLUDE_EIGEN)
CXXFLAGS = $(WARNINGS) $(STANDARD) $(OPTIMIZATION)

# Old version (Romberg static lib.)
# CPPFLAGS = -I./utility -I/usr/local/Cellar/eigen/3.3.7/include/eigen3 -I./romberg
# LDLIBS = -L./libraries -lromberg

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

outclean:
	$(RM) output_files/*.txt
	$(RM) output_files/samples/*.txt
