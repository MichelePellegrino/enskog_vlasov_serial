CXX = g++

STANDARD = -std=c++11
WARNINGS = -Wall
OPTIMIZATION = -g
CXXFLAGS = $(WARNINGS) $(STANDARD) $(OPTIMIZATION)

RESULT = ../libraries/libromberg.a

SRC = poly_interp.cpp base_interp.cpp
OBJS := $(SRC:%.cpp=%.o)
HEADERS = $(wildcard *.hpp)

interpolation.o: recipe_types.hpp interpolation.hpp

.DEFAULT_GOAL = all

.PHONY: all clean distclean

all: $(RESULT)

$(RESULT): $(OBJS) $(HEADERS)
	ar rs $@ $^

clean:
	$(RM) $(RESULT)

distclean:
	$(RM) $(RESULT)
	$(RM) *.o

cleanobj:
	$(RM) *.o
