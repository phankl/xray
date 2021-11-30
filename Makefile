#Generic c++ Makefile

CXX = g++
CXXFLAGS = -O3 -g -fopenmp -std=c++11
LINKS = -lfftw3 -lm

MAIN = xrayStructure

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

$(MAIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LINKS)

.PHONY: clean
clean:
	rm -f $(OBJ) $(MAIN)
