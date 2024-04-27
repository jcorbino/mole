MOLE = ../mole_C++

CXX ?= g++

CXXFLAGS = -O3 -fopenmp -DARMA_DONT_USE_WRAPPER

ifdef ARMA
INCPATH = -I. -I$(ARMA)/include -I$(MOLE)
else
ifeq (,$(filter clean,$(MAKECMDGOALS)))
$(error ARMA is not defined. Please define ARMA with the path to the Armadillo library)
endif
endif

LIBS = -L$(MOLE) -lmole -lopenblas

all: test1 test2 test3 test4 run

test1: test1.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test1 test1.cpp $(LIBS)

test2: test2.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test2 test2.cpp $(LIBS)

test3: test3.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test3 test3.cpp $(LIBS)

test4: test4.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test4 test4.cpp $(LIBS)

run: test1 test2 test3 test4
	LD_LIBRARY_PATH=../mole_C++ ./test1
	LD_LIBRARY_PATH=../mole_C++ ./test2
	LD_LIBRARY_PATH=../mole_C++ ./test3
	LD_LIBRARY_PATH=../mole_C++ ./test4

clean:
	rm -f test1 test2 test3 test4
