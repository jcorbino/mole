MOLE = ../mole_C++

CXX ?= g++

INCPATH = -I$(MOLE)
CXXFLAGS = -O3 -fopenmp -DARMA_DONT_USE_WRAPPER
LIBS = -L$(MOLE) -lmole -lopenblas

ifdef ARMA
INCPATH += -I. -I$(ARMA)/include
else
ifeq (,$(filter clean,$(MAKECMDGOALS)))
$(warning The path to Armadillo's library was not provided. I will look for it in the standard system directories.)
endif
endif

ifdef EIGEN
CXXFLAGS += -DEIGEN
INCPATH += -I$(EIGEN)
else
ifdef SUPERLU
LIBS += -L$(SUPERLU)/lib -Wl,-rpath,$(SUPERLU)/lib -lsuperlu
else
LIBS += -lsuperlu
endif
endif

all: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 run

test1: test1.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test1 test1.cpp $(LIBS)

test2: test2.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test2 test2.cpp $(LIBS)

test3: test3.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test3 test3.cpp $(LIBS)

test4: test4.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test4 test4.cpp $(LIBS)

test5: test5.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test5 test5.cpp $(LIBS)

test6: test6.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test6 test6.cpp $(LIBS)
test7: test7.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test7 test7.cpp $(LIBS)
test8: test8.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test8 test8.cpp $(LIBS)
test9: test9.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test9 test9.cpp $(LIBS)
test10: test10.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test10 test10.cpp $(LIBS)
test11: test11.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o test11 test11.cpp $(LIBS)

run: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11
	LD_LIBRARY_PATH=../mole_C++ ./test1
	LD_LIBRARY_PATH=../mole_C++ ./test2
	LD_LIBRARY_PATH=../mole_C++ ./test3
	LD_LIBRARY_PATH=../mole_C++ ./test4
	LD_LIBRARY_PATH=../mole_C++ ./test5
	LD_LIBRARY_PATH=../mole_C++ ./test6
	LD_LIBRARY_PATH=../mole_C++ ./test7
	LD_LIBRARY_PATH=../mole_C++ ./test8
	LD_LIBRARY_PATH=../mole_C++ ./test9
	LD_LIBRARY_PATH=../mole_C++ ./test10
	LD_LIBRARY_PATH=../mole_C++ ./test11

clean:
	rm -f test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11
