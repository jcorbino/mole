# Put here the path to Armadillo
ARMA = /home/johnny/Desktop/armadillo-10.2.1

#----------------------------------------------------------------------
MOLE = ../mole_C++

# Compiler for C++
CXX = g++

# Flags for C++ compiler
CXXFLAGS = -O3 -fopenmp

# Path to header files
INCPATH = -I. -I$(ARMA)/include -I$(MOLE)

# Path to libmole.a
LIBS = -L$(MOLE) -L$(ARMA) -lmole -larmadillo -lsuperlu

# Rules:
all: transport1D elliptic1D

transport1D: transport1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o transport1D transport1D.cpp $(LIBS)

elliptic1D: elliptic1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o elliptic1D elliptic1D.cpp $(LIBS)

clean:
	rm -f transport1D elliptic1D
