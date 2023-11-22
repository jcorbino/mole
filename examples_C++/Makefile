MOLE = ../mole_C++

# Compiler for C++
CXX = g++

# Flags for C++ compiler
CXXFLAGS = -O3 -fopenmp -DARMA_DONT_USE_WRAPPER

# Path to header files
INCPATH = -I. $(if $(ARMA),-I$(ARMA)/include) -I$(MOLE)

# Path to libmole.a
LIBS = -L$(MOLE) -lmole $(if $(ARMA),-L$(ARMA)/lib) -lopenblas

ifdef EIGEN
CXXFLAGS += -DEIGEN
INCPATH += -I$(EIGEN)
else
LIBS += -lsuperlu
endif

# Rules:
all: transport1D schrodinger1D elliptic1D elliptic2D

transport1D: transport1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o transport1D transport1D.cpp $(LIBS)

schrodinger1D: schrodinger1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o schrodinger1D schrodinger1D.cpp $(LIBS)

elliptic1D: elliptic1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o elliptic1D elliptic1D.cpp $(LIBS)

elliptic2D: elliptic2D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o elliptic2D elliptic2D.cpp $(LIBS)

clean:
	rm -f transport1D schrodinger1D elliptic1D elliptic2D
