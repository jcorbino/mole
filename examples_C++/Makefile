MOLE = ../mole_C++

CXX ?= g++

ifeq ($(DEBUG),1)
CXXFLAGS = -O0 -g
else
CXXFLAGS = -O3
endif
CXXFLAGS += -fopenmp -DARMA_DONT_USE_WRAPPER

INCPATH = -I. $(if $(ARMA), -I$(ARMA)/include) -I$(MOLE)

#LIBS = -L$(MOLE) -lmole $(if $(ARMA),-L$(ARMA)/lib -Wl,-rpath,$(ARMA)/lib) -larmadillo -lopenblas
LIBS = -L$(MOLE) -lmole

ifdef ARMA
LIBS += -L$(ARMA)/lib -Wl,-rpath,$(ARMA)/lib -larmadillo -lopenblas
else
$(error ARMA is not defined. Please define ARMA with the path to the Armadillo library.)
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
