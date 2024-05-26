CXX ?= g++

ifeq ($(DEBUG),1)
CXXFLAGS = -O0 -g
else
CXXFLAGS = -O3
endif
CXXFLAGS += -fopenmp

ifdef ARMA
INCPATH = -I. -I$(ARMA)/include
else
ifeq (,$(filter clean,$(MAKECMDGOALS)))
$(warning The path to Armadillo's library was not provided. I will look for it in the standard system directories.)
endif
endif

ifdef EIGEN
CXXFLAGS += -DEIGEN
INCPATH += -I$(EIGEN)
endif

ifndef SHARED_LIB
all: libmole.a
else
CXXFLAGS += -fPIC
all: libmole.so
endif

libmole.so: divergence.o gradient.o laplacian.o interpol.o robinbc.o utils.o
	$(CXX) -shared -o libmole.so *.o

libmole.a: divergence.o gradient.o laplacian.o interpol.o robinbc.o utils.o
	ar cqs libmole.a *.o

divergence.o: divergence.cpp divergence.h utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o divergence.o divergence.cpp

gradient.o: gradient.cpp gradient.h utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o gradient.o gradient.cpp

laplacian.o: laplacian.cpp laplacian.h divergence.h gradient.h utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o laplacian.o laplacian.cpp

interpol.o: interpol.cpp interpol.h utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o interpol.o interpol.cpp

robinbc.o: robinbc.cpp robinbc.h utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o robinbc.o robinbc.cpp

utils.o: utils.cpp utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o utils.o utils.cpp

clean:
	rm -f *.o libmole.a libmole.so
