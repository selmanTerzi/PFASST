PFASST   := $(abspath ../..)
STACK    ?= $(PFASST)/stack
CXX      := $(STACK)/bin/mpic++
CXXFLAGS ?= -std=c++11 -I$(PFASST)/include -I$(STACK)/include -I$(STACK)/include/eigen3

LDFLAGS  ?= -L$(STACK)/lib -Wl,-rpath=$(STACK)/lib -lboost_program_options -lfftw3

all: parareal_test

parareal_test: parareal_test.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)
