#CC=gcc
#CXX=g++
CXX=mpic++
OPT=-O3
CFLAGS=$(OPT)
CPPFLAGS=$(OPT) -std=c++11
CPPLFLAGS=$(OMP) -lgsl -lgslcblas

all: sim

debug: OPT=-g
debug: all

sim: sim.o ucn_gen_PCG.o count_ucn.o rand_distributions.o
	$(CXX) -o sim sim.o ucn_gen_PCG.o count_ucn.o rand_distributions.o $(CPPLFLAGS)

sim.o: sim.cpp
	$(CXX) $(CPPFLAGS) -c -o sim.o sim.cpp
	
ucn_gen_PCG.o: ucn_gen_PCG.cpp ucn_gen_PCG.hpp
	$(CXX) $(CPPFLAGS) -c -o ucn_gen_PCG.o ucn_gen_PCG.cpp
    
count_ucn.o: count_ucn.cpp count_ucn.hpp
	$(CXX) $(CPPFLAGS) -c -o count_ucn.o count_ucn.cpp

rand_distributions.o: rand_distributions.cpp rand_distributions.hpp
	$(CXX) $(CPPFLAGS) -c -o rand_distributions.o rand_distributions.cpp
    
clean:
	rm -rf *.o
	rm sim