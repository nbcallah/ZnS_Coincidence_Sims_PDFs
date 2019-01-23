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

sim_root: sim_root.o ucn_gen_PCG.o count_ucn.o rand_distributions.o make_root_tree.o
	$(CXX) -o sim_root sim_root.o ucn_gen_PCG.o count_ucn.o rand_distributions.o make_root_tree.o $(CPPLFLAGS) `root-config --libs`
    
sim_minimal: sim_minimal.o ucn_gen_PCG.o count_ucn.o rand_distributions.o
	$(CXX) -o sim_minimal sim_minimal.o ucn_gen_PCG.o count_ucn.o rand_distributions.o $(CPPLFLAGS)

sim.o: sim.cpp
	$(CXX) $(CPPFLAGS) -c -o sim.o sim.cpp

sim_root.o: sim_root.cpp
	$(CXX) $(CPPFLAGS) `root-config --cflags` -c -o sim_root.o sim_root.cpp
    
sim_minimal.o: sim_minimal.cpp
	$(CXX) $(CPPFLAGS) -c -o sim_minimal.o sim_minimal.cpp
	
ucn_gen_PCG.o: ucn_gen_PCG.cpp ucn_gen_PCG.hpp
	$(CXX) $(CPPFLAGS) -c -o ucn_gen_PCG.o ucn_gen_PCG.cpp
    
count_ucn.o: count_ucn.cpp count_ucn.hpp
	$(CXX) $(CPPFLAGS) -c -o count_ucn.o count_ucn.cpp

rand_distributions.o: rand_distributions.cpp rand_distributions.hpp
	$(CXX) $(CPPFLAGS) -c -o rand_distributions.o rand_distributions.cpp

make_root_tree.o: make_root_tree.cpp make_root_tree.hpp
	$(CXX) $(CPPFLAGS) `root-config --cflags` -c -o make_root_tree.o make_root_tree.cpp
    
clean:
	rm -rf *.o
	rm -rf sim
	rm -rf sim_root
	rm -rf sim_minimal