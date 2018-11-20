#include "ucn_gen_PCG.hpp"
#include "count_ucn.hpp"
#include "rand_distributions.hpp"
#include "make_root_tree.hpp"
#include "pcg/pcg_random.hpp"
#include "seed_source.h"
#include <vector>
#include <numeric>
#include <math.h>
#include <mpi.h>
#include "TTree.h"

int main(int argc, char** argv) {
    int ierr = MPI_Init(&argc, &argv);
    int nproc;
    int rank;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    __uint128_t seed = (((__uint128_t)seed_source[1]) << 127) | (((__uint128_t)seed_source[0]));
    //One nice thing about the PCG generator is that it has 2^127 streams,
    //which are each independent. Assign one stream per core based on core ID.
    pcg64 r(seed, (unsigned int)rank);

    std::vector<double> rates = {100, 5000, 10000};

    ucn_gen_PCG generator = ucn_gen_PCG();

    //These are the parameter central values I've identified through the fits.
    generator.mu1 = 3.168308578;
    generator.sigma1 = 0.33346646;
    generator.mu2 = 3.77243909665;
    generator.sigma2 = 0.25119938;
    generator.p_pmt1 = 0.566725424184;
    generator.p_pmt2 = 1.0 - generator.p_pmt1;
    generator.p_shrt = 0.215304093372;
    generator.t_shrt = 153.902059371e-9;
    generator.p_med = 0.38085062563;
    generator.t_med = 1864.84601579e-9;
    generator.p_long = 1.0 - generator.p_shrt - generator.p_med;
    generator.t_long = 16311.0287305e-9;
    
    
    //Generate a uniform rate by simulating a poisson process.
    double rate = 1000;
    double length = 100000/rate;
    
    std::vector<double> t0s;
    double t = nextExp(r)/rate;
    while(t < length) {
        t0s.push_back(t);
        t += nextExp(r)/rate; //interarrival time is exponentially distributed
    }
    
    std::vector<evt> raw = generator.gen_evts(r, t0s);
    
    TTree* tree = make_root_tree(raw);
    
    printf("%lld\n", tree->GetEntries());
    
    ierr = MPI_Finalize();

    return 0;
}