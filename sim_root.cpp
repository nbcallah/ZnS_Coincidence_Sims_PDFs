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
//        int rank = omp_get_thread_num();
//        int num = omp_get_num_threads();
    
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
    
    double rate = 1000;
    double length = 100000/rate;

    std::vector<double> t0s;
    double t = nextExp(r)/rate;
    while(t < length) {
        t0s.push_back(t);
        t += nextExp(r)/rate;
    }
    
    std::vector<evt> raw = generator.gen_evts(r, t0s);
    
    TTree* tree = make_root_tree(raw);
    
    tree->SaveAs("test.root");
    
//    
//
//    for(auto it = rates.begin(); it < rates.end(); it++) {
//        double rate = *it;
//        double length = 100000/rate;
//        std::vector<double> effs;
//
//        double err = 1;
//        double mean = 1;
//        int i = 0;
//        while(effs.size() < 10 || err > 1e-4*mean) {
//            std::vector<double> t0s;
//            double t = nextExp(r)/rate;
//            while(t < length) {
//                t0s.push_back(t);
//                t += nextExp(r)/rate;
//            }
//
//            std::vector<evt> raw = generator.gen_evts(r, t0s);
//            std::vector<coinc> coincs = countUCN(raw, 100e-9, 1000e-9, 8);
//
//            effs.push_back((double)coincs.size()/(double)t0s.size());
//            i++;
//
//            mean = std::accumulate(effs.begin(), effs.end(), 0.0)/effs.size();
//            err = sqrt(
//                std::accumulate(effs.begin(), effs.end(),
//                                0.0,
//                                [mean](double acc, double val){return acc + (mean-val)*(mean-val);})/(effs.size()-1)
//                )/sqrt(effs.size());
//
//            printf("%d | %d Eff: %f  (%lu/%lu); Avg: %f +/- %e\n", rank, i, (double)coincs.size()/(double)t0s.size(), coincs.size(), t0s.size(), mean, err);
//        }
//
//        printf("\n%f %f %e\n\n", rate, mean, err);
//    }
    
    ierr = MPI_Finalize();

    return 0;
}