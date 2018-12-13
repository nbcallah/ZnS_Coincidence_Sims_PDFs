#include "ucn_gen_PCG.hpp"
#include "count_ucn.hpp"
#include "rand_distributions.hpp"
#include "pcg/pcg_random.hpp"
#include "seed_source.h"
#include <vector>
#include <numeric>
#include <math.h>
#include <mpi.h>
#include <getopt.h>
#include <cstring>

int main(int argc, char** argv) {
    int ierr = MPI_Init(&argc, &argv);
    int nproc;
    int rank;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    int c;

    double window = 0.0;
    double sumWindow = 0.0;
    int numph = 0;

    while (1) {
        static struct option long_options[] = {
            {"window", required_argument, 0, 'w'},
            {"sumwindow", required_argument, 0, 's'},
            {"numph", required_argument, 0, 'n'},
            {0, 0, 0, 0},
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "w:s:n:", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1) {
            break;
        }

        switch(c) {
            case 'w':
                window = atof(optarg);
                break;
            case 's':
                sumWindow = atof(optarg);
                break;
            case 'n':
                numph = atoi(optarg);
                break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default:
                exit(1);
        }
    }
    
    if(window == 0.0 || sumWindow == 0.0 || numph == 0 || window < 0.8 || sumWindow < 0.8) {
        fprintf(stderr, "Error! Usage: ./sim --window=initial_coinc_window --sumwindow=telescope_window --numph=photon_threshold [in ns]\n");
        exit(1);
    }

    //Set the seed with some random bits from seed_source.h. We only need
    //128 bits since the state of the PCG generator is 128 bits long. We can
    //Get different unique outputs by by offseting the indices by 2*i ((0,1), (2,3), (4,5), ...)
    __uint128_t seed = (((__uint128_t)seed_source[1]) << 127) | (((__uint128_t)seed_source[0]));
    //One nice thing about the PCG generator is that it has 2^127 streams,
    //which are each independent. Assign one stream per core based on core ID.
    pcg64 r(seed, (unsigned int)rank);

    //One low, one medium, one long is sufficient for drawing a straight line
    std::vector<double> rates = {100, 5000, 10000};

    ucn_gen_PCG generator = ucn_gen_PCG();

    //These are the parameter central values I've identified through the fits.
    generator.mu1 = 3.168308578;
    generator.sigma1 = 0.33346646;
    generator.mu2 = 3.77243909665;
    generator.sigma2 = 0.25119938;
    generator.p_pmt1 = 0.56596113354854816;
    generator.p_pmt2 = 1.0 - generator.p_pmt1;
    generator.p_shrt = 0.215304093372;
    generator.t_shrt = 153.902059371e-9;
    generator.p_med = 0.38085062563;
    generator.t_med = 1864.84601579e-9;
    generator.p_long = 1.0 - generator.p_shrt - generator.p_med;
    generator.t_long = 16311.0287305e-9;

    //Want to measure efficiency by rate. Loop over rates.
    for(auto it = rates.begin(); it < rates.end(); it++) {
        double rate = *it;
        double length = 100000/rate; //Same expected counts per run. Not strictly necessary.
        std::vector<double> effs;
        std::vector<double> effsDT;
        std::vector<double> effsnopup;
        std::vector<double> effsnopupDT;

        double err = 1;
        double mean = 1;
        int i = 0;
        while(effs.size() < 10 || err > 1e-4*mean) { //Measure to 1 part per 10 thousand
            //Generate flat rate with poisson process
            std::vector<double> t0s;
            double t = nextExp(r)/rate;
            while(t < length) {
                t0s.push_back(t);
                t += nextExp(r)/rate; //interarrival time is exponentially distributed
            }

            std::vector<evt> raw = generator.gen_evts(r, t0s);
            std::vector<coinc> coincs = countUCN_pup(raw, window*1e-9, sumWindow*1e-9, numph);
            std::vector<coinc> coincs_nopup = countUCN_nopup(raw, window*1e-9, sumWindow*1e-9, numph);

            //Efficiency is the number of UCN counted divided by the
            //number of UCN created.
            effs.push_back((double)coincs.size()/(double)t0s.size());
            effsDT.push_back(sumCoincs(coincs, 1.0)/(double)t0s.size());
            effsnopup.push_back((double)coincs_nopup.size()/(double)t0s.size());
            effsnopupDT.push_back(sumCoincs(coincs_nopup, 1.0)/(double)t0s.size());
            i++;

            //calculate mean, and std. dev. / sqrt(N) for the efficiencies
            mean = std::accumulate(effs.begin(), effs.end(), 0.0)/effs.size();
            err = sqrt(
                std::accumulate(effs.begin(), effs.end(),
                                0.0,
                                [mean](double acc, double val){return acc + (mean-val)*(mean-val);})/(effs.size()-1)
                )/sqrt(effs.size());

//            printf("%d | %d Eff: %f  (%lu/%lu); Avg: %f +/- %e\n", rank, i, (double)coincs.size()/(double)t0s.size(), coincs.size(), t0s.size(), mean, err);
        }

        mean = std::accumulate(effs.begin(), effs.end(), 0.0)/effs.size();
        err = sqrt(
            std::accumulate(effs.begin(), effs.end(),
                            0.0,
                            [mean](double acc, double val){return acc + (mean-val)*(mean-val);})/(effs.size()-1)
            )/sqrt(effs.size());
//        printf("(%f, %f, %e), #pup_nocorr\n", rate, mean, err);
        printf("%f %f %e 0 0 %f %f %d\n", rate, mean, err, window, sumWindow, numph);
        
        mean = std::accumulate(effsDT.begin(), effsDT.end(), 0.0)/effsDT.size();
        err = sqrt(
            std::accumulate(effsDT.begin(), effsDT.end(),
                            0.0,
                            [mean](double acc, double val){return acc + (mean-val)*(mean-val);})/(effsDT.size()-1)
            )/sqrt(effsDT.size());
//        printf("(%f, %f, %e), #pup_corr\n", rate, mean, err);
        printf("%f %f %e 0 1 %f %f %d\n", rate, mean, err, window, sumWindow, numph);
        
        mean = std::accumulate(effsnopup.begin(), effsnopup.end(), 0.0)/effsnopup.size();
        err = sqrt(
            std::accumulate(effsnopup.begin(), effsnopup.end(),
                            0.0,
                            [mean](double acc, double val){return acc + (mean-val)*(mean-val);})/(effsnopup.size()-1)
            )/sqrt(effsnopup.size());
//        printf("(%f, %f, %e), #nopup_nocorr\n", rate, mean, err);
        printf("%f %f %e 1 0 %f %f %d\n", rate, mean, err, window, sumWindow, numph);
        
        mean = std::accumulate(effsnopupDT.begin(), effsnopupDT.end(), 0.0)/effsnopupDT.size();
        err = sqrt(
            std::accumulate(effsnopupDT.begin(), effsnopupDT.end(),
                            0.0,
                            [mean](double acc, double val){return acc + (mean-val)*(mean-val);})/(effsnopupDT.size()-1)
            )/sqrt(effsnopupDT.size());
//        printf("(%f, %f, %e), #nopup_corr\n", rate, mean, err);
        printf("%f %f %e 1 1 %f %f %d\n", rate, mean, err, window, sumWindow, numph);
    }
    
    ierr = MPI_Finalize();

    return 0;
}