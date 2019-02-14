#ifndef UCN_GEN_PCG_H
#define UCN_GEN_PCG_H

#include <vector>
#include "pcg/pcg_random.hpp"

typedef struct evt { //light-weight struct for PMT events
    int ch;
    int id; //unique ID for each event.
    double t;
    
    bool operator<(const evt& rhs) const {
        return t < rhs.t;
    }
} evt;

class ucn_gen_PCG {
    public:
        double mu1;
        double sigma1;
        double mu2;
        double sigma2;
        double p_pmt1;
        double p_pmt2;
        double p_shrt;
        double t_shrt;
        double p_med;
        double t_med;
        double p_long;
        double t_long;
        double t_trunc;
        
        std::vector<evt> gen_evts(pcg64 &r, std::vector<double> t0s);
    
        ucn_gen_PCG();
        ucn_gen_PCG(double mu1,
               double sigma1,
               double mu2,
               double sigma2,
               double p_pmt1,
               double p_shrt,
               double t_shrt,
               double p_med,
               double t_med,
               double t_long,
               double t_trunc);
};

#endif /*UCN_GEN_PCG_H*/