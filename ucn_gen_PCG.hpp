#ifndef UCN_GEN_PCG_H
#define UCN_GEN_PCG_H

#include <vector>
#include "pcg/pcg_random.hpp"

typedef struct evt {
    int ch;
    int id;
    double t;
    
    bool operator<(const evt& rhs) const {
        return t < rhs.t;
    }
} evt;

std::vector<evt> gen_evts(pcg64 &r, std::vector<double> t0s);

#endif /*UCN_GEN_PCG_H*/