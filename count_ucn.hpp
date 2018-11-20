#ifndef COUNT_UCN_H
#define COUNT_UCN_H

#include <vector>
#include "ucn_gen_PCG.hpp"

typedef struct coinc {
    double t;
    double dt;
    
//    bool operator<(const evt& rhs) const {
//        return t < rhs.t;
//    }
} coinc;

std::vector<coinc> countUCN_nopup(std::vector<evt> &events, double initialWindow, double telescopeWindow, int phCut);

std::vector<coinc> countUCN_pup(std::vector<evt> &events, double initialWindow, double telescopeWindow, int phCut);

double sumCoincs(std::vector<coinc> coincs, double binsize);

#endif /*COUNT_UCN_H*/