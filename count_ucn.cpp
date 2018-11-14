#include "count_ucn.hpp"
#include "ucn_gen_PCG.hpp"
#include <vector>
#include <math.h>

std::vector<coinc> countUCN(std::vector<evt> &events, double initialWindow, double telescopeWindow, int phCut) {
    std::vector<coinc> coincs;
    
    int i = 0;
    while(i < events.size()) {
        int j = i+1;
        while(j < events.size() && (events[j].t - events[i].t) < (initialWindow + 0.1e-9)) {
            if(events[j].ch != events[i].ch && events[j].id == events[i].id) {
                int numPh = 1;
                int k = i+1;
                while(k < events.size() && (events[k].t - events[k-1].t) < (telescopeWindow + 0.1e-9)) {
                    if(events[k].id == events[i].id) {
                        numPh++;
                    }
                    k++;
                }
                if(numPh >= phCut) {
                    coincs.push_back({events[i].t, events[k-1].t-events[i].t + (telescopeWindow + 0.1e-9)});
                    i = k-1;
                    break;
                }
            }
            j++;
        }
        i++;
    }
    
    return coincs;
}

//double sumCoincs(std::vector<coinc> coincs) {
//    double first = floor(coincs.front());
//    double last = ceil(coincs.back());
//}