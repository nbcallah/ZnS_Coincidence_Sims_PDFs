#include "count_ucn.hpp"
#include "ucn_gen_PCG.hpp"
#include "pcg/pcg_random.hpp"
#include "rand_distributions.hpp"
#include <vector>
#include <math.h>

std::vector<coinc> countUCN_nopup(std::vector<evt> &events, double initialWindow, double telescopeWindow, int phCut) {
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

std::vector<coinc> countUCN_pup(std::vector<evt> &events, double initialWindow, double telescopeWindow, int phCut) {
    std::vector<coinc> coincs;
    
    int i = 0;
    while(i < events.size()) {
        int j = i+1;
        while(j < events.size() && (events[j].t - events[i].t) < (initialWindow + 0.1e-9)) {
            if(events[j].ch != events[i].ch) {
                int numPh = 1;
                int k = i+1;
                while(k < events.size() && (events[k].t - events[k-1].t) < (telescopeWindow + 0.1e-9)) {
                    numPh++;
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

std::vector<coinc> countUCN_chris(std::vector<evt> &events, double initialWindow, double telescopeWindow, int phCut, pcg64 &r) {
    double avgPhotonRate = events.size()/(events.back().t - events.front().t);
    
    std::vector<coinc> coincs;
    
    int i = 0;
    while(i < events.size()) {
        int j = i+1;
        while(j < events.size() && (events[j].t - events[i].t) < (initialWindow + 0.1e-9)) {
            if(events[j].ch != events[i].ch) {
                int numPh = 1;
                int k = i+1;
                while(k < events.size() && (events[k].t - events[k-1].t) < (telescopeWindow + 0.1e-9)) {
                    numPh++;
                    k++;
                }
                //We'll cheat a little for Chris's method. Since we're only passing
                //uniform rates, just take the total number of PMT events over the length
                //for the photon rate. Then roll the dice on whether or not to accept an event
                //based on the expected counts.
                double expectedCounts = avgPhotonRate*(events[k-1].t-events[i].t);
//                if(expectedCounts > 1.5) { printf("%f\n", expectedCounts); }
                double phCutPrime = phCut + expectedCounts;
                if(
                    (numPh >= ceil(phCutPrime))
                    ||
                    ((numPh >= floor(phCutPrime)) && (nextUc01o(r) >= phCutPrime - floor(phCutPrime)))
                ) {
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

double sumCoincs(std::vector<coinc> coincs, double binsize) {
    double dtCounts = 0;
    double t = floor(coincs.front().t/binsize);
    double binContent = 0;
    double bindt = 0;
    for(auto it = coincs.begin(); it < coincs.end(); it++) {
        if(floor(it->t/binsize) > t) { //new bin
            t = floor(it->t/binsize);
            dtCounts += binContent/((binsize - bindt)/binsize);
            binContent = 0;
            bindt = 0;
        }
        bindt += it->dt;
        binContent += 1;
    }
    
    dtCounts += binContent/((binsize - bindt)/binsize);
    
    return dtCounts;
}