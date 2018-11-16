#include "ucn_gen_PCG.hpp"
#include "pcg/pcg_random.hpp"
#include "rand_distributions.hpp"
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <limits>
#include <math.h>

ucn_gen_PCG::ucn_gen_PCG() {
    double mu1 = 0.0;
    double sigma1 = 0.0;
    double mu2 = 0.0;
    double sigma2 = 0.0;
    double p_pmt1 = 0.0;
    double p_shrt = 0.0;
    double t_shrt = 0.0;
    double p_med = 0.0;
    double t_med = 0.0;
    double t_long = 0.0;
}

ucn_gen_PCG::ucn_gen_PCG(double mu1,
               double sigma1,
               double mu2,
               double sigma2,
               double p_pmt1,
               double p_shrt,
               double t_shrt,
               double p_med,
               double t_med,
               double t_long) {
    this->mu1 = mu1;
    this->sigma1 = sigma1;
    this->mu2 = mu2;
    this->sigma2 = sigma2;
    this->p_pmt1 = p_pmt1;
    this->p_pmt2 = 1.0 - p_pmt1;
    this->p_shrt = p_shrt;
    this->t_shrt = t_shrt;
    this->p_med = p_med;
    this->t_med = t_med;
    this->p_long = 1.0 - p_shrt - p_med;
    this->t_long = t_long;
}

std::vector<evt> ucn_gen_PCG::gen_evts(pcg64 &r, std::vector<double> t0s) {
    unsigned long numEvts = t0s.size();
    unsigned long numPhotons = 0;
    std::vector<unsigned int> phs(numEvts);
    
    for(int i = 0; i < numEvts-1; i+= 2) {
        double z1, z2;
        next2Norm(r, &z1, &z2);
        uint64_t u = r(); //Probability of Li or Alpha is 50/50, use 2 bits
        
        if(u & (uint64_t(1) << 63)) { //mu1,sigma1
            phs[i] = round(exp(z1*sigma1 + mu1));
            numPhotons += phs[i];
        }
        else { //mu2,sigma2
            phs[i] = round(exp(z1*sigma2 + mu2));
            numPhotons += phs[i];
        }
        
        if(u & (uint64_t(1) << 62)) { //mu1,sigma1
            phs[i+1] = round(exp(z2*sigma1 + mu1));
            numPhotons += phs[i+1];
        }
        else { //mu2,sigma2
            phs[i+1] = round(exp(z2*sigma2 + mu2));
            numPhotons += phs[i+1];
        }
    }
    if((numEvts % 2) == 1) {
        double z1, z2;
        next2Norm(r, &z1, &z2);
        uint64_t u = r();
        
        if(u & (uint64_t(1) << 63)) { //mu1,sigma1
            phs[numEvts-1] = round(exp(z1*sigma1 + mu1));
            numPhotons += phs[numEvts-1];
        }
        else { //mu2,sigma2
            phs[numEvts-1] = round(exp(z1*sigma2 + mu2));
            numPhotons += phs[numEvts-1];
        }
    }
    
    std::vector<evt> evts;
    evts.reserve(numPhotons);
    unsigned long nShrt = 0;
    for(int i = 0; i < numEvts; i++) {
        for(int j = 0; j < phs[i]; j++) {
            double u = nextUc01o(r);
            //Now we'll combine the probabilities of the 2 independent processes: pmt1 and tail.
            //so it's like p_pmt1*p_shrt, p_pmt1*p_med, p_pmt1*p_long,
            //p_pmt2*p_shrt, p_pmt2*p_med, p_pmt2*p_long. ch can be checked simply. Then
            //split by channel for determining tau.
            int ch = u < p_pmt1 ? 1 : 2;
            double tau = u < p_pmt1 ?
                (u < p_pmt1*p_shrt ? t_shrt : (u < (p_pmt1*p_shrt + p_pmt1*p_med) ? t_med : t_long))
                 :
                (u < p_pmt1 + p_pmt2*p_shrt ? t_shrt : (u < (p_pmt1 + p_pmt2*p_shrt + p_pmt2*p_med) ? t_med : t_long));
            if(tau == t_shrt) {
                nShrt += 1;
            }
            double expOff = nextExp(r)*tau;
            evts.push_back({ch,
                            i,
                            t0s[i] + expOff});
        }
    }
    
    std::sort(evts.begin(), evts.end());
    
    unsigned long numCleaned = evts.size();
    
    for(int i = 0; i < evts.size(); i++) {
        if(isinf(evts[i].t)) {
            continue;
        }
        int j = i+1;
        while(j < evts.size() && (evts[j].t - evts[i].t) < 15e-9) {
            if(evts[j].ch == evts[i].ch) {
                evts[j].t = -std::numeric_limits<float>::infinity();
                numCleaned--;
            }
            j++;
        }
    }
    
//    printf("%lu %lu %f\n", numPhotons, numCleaned, (double)numCleaned/(double)numPhotons);
    
    std::vector<evt> cleanedEvts;
    cleanedEvts.reserve(numCleaned);
    for(int i = 0; i < evts.size(); i++) {
        if(isfinite(evts[i].t)) {
            cleanedEvts.push_back(evts[i]);
        }
    }
    
    return cleanedEvts;
}