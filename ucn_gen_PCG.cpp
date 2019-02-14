#include "ucn_gen_PCG.hpp"
#include "pcg/pcg_random.hpp"
#include "rand_distributions.hpp"
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <limits>
#include <math.h>

ucn_gen_PCG::ucn_gen_PCG() { //Default constructer, just 0 everything.
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
    double t_trunc = 0.0;
}

ucn_gen_PCG::ucn_gen_PCG(double mu1, //Construct by passing each parameter
               double sigma1,
               double mu2,
               double sigma2,
               double p_pmt1,
               double p_shrt,
               double t_shrt,
               double p_med,
               double t_med,
               double t_long,
               double t_trunc) {
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
    this->t_trunc = t_trunc;
}

std::vector<evt> ucn_gen_PCG::gen_evts(pcg64 &r, std::vector<double> t0s) {
    unsigned long numEvts = t0s.size();
    unsigned long numPhotons = 0;
    std::vector<unsigned int> phs(numEvts);
    
    //First task is to generate pulse heights for each UCN and tally number of
    //photons.
    for(int i = 0; i < numEvts-1; i+= 2) { //We get normal RVs 2 at a time
        double z1, z2;
        next2Norm(r, &z1, &z2);
        uint64_t u = r(); //Probability of Li or Alpha is 50/50, use 2 highest bits
        
        if(u & (uint64_t(1) << 63)) { //mu1,sigma1
            phs[i] = round(exp(z1*sigma1 + mu1)); //round to simulate discreteness of photons
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
    if((numEvts % 2) == 1) { //If we had odd number of events generate last one
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
    
    //Now that we've generate the pulse heights, we will generate photons from the ZnS
    //scintillation curve.
    std::vector<evt> evts;
    evts.reserve(numPhotons); //will need numPhoton entries
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
            double expOff = nextExp(r)*tau; //offset from beginning of event is exponentially
                                            //distributed with t_shrt, t_med, or t_long
            //Push back event, with unique ID of the iterator through the PHS array.
            if(expOff < t_trunc) {
                evts.push_back({ch,
                                i,
                                t0s[i] + expOff});
            }
        }
    }
    
    //time-order events
    std::sort(evts.begin(), evts.end());
    
    //Now we need to perform a search for events which would be unseen due to deadtime.
    unsigned long numCleaned = evts.size();
    
    for(int i = 0; i < evts.size(); i++) {
        if(isinf(evts[i].t)) { //Events which have been already deadtime'd don't generate more deadtime
            continue;
        }
        int j = i+1;
        while(j < evts.size() && (evts[j].t - evts[i].t) < 15e-9) { //While we have events which could be deadtime'd
            if(evts[j].ch == evts[i].ch) { //If the channels match, then event would not be seen.
                evts[j].t = -std::numeric_limits<float>::infinity(); //Signal is -infinite time.
                numCleaned--;
            }
            j++;
        }
    }
        
    //Now we just populate a new vector with only events which would be seen by the DAQ.
    std::vector<evt> cleanedEvts;
    cleanedEvts.reserve(numCleaned);
    for(int i = 0; i < evts.size(); i++) {
        if(isfinite(evts[i].t)) {
            cleanedEvts.push_back(evts[i]);
        }
    }
    
    return cleanedEvts;
}