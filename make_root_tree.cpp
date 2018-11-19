#include <vector>
#include <math.h>
#include "TTree.h"
#include "ucn_gen_PCG.hpp"
#include "make_root_tree.hpp"

//typedef struct root_entry {
//    int channel;
//    int edge;
//    int tag;
//    int full;
//    unsigned long time;
//    double realtime;
//} root_entry;

TTree* make_root_tree(std::vector<evt> &evts) {
    TTree* t = new TTree;
    root_entry event;
    event.tag = 0;
    event.full = 0;
    
    t->SetName("tmcs_0");
    t->SetTitle("Nates_Fake_Data");
    
    t->Branch("channel", (void *)&event.channel, "channel/I");
    t->Branch("edge", (void *)&event.edge, "edge/I");
    t->Branch("tag", (void *)&event.tag, "tag/I");
    t->Branch("full", (void *)&event.full, "full/I");
    t->Branch("time", (void *)&event.time, "time/l");
    t->Branch("realtime", (void *)&event.realtime, "realtime/D");
    
    for(auto it = evts.begin(); it < evts.end(); it++) {
        event.channel = it->ch;
        event.edge = it->id;
        event.realtime = it->t;
        event.time = (unsigned long)round(it->t / 800e-12);
        t->Fill();
    }
    
    return t;
}