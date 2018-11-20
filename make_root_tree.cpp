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
    //Allocate new TTree
    TTree* t = new TTree;
    //Struct to hold each tree branch
    root_entry event;
    //We won't use tag or full
    event.tag = 0;
    event.full = 0;
    
    //Name the tree tmcs_0 as in Robby's trees
    t->SetName("tmcs_0");
    t->SetTitle("Nates_Fake_Data");
    
    //Set addresses
    t->Branch("channel", (void *)&event.channel, "channel/I");
    t->Branch("edge", (void *)&event.edge, "edge/I");
    t->Branch("tag", (void *)&event.tag, "tag/I");
    t->Branch("full", (void *)&event.full, "full/I");
    t->Branch("time", (void *)&event.time, "time/l");
    t->Branch("realtime", (void *)&event.realtime, "realtime/D");
    
    //Loop through synthetic data; fill tree for each entry
    for(auto it = evts.begin(); it < evts.end(); it++) {
        event.channel = it->ch;
        event.edge = it->id;    //edge holds a unique ID for each UCN event.
        event.realtime = it->t;
        event.time = (unsigned long)round(it->t / 800e-12);
        t->Fill();
    }
    
    return t;
}