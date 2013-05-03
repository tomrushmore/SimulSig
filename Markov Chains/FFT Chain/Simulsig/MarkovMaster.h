//
//  MarkovMaster.h
//  Simulsig
//
//  Created by Thomas Rushmore on 27/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __Simulsig__MarkovMaster__
#define __Simulsig__MarkovMaster__

#include <iostream>
#include <vector>
#include "MarkovChain.h"
#include "StateMFCCs.h"
#include "Database.h"
#include "MarkovPlayer.h"


class MarkovMaster{
public:
    MarkovMaster();
    ~MarkovMaster();
    void SetStates();
    void Train();
    void PlayChain();
    void CalcHoppedTrackLen(int num_tracks);
    
private:
    MarkovChain mkv_ch;
    MarkovPlayer mkv_pl;
    StateMFCCs state_mfccs;
    Database database;
    double** mfcc_train_matrix;
    double** fft_matrix;

    std::vector<int> trk_n_frame_hop;
    int max_h_len;
    
};

#endif /* defined(__Simulsig__MarkovMaster__) */
