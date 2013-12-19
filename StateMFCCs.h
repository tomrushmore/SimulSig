//
//  StateMFCCs.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 12/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__StateMFCCs__
#define __FFT_Chain__StateMFCCs__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "PrintfDebug.h"
#include "SetupMacros.h"

class StateMFCCs{
public:
    StateMFCCs();
    ~StateMFCCs();
    void SetupWrite(int s_num_states,int s_mfcc_num_coef, int s_fft_length,int s_num_tracks);
    void SetMFCCFFTPairs(double** mfccs,double** fft_matrix,int* best_fit_indices);
    void WriteMFCCFFTPairs();
    
    void ReadMFCCFFTStatePairs();
    double** GetMFCCMatrix();
    double** GetFFTMatrix();
    int GetNumMFCCCoef();
    int GetNumActiveClusters();
    int GetFFTLength();
    int GetNumTracks();
    int GetGrouping();
    int GetHopOverlap();
    int GetLongestTrackLength();
private:
    int num_states;
    int mfcc_num_coef;
    int fft_length;
    double** mfcc_at_state;
    double** fft_at_state;
    int num_tracks;
    int num_active_clusters;
    int grouping;
    int orig_mfcc_num_coef;
    int hop_overlap;
    int* track_frame_length;
};




#endif /* defined(__FFT_Chain__StateMFCCs__) */
