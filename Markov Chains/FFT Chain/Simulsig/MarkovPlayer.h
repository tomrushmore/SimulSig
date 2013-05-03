//
//  MarkovPlayer.h
//  Simulsig
//
//  Created by Thomas Rushmore on 02/04/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __Simulsig__MarkovPlayer__
#define __Simulsig__MarkovPlayer__

#include <iostream>
#include <ostream>

#include <vector>
#include "MarkovChain.h"
#include "fftw3.h"
#include <sndfile.hh>

#define RAMP 1

enum MKVTYPE{
    HOMOG = 0,
    NHOMOG
};

class MarkovPlayer{
public:
    MarkovPlayer();
    ~MarkovPlayer();
    void Setup(MarkovChain* s_mkv_ch,double** s_fft_matrix,int s_n_sequences,int s_max_track_length,int s_n_coef,int s_n_active_clusters,int s_fft_length,int s_grouping,int s_hop_overlap);
    void GenerateChainSequences(const MKVTYPE type);
    void GenerateSampleBuffers();
    void GenerateSignals();
    void WriteSignal();
    void WriteFinalSignals();
    void RealHopOvrAddPushback(int sample_num);
private:
    
    void BuildRamp();
    fftw_complex* IntlvdToFFTW(double *interleaved_fft);
    double* InverseFFTW();
    fftw_plan ifft;
    double* real_signal;
    std::vector<double> real_frames;
    int hop_overlap;
    int rem_frames;
    double* forward_ramp;
    
    std::vector<std::vector<int>> chain_plays;
    std::vector<int> single_chain_play;
    int n_sequences;
    int max_track_length;
    
    MarkovChain* mkv_ch;
    double** fft_matrix;
    int n_coef;
    int n_active_clusters;
    int fft_length;
    int fft_length_2_m1;
    int fft_length_m2;
    int grouping;
    fftw_complex *single_frame;
    double* single_frame_interleaved;
    double** time_dom_signals;
    int vector_index;
    std::vector<std::vector<double>> all_time_doms;
    std::vector<std::vector<double>> final_signals;
    std::vector<std::vector<double>> abs_signal;
    std::vector<double> cur_track;
    

};

#endif /* defined(__Simulsig__MarkovPlayer__) */
