//
//  Concat.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 12/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__Player__
#define __FFT_Chain__Player__

#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>

#include "Database.h"
#include "fftw3.h"
#include "HammingWindow.h"
#include "PrintfDebug.h"
#include "SetupMacros.h"

class Player{
public:
    Player();
    ~Player();
    void Setup(int fft_size,int hop_size);
    void Play();
    
private:
    double EuclidDistance(double* a,double* b);
    fftw_complex* IntlvdToFFTW(double* interleaved_fft);
    fftw_complex* VecToFFTW();
    double* InverseFFTW();
    void WriteSignal();
    void FFTFramePushback(double* fft_values);
    void RealPushback();
    void RealHopOvrAddPushback();
    
    
    Database database;
    Database database2;
    std::vector<std::pair<double, double>> frames;
    std::vector<double> real_frames;
    fftw_complex *single_frame;
    fftw_complex *all_frames;
    double* fft_interleaved;
    double* real_signal;
    int hop_overlap;
    long vector_index;
    long vector_index_offset;
    fftw_plan ifft;
    
    HammingWindow ham;
};
#endif /* defined(__FFT_Chain__Player__) */
