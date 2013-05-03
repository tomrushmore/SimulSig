//
//  MFCC.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 28/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__MFCC__
#define __FFT_Chain__MFCC__

#include <iostream>
#include <math.h>
#include <vector>
#include <Accelerate/Accelerate.h>

#include "fftw3.h"
#include "PrintfDebug.h"
#include "GlobalTypedefs.h"

class MFCC{
public:
    MFCC();
    ~MFCC();
    void Setup(const u_short_int *s_n_bands,const u_short_int s_mn_freq,const u_short_int s_mx_freq,const u_short_int *s_fft_len,const int s_sampr);
    double*  CalcMFCC(double * buffer);
    double** GetFilterBank();
    int      GetNumBands();
    u_short_int  GetVecLength();
    u_short_int* GetVelLengthAdr();
    double*  GetMF();
    std::vector<double>* LinFreqGraph();
    std::vector<double>* MelFreqGraph();
    
private:
    void BuildFilter(double* result,double* filter_bank,const int filter_index,const u_short_int size);
    static double FtoMel(double freq);
    static double MtoFreq(double freq);
    void    BuildFilterBank();

    u_short_int *n_mel_bands;
    u_short_int *fft_len;
    u_short_int fft_len_d2;
    u_short_int fft_len_d2_p1;
    u_short_int n_mel_b_d2_p1;
    int      sampr;
    double*  mel_bank_c_freqs;
    double** mel_bank;
    double** mel_applied_filter;
    double*  mfcc_temp_vec;
    double*  final_output;
    u_short_int   mn_freq;
    u_short_int   mx_freq;
    std::vector<double> lin_freq_graph;
    std::vector<double> mel_freq_graph;
    fftw_plan dct_p;
};

#endif /* defined(__FFT_Chain__MFCC__) */
