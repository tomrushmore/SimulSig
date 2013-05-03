//
//  DBWrite.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 22/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__DBWrite__
#define __FFT_Chain__DBWrite__

#include <iostream>
#include <vector>
#include "FFT.h"
#include "PrintfDebug.h"
#include "GlobalTypedefs.h"

class DBWrite{
public:
    DBWrite();
    ~DBWrite();
    void AddVectorMag(double* mag_vector);
    void AddVectorReal(double* buffer);
    void AddFFTMFCCPair(fftw_complex* fft,double *mfcc);
    void AddFFTMFCCSplit(fftw_complex* fft,double *mfcc);
    void AddNumFramesPerTrack(long num_frames);
    void SetFixedFramesLength(const u_short_int *s_fft_frm_len);
    void SetFixFFTMFCCLength(const u_short_int *s_fft_frm_len,const u_short_int *s_mfcc_frm_len);
    void AddMarkov(double states);
    void SetBufferSizes(long size);
    std::vector<long>* GetBufferSizes();

    std::vector<double>* GetMFCCVector();
    std::vector<double>* GetFFTVector();
    std::vector<long>*   GetTrackNumFramesVector();
    int GetFramesLength();
    unsigned long GetSize();
    std::vector<double>* GetMarkov();

private:
    std::vector<std::pair<double, double>> frames;
    std::vector<double> frames_single;
    std::vector<double> fft_interleaved_v;
    std::vector<double> mfcc_v;
    std::vector<double> fft_mfcc_seq_v;
    std::vector<long> num_frames_p_track;
    std::vector<double> Markovs;
    std::vector<long> buffer_sizes;
    u_short_int *fft_frm_len;
    u_short_int *mfcc_frm_len;
    
};

#endif /* defined(__FFT_Chain__DBWrite__) */
