//
//  DBWrite.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 22/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "DBWrite.h"

DBWrite::DBWrite()
{
}

DBWrite::~DBWrite()
{
    
}

void DBWrite::AddMarkov(double states)
{
    Markovs.push_back(states);
}

std::vector<double>* DBWrite::GetMarkov()
{
    return &Markovs;
}

void DBWrite::AddFFTMFCCPair(fftw_complex *fft, double *mfcc)
{
    for(int i = 0 ; i < (*fft_frm_len); i++){
        fft_mfcc_seq_v.push_back(fft[i][0]);
        fft_mfcc_seq_v.push_back(fft[i][1]);
    }
    for(int i = 0 ; i < (*mfcc_frm_len); i++){
        fft_mfcc_seq_v.push_back(mfcc[i]);
    }

}

void DBWrite::AddFFTMFCCSplit(fftw_complex *fft, double *mfcc)
{
    for(int i = 0; i < (*fft_frm_len); i++){
        fft_interleaved_v.push_back(fft[i][0]);
        fft_interleaved_v.push_back(fft[i][1]);
    }
    for(int i = 0; i < (*mfcc_frm_len); i++){
        mfcc_v.push_back(mfcc[i]);
    }
}

void DBWrite::AddNumFramesPerTrack(long num_frames)
{
    num_frames_p_track.push_back(num_frames);
}

void DBWrite::SetBufferSizes(long size)
{
    buffer_sizes.push_back(size);
}

std::vector<long>* DBWrite::GetBufferSizes()
{
    return &buffer_sizes;
}

void DBWrite::AddVectorMag(double *mag_vector)
{
    for(int i = 0 ; i < (*fft_frm_len); i++){
        frames_single.push_back(mag_vector[i]);
    }
}

void DBWrite::SetFixedFramesLength(const u_short_int *s_fft_frm_len)
{
    fft_frm_len = const_cast<u_short_int*>(s_fft_frm_len);
}

int DBWrite::GetFramesLength()
{
    return (*fft_frm_len);
}
void DBWrite::AddVectorReal(double *buffer)
{
    for(int i = 0 ; i < ((*fft_frm_len) - 1)*2; i++){
        frames_single.push_back(buffer[i]);
    }
}

void DBWrite::SetFixFFTMFCCLength(const u_short_int *s_fft_frm_len,const u_short_int *s_mfcc_frm_len)
{
    mfcc_frm_len = const_cast<u_short_int*>(s_mfcc_frm_len);
    fft_frm_len  = const_cast<u_short_int*>(s_fft_frm_len);
    printf("fft_frm_len set to %d\n",(*fft_frm_len));
}



std::vector<double>* DBWrite::GetMFCCVector()
{
    return &mfcc_v;
}

std::vector<double>* DBWrite::GetFFTVector()
{
    return &fft_interleaved_v;
}

std::vector<long>* DBWrite::GetTrackNumFramesVector()
{
    return &num_frames_p_track;
}

unsigned long DBWrite::GetSize()
{
    debug_print(("Vector Size : %ld",frames.size()));
    return frames.size();
}