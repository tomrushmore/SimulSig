//
//  DBReader.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 09/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__Database__
#define __FFT_Chain__Database__

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include "PrintfDebug.h"
#include "SetupMacros.h"

class Database{
public:
    Database();
    ~Database();
    void SetFileDirectory(const char* dir);
    int ReadFile(int file_num);
    int ReadFileGroup(int file_num,int s_grouping,int type);
    double** GetMFCCMatrix();
    double** GetFFTMatrix();
    int GetTrackFrameLengths(int index);
    int GetLongestTrackLength();
    long GetNumFrames();
    int  GetNumMFCCCoefs();
    int  GetFFTLength();
    int  GetFFTLengthGrouped();
    int GetNumTracks();
    void GroupTrackLengths();
    int GetGroupedS();
    long GetBufSizes(int idx);
    
private:
    long grouped_size;
    double** mfcc_matrix;
    double** fft_matrix;
    double** fft_matrix_two;
    int* track_frame_length;
    long* track_buffer_sizes;
    int fft_frame_length;
    int mfcc_frame_r_length;
    long int database_frame_num;
    int mfcc_num_coef;
    long int mfcc_total_num;
    int num_tracks;
    int grouping;
    int fft_length_grouped;
    int gr_size;
    int gr_coef;
    int gr_fft;
   // std::string directory;
    const std::string fft_path  = "fft_frames.txt";
    const std::string mfcc_path = "mfcc_frames.txt";
    
};
#endif /* defined(__FFT_Chain__Database__) */
