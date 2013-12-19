//
//  FeatureExtraction.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 15/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__FeatureExtraction__
#define __FFT_Chain__FeatureExtraction__

#include <limits>

#include "GlobalTypedefs.h"
#include "HammingWindow.h"
#include "AudioFileReader.h"
#include "FFT.h"
#include "FileWriteR.h"
#include "DBWrite.h"
#include "MFCC.h"
#include "PrintfDebug.h"
#include "SetupMacros.h"

#define WRITE PAIR

enum PROCESSTYPE{
    CLUSTANALYSIS,
    MARKOVANALYSIS
};

class FeatureExtraction{
public:
    FeatureExtraction();
    ~FeatureExtraction();
    bool Setup(const int s_fft_size,const int s_hop_overlap,const int n_mel_bands,const int sample_rate,const double min_freq,const double max_freq,const int samp_chunk_size,const int mfcc_grouping);
    void SetWriteFilePath(const char * write_file_path);
    void ProcessFile(const char * file_path);
    void ProcessSingleFile(const char *file_path,const PROCESSTYPE type);
    void ProcessDatabase(const char *file_path,int file_num);
    void WriteDatabase(long num_tracks,WRITER type);
    u_short_int GetFinMBLen();
    
private:
    void ConsoleDebug();
    static bool IsPowerTwo(int number);
    
    AudioFileReader file_reader;
    HammingWindow hamming_window;
    FFT fft;
    MFCC mfcc;
    FileWriteR file_writer;
    DBWrite database;
    
    u_short_int fft_len;
    u_short_int hop_ovrlp;
    u_short_int n_mfcc_coefs;
    u_short_int fin_n_mfcc_coefs;
    long int total_num_frames;
    char * directory;
    std::string hamming_directory;
    fftw_complex *fft_write;
    fftw_complex *fft_write_no_w;
    double *win_chunk;
    double *n_win_chunk;
    int sample_chunk_size;
    int sample_chunk_offset;
};

enum WRITETYPE{
    MAGN,
    MFCC,
    PAIR
};




#endif /* defined(__FFT_Chain__FeatureExtraction__) */
