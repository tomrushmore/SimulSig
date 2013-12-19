//
//  FileWriteR.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 15/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__FileWriteR__
#define __FFT_Chain__FileWriteR__

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "PrintfDebug.h"
#include "SetupMacros.h"
#include "GlobalTypedefs.h"

enum WRITER{
    WRITEFILE = 0,
    NWRITEFILE,
};

class FileWriteR{
public:
    FileWriteR();
    ~FileWriteR();
    void SetFilePath(const char* path);
    void DBWriteSetup(const char *path);
    void DBWriteFFTMFCCSplit(const std::vector<double>* fft,const std::vector<double>* mfcc,const std::vector<long>* frames_length_vec,const int mfcc_num_coef,const int fft_vals_p_frame,const long int num_frames,const int file_num,const long num_tracks,const  WRITER type,const std::vector<long>* buf_sizes);
    void DBWriteSingle(const std::vector<double>* frames);
    std::string WriteVector(const double * vector,const u_short_int *length,const char * name);
    void WriteVector(const double * vector,const int length);
  
  
    void WriteMelBanks(const double** filterbank,const int bank_length,const int num_bands);
    void WriteFFT(double * fft);
    void DBWriteFreqGraph(std::vector<double>* freqs,const int type);
    void DBWriteMarkov(const std::vector<double>*states);
    
private:
    char * file_name;
    const char * directory;
    int iter;
    
    std::ofstream *out_file;

    
};

#endif /* defined(__FFT_Chain__FileWriteR__) */
