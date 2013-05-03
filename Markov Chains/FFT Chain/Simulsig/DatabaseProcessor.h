//
//  DatabaseProcessor.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 26/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__DatabaseProcessor__
#define __FFT_Chain__DatabaseProcessor__

#include <iostream>
#include <string>
#include <vector>
#include <dirent.h>

#include "FeatureExtraction.h"
#include "PrintfDebug.h"
#include "SetupMacros.h"

class DatabaseProcessor{
public:
    DatabaseProcessor(const char* s_samp_directory,const char* s_folder,const char* s_main_directory);
    ~DatabaseProcessor();
    int FeExSetup(const int set_analysis_size,const int set_hop_overlap,const int n_mel_bands,const int sample_rate,const double min_freq,const double max_freq,const int samp_chunk_size,const int mfcc_grouping);
    void Process(WRITER type);
    
private:
    int GetFileNames();

    std::string samp_directory,main_directory;
    const char* folder;
    FeatureExtraction FeExtract;
    std::vector<std::string> file_names;
    
};
#endif /* defined(__FFT_Chain__DatabaseProcessor__) */
