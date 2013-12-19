//
//  DatabaseProcessor.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 26/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "DatabaseProcessor.h"

DatabaseProcessor::DatabaseProcessor(const char* s_directory, const char* s_folder,const char* s_main_directory)
{
    samp_directory = s_directory;
    samp_directory.append(s_folder);
    main_directory = s_main_directory;
    folder = s_folder;
    FeExtract.SetWriteFilePath(s_main_directory);
    GetFileNames();
}

DatabaseProcessor::~DatabaseProcessor()
{
    
}

int DatabaseProcessor::FeExSetup(const int set_analysis_size,const int set_hop_overlap,const int n_mel_bands,const int sample_rate,const double min_freq,const double max_freq,const int samp_chunk_size,const int mfcc_grouping)
{
    bool s_cmplt = FeExtract.Setup(set_analysis_size,set_hop_overlap,n_mel_bands,sample_rate,min_freq,max_freq,samp_chunk_size,mfcc_grouping);
    if(!s_cmplt) return 0;
    return 1;
}

int DatabaseProcessor::GetFileNames()
{
    DIR* direc = opendir(samp_directory.c_str());
    dirent* pdir;
    printf("\nFiles found in selected directory\n");
    while((pdir = readdir(direc))){
        std::string tmp_f_name = pdir->d_name;
        if(tmp_f_name.length() > 4){
            std::string file_ext = tmp_f_name.substr(tmp_f_name.length()-4,4);
            if(file_ext == ".wav"){
                printf("File : %s\n",pdir->d_name);
                file_names.push_back(tmp_f_name);
            }
        }
    }
    printf("\n");
    closedir(direc);
    if(file_names.size() == 0){
        printf("0 files in selected folder.\n");
        return 0;
    }
    return 1;
}

void DatabaseProcessor::Process(WRITER type)
{
    // Process all files in the folder and write to the MFCC/FFT database.
    if(file_names.size() > 0){
        for(int i = 0 ; i < file_names.size(); i++){
            std::string file_path = samp_directory;
            file_path.append(file_names[i]);
            FeExtract.ProcessSingleFile(file_path.c_str(),MARKOVANALYSIS);
        }
        FeExtract.WriteDatabase(file_names.size(),type);
    } else {
        printf("No files found to process.\n");
    }
}
