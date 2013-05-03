//
//  FileWriteR.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 15/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "FileWriteR.h"

FileWriteR::FileWriteR()
{
    //directory = "/Users/thomasrushmore/Documents/Goldsmiths/Third Year/Final Project/Graphs";
    iter = 0;
}

FileWriteR::~FileWriteR()
{
    
}

void FileWriteR::SetFilePath(const char *path)
{
    directory = const_cast<char *>(path);
}

void FileWriteR::DBWriteSetup(const char *path)
{
    directory = const_cast<char *>(path);
}

void FileWriteR::DBWriteFFTMFCCSplit(const std::vector<double> *fft,const std::vector<double> *mfcc,const std::vector<long>* frames_length_vec,const int mfcc_num_coef,const int fft_vals_p_frame,const long int num_frames,const int file_num,const long num_tracks,const WRITER type,const std::vector<long>* buf_sizes)
{
    printf("*------------File Writer-------------*\n");
    
    std::string dir = directory;
    std::string rtype;
    switch(type){
        case NWRITEFILE:  rtype = "/Database/Corpus/DB/"; break;
        case WRITEFILE: rtype = "/Database/Corpus/";    break;
    }
    dir.append(rtype);
    std::string suffix = ".txt";
    char buf[10];
    sprintf(buf, "%d", file_num);
    
    std::string database_summary(dir);
    database_summary.append("database_summary_");
    database_summary.append(buf);
    database_summary.append(suffix);
    std::ofstream database_sum_file(database_summary.c_str());
    printf("Writing database summary...\n");
    database_sum_file << mfcc_num_coef << " ";
    database_sum_file << fft_vals_p_frame << " ";
    database_sum_file << num_frames  << " ";
    database_sum_file << num_tracks << " ";
    database_sum_file.close();
    printf("Database summary write complete.\n");
    
    std::string fft_path(dir);
    fft_path.append("fft_frames_");
    fft_path.append(buf);
    fft_path.append(suffix);
    std::ofstream fft_file(fft_path.c_str());
    fft_file << std::fixed << std::setprecision(10);
    printf("Writing FFT database...\n");
    printf("Write FFT SIZE : %ld\n",fft->size());
    for(int i = 0; i < fft->size(); i++){
        fft_file << (*fft)[i] << " ";
    }
    fft_file.close();
    printf("FFT database write complete.\n");
    
    std::string mfcc_path(dir);
    mfcc_path.append("mfcc_frames_");
    mfcc_path.append(buf);
    mfcc_path.append(suffix);
    std::ofstream mfcc_file(mfcc_path.c_str());
    mfcc_file << std::fixed << std::setprecision(10);
    printf("Writing MFFC database\n");

    for(int i = 0 ; i < mfcc->size(); i++){
        mfcc_file << (*mfcc)[i] << " ";
    }
    mfcc_file.close();
    printf("MFCC database write complete.\n");
    
    std::string frame_lengths(dir);
    frame_lengths.append("track_frame_lengths");
    frame_lengths.append(suffix);
    std::ofstream frame_lengths_file(frame_lengths.c_str());
    printf("Writing track frame lengths database...\n");
    for(int i = 0 ; i < frames_length_vec->size(); i++){
        frame_lengths_file << (*frames_length_vec)[i] << " ";
    }
    frame_lengths_file.close();
    printf("Track frame lengths write complete.\n");
    
    std::string buf_lengths(dir);
    buf_lengths.append("buffer_sizes");
    buf_lengths.append(suffix);
    std::ofstream buf_size_file(buf_lengths.c_str());
    for(int i = 0 ; i < buf_sizes->size(); i++){
        buf_size_file << (*buf_sizes)[i] << " ";
    }
    printf("*------------------------------------*\n");
    
}

void FileWriteR::DBWriteSingle(const std::vector<double> *frames)
{
    std::string path(directory);
    path.append("/FFTs/fftwritetest.txt");
    std::ofstream out_file(path);
    out_file << std::fixed << std::setprecision(10);
    printf("Writing File...\n");
    for(int i = 0 ; i < frames->size();i++){
        out_file << (*frames)[i] << " ";
    }
    debug_print(("Number of double values : %ld\n",frames->size()));
}

void FileWriteR::WriteVector(const double *vector,const int length)
{
    std::string dir(directory);
    dir.append("example");
    char buf[10]; // enough to hold all numbers up to 64-bits
    sprintf(buf, "%d", iter);
    dir.append(buf);
    dir.append(".txt");
    std::cout << dir << std::endl;
    std::ofstream myfile;
    myfile.open (dir.c_str());
    for(int i = 0 ; i < length; i++){
        myfile << vector[i] << " ";
    }
    myfile.close();
    iter++;
}

std::string FileWriteR::WriteVector(const double *vector,const u_short_int *length,const char * name)
{
    std::string dir(directory);
    dir.append("Database/");
    dir.append(name);
    dir.append(".txt");
    std::ofstream myfile;
    myfile.open (dir.c_str());
    for(int i = 0 ; i < *length; i++){
        myfile << vector[i] << " ";
    }
    myfile.close();
    iter++;
    return dir;
}

void FileWriteR::DBWriteFreqGraph(std::vector<double>* freqs,const int type)
{
    std::string full_path = "/Users/thomasrushmore/Documents/Goldsmiths/Third Year/Final Project/Graphs/";

    type ? full_path.append("mel_spacing") : full_path.append("linear_spacing");
    full_path.append("_graph");
    full_path.append(".txt");
    std::ofstream out_file(full_path);
    out_file << std::fixed << std::setprecision(10);
    printf("Writing File...\n");
    for(int i = 0 ; i < freqs->size();i++){
        out_file << (*freqs)[i] << " ";
    }
    debug_print(("Number of double values : %ld\n",freqs->size()));
}

void FileWriteR::WriteMelBanks(const double **filterbank,const int bank_length,const int num_bands)
{
    std::string path(directory);
    path.append("/mels/");
    std::string txt = ".txt";
    int fn = 0;
    for(int i = 0 ; i < num_bands; i++){
        char buf[10];
        std::string paf = path;
        sprintf(buf, "%d", fn);
        paf.append(buf);
        paf.append(txt);
        std::ofstream out_file(paf);
        out_file << std::fixed << std::setprecision(10);
        fn++;
        for(int j = 0 ; j < bank_length; j++){
            out_file << filterbank[i][j] << " " ;
        }
        out_file.close();
    }
}

void FileWriteR::DBWriteMarkov(const std::vector<double> *states)
{
    std::string dir = MACHINEPATH;
    dir.append("/Database/Corpus/");
    std::string suffix = ".txt";
    
    std::string markov(dir);
    markov.append("markovs");
    markov.append(suffix);
    std::ofstream markov_file(markov.c_str());
    printf("Writing markov summary...\n");
    printf("Markov write size : %ld \n",states->size());
    for(int i = 0 ; i < states->size(); i++){
        markov_file << (*states)[i] << " ";
    }
    markov_file.close();
    printf("Markov_file write complete.\n");

}
