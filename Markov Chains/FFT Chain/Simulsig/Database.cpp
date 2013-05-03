//
//  DBReader.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 09/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "Database.h"

Database::Database()
{
    mfcc_total_num = 100;
    mfcc_num_coef = 40;
    database_frame_num = 10;
    mfcc_frame_r_length = mfcc_num_coef;
    mfcc_matrix = NULL;
    track_frame_length = NULL;
    fft_matrix = NULL;
    fft_matrix_two = NULL;
    grouping = 1;
}

Database::~Database()
{
    if(mfcc_matrix){
        printf("Freeing mfcc_matrix\n");
        for(int i = 0 ; i < gr_size; i++){
            free(mfcc_matrix[i]);
        }
        free(mfcc_matrix);
    }
    if(fft_matrix){
        printf("Freeing fft_matrix\n");
        for(int i = 0 ; i < gr_size; i++){
            free(fft_matrix[i]);
        }
       free(fft_matrix);
    }
    if(track_frame_length){
        free(track_frame_length);
    }
}

int Database::ReadFile(int file_num)
{
    printf("*--------------Database Read-----------------*\n");
    std::string dir = MACHINEPATH;
    dir.append("Database/Corpus/");
    std::string suffix = ".txt";
    char buf[10]; 
    sprintf(buf, "%d", file_num);
    
    std::ifstream database_setup;
    printf("Reading database setup file...\n");
    std::string db_file(dir);
    db_file.append("database_summary_");
    db_file.append(buf);
    db_file.append(suffix);
    database_setup.open(db_file.c_str());
    database_setup >> mfcc_num_coef;
    database_setup >> fft_frame_length;
    fft_frame_length += 2;
    database_setup >> database_frame_num;
    database_setup >> num_tracks;
    database_setup.close();
    printf("Complete\n");

    std::ifstream mfccFile;
    printf("Reading mfcc file...\n");
    std::string mf_file(dir);
    mf_file.append("mfcc_frames_");
    mf_file.append(buf);
    mf_file.append(suffix);
    mfccFile.open(mf_file);
    mfcc_matrix = (double**)malloc(sizeof(double*)*database_frame_num);
    for(int i = 0; i < database_frame_num; i++){
        mfcc_matrix[i] = (double*)malloc(sizeof(double)*mfcc_num_coef);
        for(int j = 0 ; j < mfcc_num_coef;j++){
            mfccFile >> mfcc_matrix[i][j];
        }
    }
    mfccFile.close();
    printf("Complete\n");
    
    std::ifstream fftFile;
    printf("Reading fft frames file...\n");
    std::string ft_file(dir);
    ft_file.append("fft_frames_");
    ft_file.append(buf);
    ft_file.append(suffix);
    fftFile.open(ft_file);
    fft_matrix = (double**)malloc(sizeof(double*)*database_frame_num);
    for(int i = 0 ; i < database_frame_num; i++){
        fft_matrix[i] = (double*)malloc(sizeof(double)*fft_frame_length);
        for(int j = 0; j < fft_frame_length; j++){
            fftFile >> fft_matrix[i][j];
        }
    }

    fftFile.close();
    printf("Complete\n");
    
    debug_print(("Reading track frame lengths..."));
    std::string track_frame_l_file(MACHINEPATH);
    track_frame_l_file.append("Database/Corpus/track_frame_lengths.txt");
    std::ifstream tr_frm_l_file;
    tr_frm_l_file.open(track_frame_l_file);
    track_frame_length = (int*)malloc(sizeof(int)*num_tracks);
    for(int i = 0 ; i < num_tracks; i++){
        tr_frm_l_file >> track_frame_length[i];
    }
    tr_frm_l_file.close();
    debug_print(("\tComplete\n"));
    
    debug_print(("Reading buf lengths..."));
    std::string buf_s_file_p(MACHINEPATH);
    buf_s_file_p.append("Database/Corpus/buffer_sizes.txt");
    std::ifstream buf_s_file;
    buf_s_file.open(buf_s_file_p);
    track_buffer_sizes = (long*)malloc(sizeof(long)*num_tracks);
    for(int i = 0 ; i < num_tracks; i++){
        buf_s_file >> track_buffer_sizes[i];
    }
    buf_s_file.close();
    debug_print(("\tComplete\n"));
    debug_print(("Number MFCC Coefficients     : %d\n",mfcc_num_coef));
    debug_print(("Number values p FFT frame    : %d\n",fft_frame_length));
    debug_print(("Total num frames in database : %ld\n\n",database_frame_num));
    printf("*--------------------------------------------*\n");

    gr_size = (int)database_frame_num;
    gr_coef = mfcc_num_coef;
    gr_fft = fft_frame_length;
    return 0;
}


int Database::ReadFileGroup(int file_num,int s_grouping,int type)
{
    grouping = s_grouping;
    
    // grouping must be > 0
    printf("\n*--------------Database Group Read-----------------*\n");

    std::string dir = MACHINEPATH;
    switch(type){
        case 0: dir.append("Database/Corpus/"); break;
        case 1: dir.append("Database/Corpus/DB/"); break;
    }
    //dir.append("Database/Corpus/");
    std::string suffix = ".txt";
    char buf[10]; // enough to hold all numbers up to 64-bits
    sprintf(buf, "%d", file_num);
    
    std::ifstream database_setup;
    printf("Reading database setup file...");
    std::string db_file(dir);
    db_file.append("database_summary_");
    db_file.append(buf);
    db_file.append(suffix);
    database_setup.open(db_file.c_str());
    database_setup >> mfcc_num_coef;
    database_setup >> fft_frame_length;
    fft_frame_length += 2;
    database_setup >> database_frame_num;
    database_setup >> num_tracks;
    database_setup.close();
    printf("\tComplete\n");
    
    std::ifstream mfccFile;
    printf("Reading mfcc file...");
    std::string mf_file(dir);
    mf_file.append("mfcc_frames_");
    mf_file.append(buf);
    mf_file.append(suffix);
    mfccFile.open(mf_file);
    
    grouped_size = (long)floor(database_frame_num/grouping);
    int grouped_coef = mfcc_num_coef * grouping;    
    gr_size = (int)grouped_size;
    gr_coef = grouped_coef;

    mfcc_matrix = (double**)malloc(sizeof(double*)*grouped_size);
    for(int i = 0; i < grouped_size; i++){
        mfcc_matrix[i] = (double*)malloc(sizeof(double)*grouped_coef);
    }
    for(int i = 0; i < grouped_size  ; i++){
        for(int j = 0 ; j < grouped_coef;j++){
            mfccFile >> mfcc_matrix[i][j];
        }
    }
    mfccFile.close();
    printf("\tComplete\n");
    
    std::ifstream fftFile;
    printf("Reading fft frames file...");
    std::string ft_file(dir);
    ft_file.append("fft_frames_");
    ft_file.append(buf);
    ft_file.append(suffix);
    fftFile.open(ft_file);
    fft_length_grouped = fft_frame_length * grouping;
    printf("fft len group : %d\n",fft_length_grouped);
    fft_matrix = (double**)malloc(sizeof(double*)*grouped_size);
    for(int i = 0 ; i < grouped_size; i++){
        fft_matrix[i] = (double*)malloc(sizeof(double)*fft_length_grouped);
    }
    for(int i = 0 ; i < grouped_size; i++){
        for(int j = 0; j < fft_length_grouped; j++){
            fftFile >> fft_matrix[i][j];
           // printf("%f ",fft_matrix[i][j]);
        }
    }
    fftFile.close();
    printf("\tComplete\n");
    
    gr_size = (int)grouped_size;
    gr_fft = fft_length_grouped;
    debug_print(("Reading track frame lengths..."));
    std::string track_frame_l_file(dir);
    track_frame_l_file.append("track_frame_lengths.txt");
    std::ifstream tr_frm_l_file;
    tr_frm_l_file.open(track_frame_l_file);
    track_frame_length = (int*)malloc(sizeof(int)*num_tracks);
    for(int i = 0 ; i < num_tracks; i++){
        tr_frm_l_file >> track_frame_length[i];
        printf("Track Length : %d \n",track_frame_length[i]);
    }
    tr_frm_l_file.close();
    
    debug_print(("Reading buf lengths..."));
    std::string buf_s_file_p(dir);
    buf_s_file_p.append("buffer_sizes.txt");
    std::ifstream buf_s_file;
    buf_s_file.open(buf_s_file_p);
    track_buffer_sizes = (long*)malloc(sizeof(long)*num_tracks);
    for(int i = 0 ; i < num_tracks; i++){
        buf_s_file >> track_buffer_sizes[i];
    }
    buf_s_file.close();
    debug_print(("\tComplete\n"));

    
    debug_print(("\tComplete\n"));
    printf("Num Tracks : %d\n",num_tracks);
    debug_print(("Number MFCC Coefficients     : %d\n",mfcc_num_coef));
    debug_print(("Number values p FFT frame    : %d\n",fft_frame_length));
    debug_print(("Total num frames in database : %ld\n",database_frame_num));
    debug_print(("Grouped MFCC num coefs       : %d\n",grouped_coef));
    debug_print(("Num frames when grouped      : %ld\n",grouped_size));
    GroupTrackLengths();
    printf("*--------------------------------------------------*\n");

    return 0;
}

int Database::GetFFTLengthGrouped(){
    return fft_length_grouped;
}

int Database::GetGroupedS()
{
    return (int)grouped_size;
}

void Database::SetFileDirectory(const char *dir)
{
    // do some error checks on the incoming dir
    //directory = dir;
}

int Database::GetFFTLength()
{
    return fft_frame_length;
}

long Database::GetNumFrames()
{
    return database_frame_num;
}

int Database::GetNumMFCCCoefs()
{
   return mfcc_num_coef;
}

double** Database::GetMFCCMatrix()
{
    return mfcc_matrix;
}

double** Database::GetFFTMatrix()
{
    return fft_matrix;
}

int Database::GetNumTracks()
{
    return num_tracks;
}

int Database::GetTrackFrameLengths(int index)
{
    return track_frame_length[index];
}

void Database::GroupTrackLengths()
{
    for(int i = 0 ; i < num_tracks; i++){
        track_frame_length[i] /= grouping;
    }
}

long Database::GetBufSizes(int idx)
{
    return track_buffer_sizes[idx];
}

int Database::GetLongestTrackLength()
{
    int largest = GetTrackFrameLengths(0);
    for(int i = 1 ; i < num_tracks; i++){
        if(GetTrackFrameLengths(i)  > largest){
            largest = GetTrackFrameLengths(i);
        }
    }
    return largest;
}