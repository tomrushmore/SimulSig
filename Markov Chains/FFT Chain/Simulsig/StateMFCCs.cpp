//
//  StateMFCCs.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 12/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "StateMFCCs.h"

StateMFCCs::StateMFCCs()
{
    num_active_clusters = 0;
    mfcc_at_state = NULL;
    fft_at_state = NULL;
    grouping = 0;
}

StateMFCCs::~StateMFCCs()
{
    if(mfcc_at_state){
        free(mfcc_at_state);
    }
    if(fft_at_state){
        free(fft_at_state);
    }
}

void StateMFCCs::SetupWrite(int s_num_states, int s_mfcc_num_coef,int s_fft_length,int s_num_tracks)
{
    num_states = s_num_states;
    mfcc_num_coef = s_mfcc_num_coef;
    fft_length = s_fft_length;
    num_tracks = s_num_tracks;
    mfcc_at_state = (double**)malloc(sizeof(double*)*num_states);
    for(int i = 0 ; i < num_states; i++){
        mfcc_at_state[i] = (double*)malloc(sizeof(double)*mfcc_num_coef);
    }
    fft_at_state = (double**)malloc(sizeof(double*)*num_states);
    for(int i = 0 ; i < num_states; i++){
        fft_at_state[i] = (double*)malloc(sizeof(double)*fft_length*MFCCGROUPING);
    }
    debug_print(("Markov Chain Initialization\n"));
    debug_print(("Num states : %d \n",num_states));
    debug_print(("Grouped FFT Length   : %d \n",fft_length*MFCCGROUPING));
    debug_print(("Grouped MFCC Length  : %d \n",mfcc_num_coef));

}

void StateMFCCs::SetMFCCFFTPairs(double **mfccs, double **fft_matrix,int *best_fit_indices)
{
    for(int i = 0; i < num_states; i++){
        int selected_pair = best_fit_indices[i];
        printf("Best MFCC Fit : %d \n",selected_pair);
        mfcc_at_state[i] = mfccs[selected_pair];
        fft_at_state[i] = fft_matrix[selected_pair];
        
    }
}

void StateMFCCs::WriteMFCCFFTPairs()
{
    std::string dir = MACHINEPATH;
    dir.append("Database/MFCC FFT States/");
    std::string fft_path(dir);
    printf("Writing clusters into : %s \n",dir.c_str());
    fft_path.append("fft_state_frames.txt");
    std::string mfcc_path(dir);
    mfcc_path.append("mfcc_state_frames.txt");
    std::string database_summary(dir);
    database_summary.append("state_summary.txt");
    
    std::ofstream database_sum_file(database_summary.c_str());
    printf("Writing state summary...");
    database_sum_file << num_states << " ";
    database_sum_file << mfcc_num_coef << " ";
    database_sum_file << fft_length << " ";
    database_sum_file << num_tracks << " ";
    database_sum_file << MFCCGROUPING << " ";
    database_sum_file << fft_length - HOPOVRLAPMS << " ";
    database_sum_file.close();
    printf("\tComplete\n");
    
    std::ofstream fft_file(fft_path.c_str());
    fft_file << std::fixed << std::setprecision(10);
    printf("Writing FFT state frames...");
    
    for(int i = 0; i < num_states; i++){
        for(int j = 0; j < fft_length*MFCCGROUPING; j++){
            fft_file << fft_at_state[i][j] << " ";
        }
    }
    fft_file.close();
    printf("\tComplete\n");
    
    std::ofstream mfcc_file(mfcc_path.c_str());
    mfcc_file << std::fixed << std::setprecision(10);
    printf("Writing MFFC state frames...");
    for(int i = 0 ; i < num_states; i++){
        for(int j = 0; j < mfcc_num_coef; j++){
            mfcc_file << mfcc_at_state[i][j] << " ";
        }
    }
    mfcc_file.close();
    printf("\tComplete\n\n");
}

void StateMFCCs::ReadMFCCFFTStatePairs()
{
    printf("*------------State Pairs---------------*\n");

    debug_print(("Reading State Summary file..."));
    std::string state_file(MACHINEPATH);
    state_file.append("Database/MFCC FFT States/state_summary.txt");
    std::ifstream state_summary;
    state_summary.open(state_file);
    state_summary >> num_active_clusters;
    state_summary >> mfcc_num_coef;
    state_summary >> fft_length;
    state_summary >> num_tracks;
    state_summary >> grouping;
    state_summary >> hop_overlap;
    state_summary.close();
    debug_print(("\tComplete\n"));
    debug_print(("Number of active clusters   : %d\n",num_active_clusters));
    debug_print(("Number of MFCC coef p group : %d\n",mfcc_num_coef));
    debug_print(("FFT Length                  : %d\n",fft_length));
    debug_print(("Num Tracks                  : %d\n",num_tracks));
    debug_print(("Grouping                    : %d\n",grouping));
    debug_print(("Hop Overlap Ms              : %d\n",hop_overlap));    
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
    
    debug_print(("Reading clustered MFCC Vectors..."));
    std::string mfcc_file(MACHINEPATH);
    mfcc_file.append("Database/MFCC FFT States/mfcc_state_frames.txt");
    std::ifstream mfcc_reader;
    mfcc_reader.open(mfcc_file);
    mfcc_at_state = (double**)malloc(sizeof(double*)*num_active_clusters);
    for(int i = 0 ; i < num_active_clusters; i++){
        mfcc_at_state[i] = (double*)malloc(sizeof(double)*mfcc_num_coef);
        for(int j = 0 ; j < mfcc_num_coef; j++){
            mfcc_reader >> mfcc_at_state[i][j];
        }
    }
    mfcc_reader.close();
    debug_print(("\tComplete\n"));
    
    debug_print(("Reading clustered FFT Vectors..."));
    std::string fft_file(MACHINEPATH);
    fft_file.append("Database/MFCC FFT States/fft_state_frames.txt");
    std::ifstream fft_reader;
    fft_reader.open(fft_file);
    fft_at_state = (double**)malloc(sizeof(double*)*num_active_clusters);
    //printf("CHOCK : %d\n",fft_length*grouping);
    for(int i = 0 ; i < num_active_clusters; i++){
        fft_at_state[i] = (double*)malloc(sizeof(double)*fft_length*grouping);
        for(int j = 0 ; j < fft_length*grouping; j++){
            fft_reader >> fft_at_state[i][j];
        }
    }
    fft_reader.close();
    debug_print(("\tComplete\n"));
    printf("*--------------------------------------*\n");

}

int StateMFCCs::GetLongestTrackLength()
{
    int largest = track_frame_length[0];
    for(int i = 1 ; i < num_tracks; i++){
        if(track_frame_length[i]  > largest){
            largest = track_frame_length[i];
        }
    }
    for(int i = 0 ; i < num_tracks; i++){
        // printf("LENGTH : %d \n",GetTrackFrameLengths(i));
    }
    return largest;
}

double** StateMFCCs::GetMFCCMatrix()
{
    if(mfcc_at_state){
        return mfcc_at_state;
    } else {
        return 0;
    }
//return mfcc_at_state;
}

double** StateMFCCs::GetFFTMatrix()
{
    if(fft_at_state){
        return fft_at_state;
    } else {
        return 0;
    }
}

int StateMFCCs::GetNumMFCCCoef()
{
    return mfcc_num_coef;
}

int StateMFCCs::GetNumActiveClusters()
{
    return num_active_clusters;
}

int StateMFCCs::GetFFTLength()
{
    return fft_length;
}

int StateMFCCs::GetNumTracks()
{
    return num_tracks;
}

int StateMFCCs::GetGrouping()
{
    return grouping;
}

int StateMFCCs::GetHopOverlap()
{
    return hop_overlap;
}
